#include "FastAdamsKernel.h"

#include "FastCudaSolver.h"
#include "Meteorites.Core/Constants.h"
#include "Meteorites.Core/Solvers/Adams.h"

namespace
{

#ifdef DISPLAY_WARP_DIVERGENCE
__device__ inline void PrintWarpMask(const char *message)
{
  unsigned int mask = __activemask();
  if(threadIdx.x == 0 && blockIdx.x == 0)
  {
    printf(message);
    for (int i = 31; i >= 0; i--)
    {
      printf("%d", (mask >> i) & 1);
      if (i % 8 == 0 && i != 0) printf(" ");
    }
    printf("\n");
  }
}
#endif

__global__ void RestoreMeteoroidsKernel_(VirtualMeteoroid *meteoroids, const curandState *states,
                                         const size_t size, const real v0, const real h0)
{
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid < size)
  {
    curandState state = states[tid];
        
    meteoroids[tid] = VirtualMeteoroid(1e5 + curand_uniform(&state) * (5e6 - 1e5),
                                       0.1 + curand_uniform(&state) * (0.9 - 0.1),
                                       2000.0 + curand_uniform(&state) * (5000.0 - 2000.0),
                                       0.5 + curand_uniform(&state) * (2.5 - 0.5),
                                       0.0 + curand_uniform(&state) * (0.25 - 0.0),
                                       10.0 + curand_uniform(&state) * (500.0 - 10.0),
                                       v0, h0,
                                       curand_uniform(&state) * (M_PI / 2));
  }
}

template <unsigned int STEPS>
__device__ void InitContext(Adams::Layer *steps, const VirtualMeteoroid &meteoroid, real dt)
{
  Adams::Unchangeable params(meteoroid);
  Adams::SetLayer(steps[STEPS], params,
                  meteoroid.V0, meteoroid.Gamma0, meteoroid.h0, meteoroid.l0, meteoroid.M0);

  Adams::OneStepIteration(steps[STEPS - 1], steps[STEPS], params, dt);

  if constexpr (STEPS >= 2)
  {
    Adams::TwoStepIteration(steps[STEPS - 2], steps[STEPS - 1],
                            steps[STEPS], params, dt);
  }

  if constexpr (STEPS >= 3)
  {
    Adams::ThreeStepIteration(steps[STEPS - 3], steps[STEPS - 2],
                              steps[STEPS - 1], steps[STEPS], params, dt);
  }
}

template <unsigned int STEPS>
__device__ inline bool
AdamsStep(const real *timestamps, const size_t n_timestamps, const real dt, const real timeout,
          Adams::Layer *steps, real *V_arg, real *h_arg, size_t &timestamp,
          size_t &nxt, const Adams::Unchangeable &params, real &t)
{
  if (timestamp < n_timestamps && t >= timestamps[timestamp])
  {
    const auto &step = steps[(nxt + 1) % (STEPS + 1)];
    V_arg[timestamp] = step.V;
    h_arg[timestamp] = step.h;
    timestamp++;
  }

  Adams::Iteration<STEPS>(steps, params, nxt, dt);
  auto M = steps[nxt].M;
  auto h = steps[nxt].h;

  t += dt;

  nxt = (nxt + STEPS) % (STEPS + 1);

  if (M <= (real)0.01 || h <= (real)0.0 || t >= timeout)
  {
    return false;
  }
  return true;
}

__device__ inline real
GPUL2Compute(const real *v, const real *h, const real *ref_v, const real *ref_h, const size_t size)
{
  real v_sum = 0.0, h_sum = 0.0;

  for (int i = 0; i < size; i++)
  {
    real dv = (ref_v[i] - v[i]) / ref_v[0];
    v_sum += dv * dv;

    real dh = (ref_h[i] - h[i]) / ref_h[0];
    h_sum += dh * dh;
  }
  return (sqrt(v_sum) + sqrt(h_sum)) / sqrt(size);
}

template<size_t BEST_METEOROIDS_BUFFER_SIZE>
__device__ void
InsertToTopSmallest(real *devs, curandState *states,
                    const real dev, const curandState state,
                    const real border_dev = std::numeric_limits<real>::max())
{
  if (dev > devs[BEST_METEOROIDS_BUFFER_SIZE - 1] || dev > border_dev)
  { return; }

  int pos = 0;
  while (pos < BEST_METEOROIDS_BUFFER_SIZE && dev > devs[pos])
  { pos++; }

  if (pos < BEST_METEOROIDS_BUFFER_SIZE)
  {
    for (int i = BEST_METEOROIDS_BUFFER_SIZE - 1; i > pos; i--)
    {
      devs[i] = devs[i-1];
      states[i] = states[i-1];
    }
    devs[pos] = dev;
    states[pos] = state;
  }
}

__device__ inline void
GenerateCase(VirtualMeteoroid &curr_meteoroid, curandState &curr_state, const real v0, const real h0)
{
  curr_meteoroid.H = 1e5 + curand_uniform(&curr_state) * (5e6 - 1e5);
  curr_meteoroid.Ch = 0.1 + curand_uniform(&curr_state) * (0.9 - 0.1);
  curr_meteoroid.Rho = 2000.0 + curand_uniform(&curr_state) * (5000.0 - 2000.0);
  curr_meteoroid.Cd = 0.5 + curand_uniform(&curr_state) * (2.5 - 0.5);
  curr_meteoroid.Cl = 0.0 + curand_uniform(&curr_state) * (0.25 - 0.0);
  curr_meteoroid.M0 = 10.0 + curand_uniform(&curr_state) * (500.0 - 10.0);
  curr_meteoroid.V0 = v0;
  curr_meteoroid.h0 = h0;
  curr_meteoroid.Gamma0 = curand_uniform(&curr_state) * (M_PI / 2);
}

//  WARP DIVERGENCE SCHEME:

//  ADAMS_STEPS  ||||||||
//               ||||||||
//               ||  ||||
//               |   ||
//                    |
//  UPDATE_CASE  ||||||||
//               ||||||||
//  ADAMS_STEPS  ||||||||

template <unsigned int STEPS, size_t BEST_METEOROIDS_BUFFER_SIZE>
__global__ void
FastAdamsKernel(const uint64_t *seeds, const size_t n_problems,
                const real dt, const real timeout,
                const real *timestamps, const size_t n_timestamps,
                real *functional_args, const real *functional_args_references,
                real *best_deviations, curandState *best_curand_states,
                const real border_deviation,
                const size_t meteoroids_per_thread)
{
  const size_t tid = blockIdx.x * blockDim.x + threadIdx.x;

  if (tid + 1 > n_problems / meteoroids_per_thread)
  { return; }

  curandState base_state;
  curand_init(seeds[tid], 0, 0, &base_state);
  curandState curr_state = base_state;
  VirtualMeteoroid curr_meteoroid;
  
  const real *ref_v_args = functional_args_references,
             *ref_h_args = functional_args_references + n_timestamps;
  real       *v_args = functional_args + tid * n_timestamps,
             *h_args = functional_args + tid * n_timestamps + gridDim.x * blockDim.x * n_timestamps;

  curandState local_best_curand_states[BEST_METEOROIDS_BUFFER_SIZE];
  real        local_best_deviations[BEST_METEOROIDS_BUFFER_SIZE];
  for (int i = 0; i < BEST_METEOROIDS_BUFFER_SIZE; i++)
  { local_best_deviations[i] = std::numeric_limits<real>::max(); }

  Adams::Layer steps[STEPS + 1];
  real         t;
  size_t       timestamp;
  size_t       nxt;
  real         curr_deviation;

  for (size_t i = 0; i < meteoroids_per_thread; i++)
  {
    GenerateCase(curr_meteoroid, curr_state, ref_v_args[0], ref_h_args[0]);

    for (size_t i = 0; i < n_timestamps; i++)
    {
      v_args[i] = 0;
      h_args[i] = 0;
    }
    t = dt * STEPS;
    timestamp = 0;
    nxt = STEPS;
    InitContext<STEPS>(steps, curr_meteoroid, dt);

    while (true)
    {
#ifdef DISPLAY_WARP_DIVERGENCE
      PrintWarpMask("ADAMS_STEP: ");
#endif
      if (!AdamsStep<STEPS>(timestamps, n_timestamps, dt, timeout, steps, v_args, h_args,
                            timestamp, nxt, curr_meteoroid, t)) 
      { break; }
    }
#ifdef DISPLAY_WARP_DIVERGENCE
    PrintWarpMask("UPDATE: ");
#endif
    curr_deviation = GPUL2Compute(v_args, h_args, ref_v_args, ref_h_args, n_timestamps);
    InsertToTopSmallest<BEST_METEOROIDS_BUFFER_SIZE>
        (local_best_deviations, local_best_curand_states,
         curr_deviation, base_state, border_deviation);
    base_state = curr_state;
  }

  for (int i = 0; i < BEST_METEOROIDS_BUFFER_SIZE; i++)
  {
    best_deviations[tid * BEST_METEOROIDS_BUFFER_SIZE + i] = local_best_deviations[i];
    best_curand_states[tid * BEST_METEOROIDS_BUFFER_SIZE + i] = local_best_curand_states[i];
  }
}

//  WARP DIVERGENCE SCHEME:

//  ADAMS_STEP   ||||||||
//               ||||||||
//  UPDATE_CASE       |
//                    |
//  ADAMS_STEP   ||||||||
//               ||||||||
//               ||||||||
//  UPDATE_CASE    |   |
//                 |   |

template <unsigned int STEPS, size_t BEST_METEOROIDS_BUFFER_SIZE>
__global__ void
FastAdamsBalancedKernel(const uint64_t *seeds, const size_t n_problems,
                        const real dt, const real timeout,
                        const real *timestamps, const size_t n_timestamps,
                        real *functional_args, const real *functional_args_references,
                        real *best_deviations, curandState *best_curand_states,
                        const real border_deviation,
                        const size_t meteoroids_per_thread)
{
  const size_t tid = blockIdx.x * blockDim.x + threadIdx.x;

  if (tid + 1 > n_problems / meteoroids_per_thread)
  { return; }

  __shared__ uint32_t meteoroids_counter;
  if (threadIdx.x == 0)
  { meteoroids_counter = 0; }

  uint32_t meteoroids_per_block;
  if (blockIdx.x == n_problems / meteoroids_per_thread * blockDim.x)
  { meteoroids_per_block = n_problems % meteoroids_per_thread * blockDim.x; }
  else
  { meteoroids_per_block = meteoroids_per_thread * blockDim.x; }

  curandState base_state;
  curand_init(seeds[tid], 0, 0, &base_state);
  curandState curr_state = base_state;
  VirtualMeteoroid curr_meteoroid;
  
  const real *ref_v_args = functional_args_references,
             *ref_h_args = functional_args_references + n_timestamps;
  real       *v_args = functional_args + tid * n_timestamps,
             *h_args = functional_args + tid * n_timestamps + gridDim.x * blockDim.x * n_timestamps;

  curandState local_best_curand_states[BEST_METEOROIDS_BUFFER_SIZE];
  real        local_best_deviations[BEST_METEOROIDS_BUFFER_SIZE];
  for (int i = 0; i < BEST_METEOROIDS_BUFFER_SIZE; i++)
  { local_best_deviations[i] = std::numeric_limits<real>::max(); }

  Adams::Layer steps[STEPS + 1];
  real         t;
  size_t       timestamp;
  size_t       nxt;
  real         curr_deviation;

  GenerateCase(curr_meteoroid, curr_state, ref_v_args[0], ref_h_args[0]);

  for (size_t i = 0; i < n_timestamps; i++)
  {
    v_args[i] = 0;
    h_args[i] = 0;
  }
  t = dt * STEPS;
  timestamp = 0;
  nxt = STEPS;
  InitContext<STEPS>(steps, curr_meteoroid, dt);

  while (true)
  {
#ifdef DISPLAY_WARP_DIVERGENCE
    PrintWarpMask("ADAMS_STEP: ");
#endif
    if (AdamsStep<STEPS>(timestamps, n_timestamps, dt, timeout, steps, v_args, h_args,
                         timestamp, nxt, curr_meteoroid, t))
    { continue; }
    else
    {
#ifdef DISPLAY_WARP_DIVERGENCE
      PrintWarpMask("UPDATE: ");
#endif
      curr_deviation = GPUL2Compute(v_args, h_args, ref_v_args, ref_h_args, n_timestamps);
      InsertToTopSmallest<BEST_METEOROIDS_BUFFER_SIZE>
          (local_best_deviations, local_best_curand_states,
           curr_deviation, base_state, border_deviation);
      base_state = curr_state;

      if (atomicAdd(&meteoroids_counter, 1) >= meteoroids_per_block)
      { break; }

      GenerateCase(curr_meteoroid, curr_state, ref_v_args[0], ref_h_args[0]);

      for (size_t i = 0; i < n_timestamps; i++)
      {
        v_args[i] = 0;
        h_args[i] = 0;
      }
      t = dt * STEPS;
      timestamp = 0;
      nxt = STEPS;
      InitContext<STEPS>(steps, curr_meteoroid, dt);
    }
  }
  for (int i = 0; i < BEST_METEOROIDS_BUFFER_SIZE; i++)
  {
    best_deviations[tid * BEST_METEOROIDS_BUFFER_SIZE + i] = local_best_deviations[i];
    best_curand_states[tid * BEST_METEOROIDS_BUFFER_SIZE + i] = local_best_curand_states[i];
  }
}

} // unnamed namespace

template <unsigned int STEPS, size_t BEST_METEOROIDS_BUFFER_SIZE>
void FastAdamsKernel(const uint64_t *seeds, const size_t n_problems,
                     const real dt, const real timeout,
                     const real *timestamps, size_t n_timestamps,
                     real *functional_args, const real *functional_args_references,
                     real *best_deviations, curandState *best_curand_states, const real border_deviation,
                     const size_t meteoroids_per_thread,
                     const size_t blocks_num, const size_t threads_per_block)
{
  assert(n_problems > 0);
  assert(threads_per_block > 0);
  assert(blocks_num > 0);
  assert(n_problems % meteoroids_per_thread == 0);

  FastAdamsKernel<STEPS, BEST_METEOROIDS_BUFFER_SIZE>
      <<<blocks_num, threads_per_block>>>(seeds, n_problems,
                                          dt, timeout,
                                          timestamps, n_timestamps,
                                          functional_args, functional_args_references,
                                          best_deviations, best_curand_states, border_deviation,
                                          meteoroids_per_thread);
  HANDLE_ERROR(cudaGetLastError());
}

// template-specified functions compiles only this way (C++ moment)
template void FastAdamsKernel<1u, FastCudaSolverConfig::best_meteoroids_buffer_size>
    (const uint64_t *seeds, const size_t n_problems,
     const real dt, const real timeout,
     const real *timestamps, size_t n_timestamps,
     real *functional_args, const real *functional_args_references,
     real *best_deviations, curandState *best_curand_states, const real border_deviation,
     const size_t meteoroids_per_thread,
     const size_t blocks_num, const size_t threads_per_block);
template void FastAdamsKernel<2u, FastCudaSolverConfig::best_meteoroids_buffer_size>
    (const uint64_t *seeds, const size_t n_problems,
     const real dt, const real timeout,
     const real *timestamps, size_t n_timestamps,
     real *functional_args, const real *functional_args_references,
     real *best_deviations, curandState *best_curand_states, const real border_deviation,
     const size_t meteoroids_per_thread,
     const size_t blocks_num, const size_t threads_per_block);
template void FastAdamsKernel<3u, FastCudaSolverConfig::best_meteoroids_buffer_size>
    (const uint64_t *seeds, const size_t n_problems,
     const real dt, const real timeout,
     const real *timestamps, size_t n_timestamps,
     real *functional_args, const real *functional_args_references,
     real *best_deviations, curandState *best_curand_states, const real border_deviation,
     const size_t meteoroids_per_thread,
     const size_t blocks_num, const size_t threads_per_block);

void RestoreMeteoroidsKernel(VirtualMeteoroid *meteoroids, const curandState *states,
                             const size_t size, const real v0, const real h0)
{
  RestoreMeteoroidsKernel_<<<(size - 1) / 64 + 1, 64>>>(meteoroids, states, size, v0, h0);
  HANDLE_ERROR(cudaGetLastError());
}