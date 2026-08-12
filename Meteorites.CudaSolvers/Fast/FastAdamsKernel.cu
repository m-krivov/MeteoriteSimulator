#include "FastAdamsKernel.h"

#include "FastCudaSolver.h"
#include "Meteorites.Core/Constants.h"
#include "Meteorites.Core/Solvers/Adams.h"
#include "CudaMeteoroidsGenerator.h"

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



template <uint32_t STEPS>
__device__ void InitContext(Adams::Layer &curr_layer, Adams::Layer *steps,
                            const VirtualMeteoroid &meteoroid,
                            real dt)
{
  curr_layer = { meteoroid.V0, meteoroid.Gamma0, meteoroid.h0, meteoroid.l0, meteoroid.M0 };
  Adams::OneStepIteration(curr_layer, steps[0], meteoroid, dt);

  if constexpr (STEPS >= 2)
  {
    Adams::TwoStepIteration(curr_layer, steps[1], steps[0], meteoroid, dt);
  }

  if constexpr (STEPS >= 3)
  {
    Adams::ThreeStepIteration(curr_layer, steps[2], steps[1], steps[0], meteoroid, dt);
  }
}

template <uint32_t STEPS, typename... LAYERS>
__device__ inline bool
AdamsStep(const real *timestamps, const uint32_t n_timestamps,
          const real dt, const real timeout,
          TrajectoryPoint *points, uint32_t &timestamp,
          const VirtualMeteoroid &meteoroid, real &t,
          Adams::Layer &curr_layer, Adams::Layer &curr_step, const LAYERS&... prev_steps)
{
  static_assert(sizeof...(LAYERS) == STEPS - 1);
  static_assert((std::is_same_v<LAYERS, Adams::Layer> && ...));

  if (timestamp < n_timestamps && t >= timestamps[timestamp])
  {
    points[timestamp] = { curr_layer.V, curr_layer.h };
    timestamp++;
  }

  auto t_ = std::forward_as_tuple(prev_steps...);
  if constexpr (STEPS == 1)
  {
    Adams::OneStepIteration(curr_layer, curr_step, meteoroid, dt);
  }

  else if constexpr (STEPS == 2)
  {
    Adams::TwoStepIteration(curr_layer, curr_step, std::get<0>(t_), meteoroid, dt);
  }

  else if constexpr (STEPS == 3)
  {
    Adams::ThreeStepIteration(curr_layer, curr_step, std::get<0>(t_), std::get<1>(t_), meteoroid, dt);
  }

  t += dt;

  if (curr_layer.M <= (real)0.01 || curr_layer.h <= (real)0.0 || t >= timeout)
  {
    return false;
  }
  return true;
}

// Performs STEPS Adams steps to unroll the loop over the local steps[STEPS] array
template <uint32_t STEPS>
__device__ inline bool
AdamsCycle(const real *timestamps, const uint32_t n_timestamps,
           const real dt, const real timeout,
           Adams::Layer &curr_layer, Adams::Layer *steps,
           TrajectoryPoint *points, uint32_t &timestamp,
           const VirtualMeteoroid &meteoroid, real &t)
{
  if constexpr (STEPS == 1)
  {
    if (!AdamsStep<STEPS>(timestamps, n_timestamps, dt, timeout, points, timestamp, meteoroid, t,
                          curr_layer, steps[0]))
    { return false; }

    return true;
  }

  else if constexpr (STEPS == 2)
  {
    if (!AdamsStep<STEPS>(timestamps, n_timestamps, dt, timeout, points, timestamp, meteoroid, t,
                          curr_layer, steps[0], steps[1]))
    { return false; }

    if (!AdamsStep<STEPS>(timestamps, n_timestamps, dt, timeout, points, timestamp, meteoroid, t,
                          curr_layer, steps[1], steps[0]))
    { return false; }

    return true;
  }

  else if constexpr (STEPS == 3)
  {
    if (!AdamsStep<STEPS>(timestamps, n_timestamps, dt, timeout, points, timestamp, meteoroid, t,
                          curr_layer, steps[0], steps[2], steps[1]))
    { return false; }

    if (!AdamsStep<STEPS>(timestamps, n_timestamps, dt, timeout, points, timestamp, meteoroid, t,
                          curr_layer, steps[1], steps[0], steps[2]))
    { return false; }

    if (!AdamsStep<STEPS>(timestamps, n_timestamps, dt, timeout, points, timestamp, meteoroid, t,
                          curr_layer, steps[2], steps[1], steps[0]))
    { return false; }

    return true;
  }
}

__device__ inline real
ComputeL2(const TrajectoryPoint *points, const TrajectoryPoint *ref_points,
          const uint32_t n_timestamps)
{
  real v_sum = 0.0, h_sum = 0.0;

  for (uint32_t i = 0; i < n_timestamps; i++)
  {
    real dv = (ref_points[i].v - points[i].v) / ref_points[i].v;
    v_sum += dv * dv;

    real dh = (ref_points[i].h - points[i].h) / ref_points[i].h;
    h_sum += dh * dh;
  }
  return (sqrt(v_sum) + sqrt(h_sum)) / sqrt(n_timestamps);
}

template<uint32_t BESTS_PER_THREAD>
__device__ void
InsertToTopSmallest(MeteoroidDeviation *devs, const MeteoroidDeviation &dev,
                    const real threshold_dev = std::numeric_limits<real>::max())
{
  if (dev.dev > devs[BESTS_PER_THREAD - 1].dev || dev.dev > threshold_dev)
  { return; }

  uint32_t pos = 0;
  while (pos < BESTS_PER_THREAD && dev.dev > devs[pos].dev)
  { pos++; }

  if (pos < BESTS_PER_THREAD)
  {
    for (uint32_t i = BESTS_PER_THREAD - 1; i > pos; i--)
    {
      devs[i] = devs[i-1];
    }
    devs[pos] = dev;
  }
}



//  WARP DIVERGENCE SCHEME:

//  ADAMS_STEPS  ||||||||
//               ||||||||
//               ||  ||||
//               |   ||
//                    |
//  UPDATE       ||||||||
//               ||||||||
//  ADAMS_STEPS  ||||||||

template <uint32_t STEPS, uint32_t BESTS_PER_THREAD>
__global__ void
FastAdamsKernel(const uint64_t *seeds,
                const real *timestamps, const uint32_t n_timestamps,
                const TrajectoryPoint *ref_points, TrajectoryPoint *context_points,
                const real dt, const real timeout,
                MeteoroidDeviation *global_best_meteoroids, const real threshold_dev,
                const uint32_t meteoroids_per_thread)
{
  const uint32_t tid = blockIdx.x * blockDim.x + threadIdx.x;

  CudaMeteoroidsGenerator generator(seeds[tid]);

  MeteoroidDeviation local_best_meteoroids[BESTS_PER_THREAD];
  for (uint32_t i = 0; i < BESTS_PER_THREAD; i++)
  {
    local_best_meteoroids[i].dev = std::numeric_limits<real>::max();
  }

  // Meteoroid we are currently simulating
  VirtualMeteoroid  meteoroid;
  // Trajectory points for currently simulated meteoroid
  TrajectoryPoint *points = context_points + tid * n_timestamps;

  Adams::Layer curr_layer;
  Adams::Layer steps[STEPS];
  real         t;
  uint32_t     timestamp;

  // Main loop
  for (uint32_t i = 0; i < meteoroids_per_thread; i++)
  {
    // Init context
    {
    for (uint32_t i = 0; i < n_timestamps; i++) { points[i] = {0.0, 0.0}; }

    generator.Next(meteoroid, ref_points[0]);

    t = dt * STEPS;
    timestamp = 0;
    InitContext<STEPS>(curr_layer, steps, meteoroid, dt);
    }

    while (true)
    {
#ifdef DISPLAY_WARP_DIVERGENCE
      PrintWarpMask("ADAMS_STEP: ");
#endif
      if (!AdamsCycle<STEPS>(timestamps, n_timestamps, dt, timeout, curr_layer, steps,
                             points, timestamp, meteoroid, t))
      { break; }
    }
#ifdef DISPLAY_WARP_DIVERGENCE
    PrintWarpMask("UPDATE: ");
#endif
    real dev = ComputeL2(points, ref_points, n_timestamps);
    InsertToTopSmallest<BESTS_PER_THREAD>(local_best_meteoroids, { dev, generator.GetOffset() });
  }

  for (uint32_t i = 0; i < BESTS_PER_THREAD; i++)
  {
    global_best_meteoroids[tid * BESTS_PER_THREAD + i] = local_best_meteoroids[i];
  }
}

//  WARP DIVERGENCE SCHEME:

//  ADAMS_STEPS  ||||||||
//               ||||||||
//  UPDATE            |
//                    |
//  ADAMS_STEPS  ||||||||
//               ||||||||
//               ||||||||
//  UPDATE         |   |
//                 |   |

template <uint32_t STEPS, uint32_t BESTS_PER_THREAD>
__global__ void
FastAdamsBalancedKernel(const uint64_t *seeds,
                        const real *timestamps, const uint32_t n_timestamps,
                        const TrajectoryPoint *ref_points, TrajectoryPoint *context_points,
                        const real dt, const real timeout,
                        MeteoroidDeviation *global_best_meteoroids, const real threshold_dev,
                        const uint32_t meteoroids_per_thread)
{
  const uint32_t tid = blockIdx.x * blockDim.x + threadIdx.x;

  CudaMeteoroidsGenerator generator(seeds[tid]);

  MeteoroidDeviation local_best_meteoroids[BESTS_PER_THREAD];
  for (uint32_t i = 0; i < BESTS_PER_THREAD; i++)
  {
    local_best_meteoroids[i].dev = std::numeric_limits<real>::max();
  }

  // Meteoroid we are currently simulating
  VirtualMeteoroid  meteoroid;
  // Trajectory points for currently simulated meteoroid
  TrajectoryPoint *points = context_points + tid * n_timestamps;

  Adams::Layer curr_layer;
  Adams::Layer steps[STEPS];
  real         t;
  uint32_t     timestamp;

  __shared__ uint32_t meteoroids_counter;
  if (threadIdx.x == 0) meteoroids_counter = 0;

  // Init context
  {
  for (uint32_t i = 0; i < n_timestamps; i++) { points[i] = {0.0, 0.0}; }

  generator.Next(meteoroid, ref_points[0]);

  t = dt * STEPS;
  timestamp = 0;
  InitContext<STEPS>(curr_layer, steps, meteoroid, dt);
  }

  // Main loop
  while (true)
  {
#ifdef DISPLAY_WARP_DIVERGENCE
    PrintWarpMask("ADAMS_STEP: ");
#endif
    if (AdamsCycle<STEPS>(timestamps, n_timestamps, dt, timeout, curr_layer, steps,
                          points, timestamp, meteoroid, t))
    { continue; }
    else
    {
#ifdef DISPLAY_WARP_DIVERGENCE
      PrintWarpMask("UPDATE: ");
#endif
      real dev = ComputeL2(points, ref_points, n_timestamps);
      InsertToTopSmallest<BESTS_PER_THREAD>(local_best_meteoroids, { dev, generator.GetOffset() });

      if (atomicAdd(&meteoroids_counter, 1) >= meteoroids_per_thread * blockDim.x)
      { break; }

      // Init context
      {
      for (uint32_t i = 0; i < n_timestamps; i++) { points[i] = {0.0, 0.0}; }

      generator.Next(meteoroid, ref_points[0]);

      t = dt * STEPS;
      timestamp = 0;
      InitContext<STEPS>(curr_layer, steps, meteoroid, dt);
      }
    }
  }

  for (uint32_t i = 0; i < BESTS_PER_THREAD; i++)
  {
    global_best_meteoroids[tid * BESTS_PER_THREAD + i] = local_best_meteoroids[i];
  }
}

} // unnamed namespace

template <uint32_t STEPS, uint32_t BEST_METEOROIDS_PER_THREAD>
void FastAdamsKernel(const uint64_t *seeds,
                     const real *timestamps, const uint32_t n_timestamps,
                     const TrajectoryPoint *reference_points, TrajectoryPoint *context_points,
                     const real dt, const real timeout,
                     MeteoroidDeviation *global_best_meteoroids, const real threshold_dev,
                     const uint32_t meteoroids_per_thread,

                     const size_t blocks_num, const size_t threads_per_block)
{
  assert(threads_per_block > 0);
  assert(blocks_num > 0);

  FastAdamsKernel<STEPS, BEST_METEOROIDS_PER_THREAD><<<blocks_num, threads_per_block>>>
    (seeds, timestamps, n_timestamps, reference_points, context_points, dt, timeout,
     global_best_meteoroids, threshold_dev, meteoroids_per_thread);

  HANDLE_ERROR(cudaDeviceSynchronize());
}

// Explicit template instantiations are required (C++ moment)
template void FastAdamsKernel<1u, FastCudaSolverConfig::best_meteoroids_per_thread>
    (const uint64_t *seeds,
     const real *timestamps, const uint32_t n_timestamps,
     const TrajectoryPoint *reference_points, TrajectoryPoint *context_points,
     const real dt, const real timeout,
     MeteoroidDeviation *global_best_meteoroids, const real threshold_dev,
     const uint32_t meteoroids_per_thread,
     const size_t blocks_num, const size_t threads_per_block);
template void FastAdamsKernel<2u, FastCudaSolverConfig::best_meteoroids_per_thread>
    (const uint64_t *seeds,
     const real *timestamps, const uint32_t n_timestamps,
     const TrajectoryPoint *reference_points, TrajectoryPoint *context_points,
     const real dt, const real timeout,
     MeteoroidDeviation *global_best_meteoroids, const real threshold_dev,
     const uint32_t meteoroids_per_thread,
     const size_t blocks_num, const size_t threads_per_block);
template void FastAdamsKernel<3u, FastCudaSolverConfig::best_meteoroids_per_thread>
    (const uint64_t *seeds,
     const real *timestamps, const uint32_t n_timestamps,
     const TrajectoryPoint *reference_points, TrajectoryPoint *context_points,
     const real dt, const real timeout,
     MeteoroidDeviation *global_best_meteoroids, const real threshold_dev,
     const uint32_t meteoroids_per_thread,
     const size_t blocks_num, const size_t threads_per_block);