#include "AdamsMethod.h"

#include "Meteorites.Core/ResultFormatters.h"

#include "CudaStructures.h"

template <uint32_t STEPS>
__device__ void
UniStepAdams(Layer *steps, CudaCase& params, size_t nxt, real dt);

__device__ void
OneStepAdams(Layer& res, const Layer& f0, real dt, CudaCase& params)
{
  res.Set(f0.V_ + f0.fV_ * dt, f0.gamma_ + f0.fgamma_ * dt, f0.h_ + f0.fh_ * dt, f0.l_ + f0.fl_ * dt,
          f0.M_ + f0.fM_ * dt, params);
}

template<>
__device__ void
UniStepAdams<1>(Layer *steps, CudaCase& params, size_t nxt, real dt)
{
  OneStepAdams(steps[nxt], steps[(nxt + 1) & 1], dt, params);
}

__device__ void
TwoStepAdams(Layer& res, const Layer& f1, const Layer& f0, real dt, CudaCase& params)
{
  const auto c1 = (real)1.5;
  const auto c0 = -(real)0.5;

  res.Set(f1.V_ + (c1 * f1.fV_ + c0 * f0.fV_) * dt, f1.gamma_ + (c1 * f1.fgamma_ + c0 * f0.fgamma_) * dt,
          f1.h_ + (c1 * f1.fh_ + c0 * f0.fh_) * dt, f1.l_ + (c1 * f1.fl_ + c0 * f0.fl_) * dt,
          f1.M_ + (c1 * f1.fM_ + c0 * f0.fM_) * dt, params);
}

template <>
__device__ void
UniStepAdams<2>(Layer *steps, CudaCase& params, size_t nxt, real dt)
{
  TwoStepAdams(steps[nxt], steps[(nxt + 1) % 3], steps[(nxt + 2) % 3], dt, params);
}

__device__ void
ThreeStepAdams(Layer& res, const Layer& f2, const Layer& f1, const Layer& f0, real dt, CudaCase& params)
{
  const auto c2 = (real)23 / 12;
  const auto c1 = -(real)16 / 12;
  const auto c0 = (real)5 / 12;

  res.Set(f2.V_ + (c2 * f2.fV_ + c1 * f1.fV_ + c0 * f0.fV_) * dt,
          f2.gamma_ + (c2 * f2.fgamma_ + c1 * f1.fgamma_ + c0 * f0.fgamma_) * dt,
          f2.h_ + (c2 * f2.fh_ + c1 * f1.fh_ + c0 * f0.fh_) * dt,
          f2.l_ + (c2 * f2.fl_ + c1 * f1.fl_ + c0 * f0.fl_) * dt,
          f2.M_ + (c2 * f2.fM_ + c1 * f1.fM_ + c0 * f0.fM_) * dt, params);
}

template <>
__device__ void
UniStepAdams<3>(Layer *steps, CudaCase& params, size_t nxt, real dt)
{
  ThreeStepAdams(steps[nxt], steps[(nxt + 1) & 3], steps[(nxt + 2) & 3], steps[(nxt + 3) & 3], dt, params);
}

template <uint32_t STEPS>
__device__ void
InitSteps(Layer *steps, CudaCase& curr_case, real dt)
{
  steps[STEPS].Set(curr_case.v0_, curr_case.gamma0_, curr_case.h0_, 0.0, // Case::l0()
                     curr_case.m0_, curr_case);

  OneStepAdams(steps[STEPS - 1], steps[STEPS], dt, curr_case);

  if (STEPS >= 2)
  {
    TwoStepAdams(steps[STEPS - 2], steps[STEPS - 1], steps[STEPS], dt, curr_case);
  }
  if (STEPS >= 3)
  {
    ThreeStepAdams(steps[STEPS - 3], steps[STEPS - 2], steps[STEPS - 1], steps[STEPS], dt, curr_case);
  }
}

template <uint32_t STEPS>
__device__ inline IResultFormatter::Reason
AdamsStep(const real *timestamps, size_t n_timestamps, real &dt, real &timeout,
          Layer *steps, real *V_arg, real *h_arg, size_t &timestamp,
          size_t &nxt, CudaCase &current_case, real &t)
{
  if (timestamp < n_timestamps && t >= timestamps[timestamp])
  {
    const auto &step = steps[(nxt + 1) % (STEPS + 1)];
    V_arg[timestamp] = step.V_;
    h_arg[timestamp] = step.h_;
    timestamp++;
  }

  UniStepAdams<STEPS>(steps, current_case, nxt, dt);
  auto M = steps[nxt].M_;
  auto h = steps[nxt].h_;

  t += dt;

  nxt = (nxt + STEPS) % (STEPS + 1);

  if (M <= (real)0.01)
  {
    return IResultFormatter::Reason::Burnt;
  }
  if (h <= (real)0.0)
  {
    return IResultFormatter::Reason::Collided;
  }
  if (t >= timeout)
  {
    return IResultFormatter::Reason::Timeouted;
  }
  return IResultFormatter::Reason::NA;
}

__device__ inline real
GPUL2Compute(real *v, real *h, const real *ref_v, const real *ref_h, int size)
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

__device__ inline void
InsertToTopSmallest(real *devs, curandState *states,
                    real dev, curandState state,
                    real border_dev = std::numeric_limits<real>::max())
{
    if (dev > devs[BEST_CASES_BUFFER_SIZE-1] || dev > border_dev) return;
    int pos = 0;
    while (pos < BEST_CASES_BUFFER_SIZE && dev > devs[pos]) pos++;

    if (pos < BEST_CASES_BUFFER_SIZE)
    {
        for (int i = BEST_CASES_BUFFER_SIZE-1; i > pos; i--)
        {
            devs[i] = devs[i-1];
            states[i] = states[i-1];
        }
        devs[pos] = dev;
        states[pos] = state;
    }
}

__device__ inline void
GenerateCase(CudaCase &curr_case, curandState &curr_state, const real v0, const real h0)
{
  curr_case.H_ = 1e5 + curand_uniform(&curr_state) * (5e6 - 1e5);
  curr_case.ch_ = 0.1 + curand_uniform(&curr_state) * (0.9 - 0.1);
  curr_case.rho_ = 2000.0 + curand_uniform(&curr_state) * (5000.0 - 2000.0);
  curr_case.cd_ = 0.5 + curand_uniform(&curr_state) * (2.5 - 0.5);
  curr_case.cl_ = 0.0 + curand_uniform(&curr_state) * (0.25 - 0.0);
  curr_case.m0_ = 10.0 + curand_uniform(&curr_state) * (500.0 - 10.0);
  curr_case.v0_ = v0;// + (curand_uniform(&curr_state) - 0.5) * v0 * 0.001;
  curr_case.h0_ = h0;// + (curand_uniform(&curr_state) - 0.5) * h0 * 0.001;
  curr_case.gamma0_ = curand_uniform(&curr_state) * (M_PI / 2);
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

template <uint32_t STEPS>
__global__ void
AdamsMethodKernel(const uint64_t *seeds,
                  const real *timestamps, size_t n_timestamps,
                  const real *ref_v, const real *ref_h,
                  real dt, real timeout,
                  real *best_cases_devs, curandState *best_cases_states, real border_dev,
                  real *v_args_arr, real *h_args_arr,
                  uint32_t CASES_PER_THREAD)
{
  const uint32_t tid = blockIdx.x * blockDim.x + threadIdx.x;

  curandState base_state;
  curand_init(seeds[tid], 0, 0, &base_state);
  curandState curr_state = base_state;
  CudaCase    curr_case;
  const real  v0 = ref_v[0], h0 = ref_h[0];

  curandState local_best_cases_states[BEST_CASES_BUFFER_SIZE];
  real        local_best_cases_devs[BEST_CASES_BUFFER_SIZE];
  for (int i = 0; i < BEST_CASES_BUFFER_SIZE; i++)
  {
    local_best_cases_devs[i] = std::numeric_limits<double>::max();
  }
  real curr_dev;

  Layer  steps[STEPS + 1];
  real  *v_args = v_args_arr + tid * n_timestamps,
        *h_args = h_args_arr + tid * n_timestamps;
  real   t;
  size_t timestamp;
  size_t nxt;

  for (uint32_t i = 0; i < CASES_PER_THREAD; i++)
  {
    GenerateCase(curr_case, curr_state, v0, h0);

    t = dt * STEPS;
    timestamp = 0;
    nxt = STEPS;
    InitSteps<STEPS>(steps, curr_case, dt);

    while (true)
    {
#ifdef DISPLAY_WARP_DIVERGENCE
      unsigned int mask = __activemask();
      if(threadIdx.x == 0)
      {
        printf("ADAMS_STEP: ");
        for (int i = 31; i >= 0; i--) {
            printf("%d", (mask >> i) & 1);
            if (i % 8 == 0 && i != 0) printf(" ");
        }
        printf("\n");
      }
#endif
      if (AdamsStep<STEPS>(timestamps, n_timestamps, dt, timeout, steps, v_args, h_args,
                           timestamp, nxt, curr_case, t) != IResultFormatter::Reason::NA) break;
    }
#ifdef DISPLAY_WARP_DIVERGENCE
    unsigned int mask = __activemask();
    if(threadIdx.x == 0)
    {
      printf("UPDATE: ");
      for (int i = 31; i >= 0; i--) {
          printf("%d", (mask >> i) & 1);
          if (i % 8 == 0 && i != 0) printf(" ");
      }
      printf("\n");
    }
#endif
    curr_dev = GPUL2Compute(&v_args[0], &h_args[0], ref_v, ref_h, n_timestamps);
    InsertToTopSmallest(local_best_cases_devs, local_best_cases_states, curr_dev, base_state, border_dev);
    base_state = curr_state;
  }

  for (int i = 0; i < BEST_CASES_BUFFER_SIZE; i++)
  {
    best_cases_devs[tid * BEST_CASES_BUFFER_SIZE + i] = local_best_cases_devs[i];
    best_cases_states[tid * BEST_CASES_BUFFER_SIZE + i] = local_best_cases_states[i];
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

template <uint32_t STEPS>
__global__ void
AdamsMethodBalancedKernel(const uint64_t *seeds,
                          const real *timestamps, size_t n_timestamps,
                          const real *ref_v, const real *ref_h,
                          real dt, real timeout,
                          real *best_cases_devs, curandState *best_cases_states, real border_dev,
                          real *v_args_arr, real *h_args_arr,
                          uint32_t CASES_PER_THREAD)
{
  __shared__ uint32_t cases_counter;
  if (threadIdx.x == 0) cases_counter = blockDim.x;
  uint32_t case_number = blockDim.x;

  const uint32_t tid = blockIdx.x * blockDim.x + threadIdx.x;

  curandState base_state;
  curand_init(seeds[tid], 0, 0, &base_state);
  curandState curr_state = base_state;
  CudaCase    curr_case;
  const real  v0 = ref_v[0], h0 = ref_h[0];

  curandState local_best_cases_states[BEST_CASES_BUFFER_SIZE];
  real        local_best_cases_devs[BEST_CASES_BUFFER_SIZE];
  for (int i = 0; i < BEST_CASES_BUFFER_SIZE; i++)
  {
    local_best_cases_devs[i] = std::numeric_limits<double>::max();
  }
  real curr_dev;

  Layer  steps[STEPS + 1];
  real  *v_args = v_args_arr + tid * n_timestamps,
        *h_args = h_args_arr + tid * n_timestamps;
  real   t;
  size_t timestamp;
  size_t nxt;

  GenerateCase(curr_case, curr_state, v0, h0);

  t = dt * STEPS;
  timestamp = 0;
  nxt = STEPS;
  InitSteps<STEPS>(steps, curr_case, dt);

  while (true)
  {
#ifdef DISPLAY_WARP_DIVERGENCE
    unsigned int mask = __activemask();
    if(threadIdx.x == 0)
    {
      printf("ADAMS_STEP: ");
      for (int i = 31; i >= 0; i--) {
          printf("%d", (mask >> i) & 1);
          if (i % 8 == 0 && i != 0) printf(" ");
      }
      printf("\n");
    }
#endif
    if (AdamsStep<STEPS>(timestamps, n_timestamps, dt, timeout, steps, v_args, h_args,
                         timestamp, nxt, curr_case, t) == IResultFormatter::Reason::NA)
    {
      continue;
    }
    else
    {
#ifdef DISPLAY_WARP_DIVERGENCE
      unsigned int mask = __activemask();
      if(threadIdx.x == 0)
      {
        printf("UPDATE: ");
        for (int i = 31; i >= 0; i--) {
            printf("%d", (mask >> i) & 1);
            if (i % 8 == 0 && i != 0) printf(" ");
        }
        printf("\n");
      }
#endif
      curr_dev = GPUL2Compute(&v_args[0], &h_args[0], ref_v, ref_h, n_timestamps);
      InsertToTopSmallest(local_best_cases_devs, local_best_cases_states, curr_dev, base_state, border_dev);
      base_state = curr_state;

      case_number = atomicAdd(&cases_counter, 1);
      if (case_number >= blockDim.x * CASES_PER_THREAD) break;

      GenerateCase(curr_case, curr_state, v0, h0);

      t = dt * STEPS;
      timestamp = 0;
      nxt = STEPS;
      InitSteps<STEPS>(steps, curr_case, dt);
    }
  }
  for (int i = 0; i < BEST_CASES_BUFFER_SIZE; i++)
  {
    best_cases_devs[tid * BEST_CASES_BUFFER_SIZE + i] = local_best_cases_devs[i];
    best_cases_states[tid * BEST_CASES_BUFFER_SIZE + i] = local_best_cases_states[i];
  }
}

template <uint32_t STEPS>
void CudaAdamsMethodLauncher(const uint64_t *seeds,
                             const real *timestamps, size_t n_timestamps,
                             const real *ref_v, const real *ref_h,
                             real dt, real timeout,
                             real *best_cases_devs, curandState *best_cases_states, real border_dev,
                             real *v_args_arr, real *h_args_arr,
                             uint32_t CASES_PER_THREAD,
                             uint32_t BLOCKS_NUM, uint32_t THREADS_PER_BLOCK)
{

  AdamsMethodBalancedKernel<STEPS><<<BLOCKS_NUM, THREADS_PER_BLOCK>>>
    (seeds,
     timestamps, n_timestamps,
     ref_v, ref_h,
     dt, timeout,
     best_cases_devs, best_cases_states, border_dev,
     v_args_arr, h_args_arr,
     CASES_PER_THREAD);
}

// template-specified functions compiles only this way
template void CudaAdamsMethodLauncher<1u>(const uint64_t *seeds,
                                          const real *timestamps, size_t n_timestamps,
                                          const real *ref_v, const real *ref_h,
                                          real dt, real timeout,
                                          real *best_cases_devs, curandState *best_cases_states, real border_dev,
                                          real *v_args_arr, real *h_args_arr,
                                          uint32_t CASES_PER_THREAD,
                                          uint32_t BLOCKS_NUM, uint32_t THREADS_PER_BLOCK);
template void CudaAdamsMethodLauncher<2u>(const uint64_t *seeds,
                                          const real *timestamps, size_t n_timestamps,
                                          const real *ref_v, const real *ref_h,
                                          real dt, real timeout,
                                          real *best_cases_devs, curandState *best_cases_states, real border_dev,
                                          real *v_args_arr, real *h_args_arr,
                                          uint32_t CASES_PER_THREAD,
                                          uint32_t BLOCKS_NUM, uint32_t THREADS_PER_BLOCK);
template void CudaAdamsMethodLauncher<3u>(const uint64_t *seeds,
                                          const real *timestamps, size_t n_timestamps,
                                          const real *ref_v, const real *ref_h,
                                          real dt, real timeout,
                                          real *best_cases_devs, curandState *best_cases_states, real border_dev,
                                          real *v_args_arr, real *h_args_arr,
                                          uint32_t CASES_PER_THREAD,
                                          uint32_t BLOCKS_NUM, uint32_t THREADS_PER_BLOCK);
