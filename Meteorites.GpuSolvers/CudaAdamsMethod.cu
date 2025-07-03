#include "CudaAdamsMethod.h"

#include "GPUParameters.h"
#include "CudaStructures.h"

template <unsigned int STEPS, unsigned int ITERS>
__device__ void UniStepAdams(ThreadContext<STEPS, ITERS>& c, size_t nxt, real dt);

__device__ void OneStepAdams(Layer& res, const Layer& f0, real dt, CudaCase& params)
{
  res.Set(f0.V_ + f0.fV_ * dt, f0.gamma_ + f0.fgamma_ * dt, f0.h_ + f0.fh_ * dt, f0.l_ + f0.fl_ * dt,
          f0.M_ + f0.fM_ * dt, params);
}

template <> __device__ void UniStepAdams(ThreadContext<1, ITERS_PER_KERNEL>& c, size_t nxt, real dt)
{
  OneStepAdams(c.steps[nxt], c.steps[(nxt + 1) & 1], dt, c.params);
}

__device__ void TwoStepAdams(Layer& res, const Layer& f1, const Layer& f0, real dt, CudaCase& params)
{
  const auto c1 = (real)1.5;
  const auto c0 = -(real)0.5;

  res.Set(f1.V_ + (c1 * f1.fV_ + c0 * f0.fV_) * dt, f1.gamma_ + (c1 * f1.fgamma_ + c0 * f0.fgamma_) * dt,
          f1.h_ + (c1 * f1.fh_ + c0 * f0.fh_) * dt, f1.l_ + (c1 * f1.fl_ + c0 * f0.fl_) * dt,
          f1.M_ + (c1 * f1.fM_ + c0 * f0.fM_) * dt, params);
}

template <> __device__ void UniStepAdams(ThreadContext<2, ITERS_PER_KERNEL>& c, size_t nxt, real dt)
{
  TwoStepAdams(c.steps[nxt], c.steps[(nxt + 1) % 3], c.steps[(nxt + 2) % 3], dt, c.params);
}

__device__ void ThreeStepAdams(Layer& res, const Layer& f2, const Layer& f1, const Layer& f0, real dt, CudaCase& params)
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

template <> __device__ void UniStepAdams(ThreadContext<3, ITERS_PER_KERNEL>& c, size_t nxt, real dt)
{
  ThreeStepAdams(c.steps[nxt], c.steps[(nxt + 1) & 3], c.steps[(nxt + 2) & 3], c.steps[(nxt + 3) & 3], dt, c.params);
}

template <unsigned int STEPS, unsigned int ITERS>
__device__ void InitContext(ThreadContext<STEPS, ITERS>& c, CudaCase& curr_case, real dt, size_t grid_idx,
                            Record*& curr_record)
{
  c.params = curr_case;
  c.steps[STEPS].Set(curr_case.v0_, curr_case.gamma0_, curr_case.h0_, 0.0, // Case::l0()
                     curr_case.m0_, c.params);
  *curr_record = {
      0.0, c.steps[STEPS].M_, c.steps[STEPS].V_, c.steps[STEPS].h_, c.steps[STEPS].l_, c.steps[STEPS].gamma_};
  curr_record++;

  OneStepAdams(c.steps[STEPS - 1], c.steps[STEPS], dt, c.params);
  *curr_record = {dt,
                  c.steps[STEPS - 1].M_,
                  c.steps[STEPS - 1].V_,
                  c.steps[STEPS - 1].h_,
                  c.steps[STEPS - 1].l_,
                  c.steps[STEPS - 1].gamma_};
  curr_record++;
  if (STEPS >= 2) {
    TwoStepAdams(c.steps[STEPS - 2], c.steps[STEPS - 1], c.steps[STEPS], dt, c.params);
    *curr_record = {dt * 2,
                    c.steps[STEPS - 2].M_,
                    c.steps[STEPS - 2].V_,
                    c.steps[STEPS - 2].h_,
                    c.steps[STEPS - 2].l_,
                    c.steps[STEPS - 2].gamma_};
    curr_record++;
  }

  if (STEPS >= 3) {
    ThreeStepAdams(c.steps[STEPS - 3], c.steps[STEPS - 2], c.steps[STEPS - 1], c.steps[STEPS], dt, c.params);
    *curr_record = {dt * 3,
                    c.steps[STEPS - 3].M_,
                    c.steps[STEPS - 3].V_,
                    c.steps[STEPS - 3].h_,
                    c.steps[STEPS - 3].l_,
                    c.steps[STEPS - 3].gamma_};
    curr_record++;
  }
}

template <unsigned int STEPS, unsigned int ITERS>
__global__ void AdamsKernel(ThreadContext<STEPS, ITERS>* contexts, CudaCase* problems, real dt, real timeout,
                            real* timestamps, size_t n_timestamps, real* functional_args, Record* records,
                            unsigned int* active_threads)
{
  size_t grid_idx = blockIdx.x * blockDim.x + threadIdx.x;
  ThreadContext<STEPS, ITERS>& c = contexts[grid_idx];

  if (c.ended) {
    return;
  }
  if (c.curr_case_num == CASES_PER_THREAD) {
    c.ended = true;
    atomicSub(active_threads, 1);
    return;
  }

  CudaCase& current_case = problems[grid_idx * CASES_PER_THREAD + c.curr_case_num];
  Record* curr_record = records + ITERS_PER_KERNEL * grid_idx;
  real t;
  real *V_arg, *h_arg;
  size_t nxt, timestamp;
  unsigned int iters_count;

  // restore context
  if (c.t == 0.0) { // new case
    InitContext(c, current_case, dt, grid_idx, curr_record);
    t = dt * (real)STEPS;
    nxt = STEPS;
    timestamp = 0;
    iters_count = STEPS + 1;
  } else { // old case
    t = c.t;
    nxt = c.nxt;
    timestamp = c.timestamp;
    iters_count = 0;
  }
  V_arg = functional_args + (grid_idx * n_timestamps * 2 * CASES_PER_THREAD) + (n_timestamps * 2 * c.curr_case_num);
  h_arg = functional_args + (grid_idx * n_timestamps * 2 * CASES_PER_THREAD) + (n_timestamps * 2 * c.curr_case_num) +
          n_timestamps;

  // main loop
  while (iters_count < ITERS) {
    // If necessery, update the functional's arguments
    if (timestamp < n_timestamps && t >= timestamps[timestamp]) {
      const auto& step = c.steps[(nxt + 1) % (STEPS + 1)];
      V_arg[timestamp] = step.V_;
      h_arg[timestamp] = step.h_;
      timestamp++;
    }

    // Compute values for the next step, store them
    UniStepAdams<STEPS, ITERS>(c, nxt, dt);
    auto M = c.steps[nxt].M_;
    auto h = c.steps[nxt].h_;

    t += dt;
    *curr_record = {t, M, c.steps[nxt].V_, h, c.steps[nxt].l_, c.steps[nxt].gamma_};

    // Check, should we stop the simulation?
    if (M <= (real)0.01 || h <= (real)0.0 || t >= timeout) {
      curr_record->t_ = 0.0; // stop marker
      c.t = 0.0;
      c.curr_case_num++;
      return;
    }
    iters_count++;
    curr_record++;
    nxt = (nxt + STEPS) % (STEPS + 1);
  }
  // update context
  c.t = t;
  c.nxt = nxt;
  c.timestamp = timestamp;
};

template <unsigned int STEPS, unsigned int ITERS>
void CudaLauncher(ThreadContext<STEPS, ITERS>* dev_thread_sandbox_arr, CudaCase* dev_problems,
                  size_t promblems_vector_size, real dt, real timeout, real* dev_timestamps, size_t n_timestamps,
                  real* dev_functional_args, Record* dev_records, unsigned int* dev_active_threads,
                  unsigned int launching_blocks, unsigned int launching_threads_per_block)
{

  AdamsKernel<STEPS, ITERS><<<launching_blocks, launching_threads_per_block>>>(
      dev_thread_sandbox_arr, dev_problems, dt, timeout, dev_timestamps, n_timestamps, dev_functional_args, dev_records,
      dev_active_threads);
}

// template-specified functions compiles only this way
template void CudaLauncher<1u, ITERS_PER_KERNEL>(ThreadContext<1u, ITERS_PER_KERNEL>* dev_thread_sandbox_arr,
                                                 CudaCase* dev_problems, size_t promblems_vector_size, real dt,
                                                 real timeout, real* dev_timestamps, size_t n_timestamps,
                                                 real* dev_functional_args, Record* dev_records,
                                                 unsigned int* dev_active_threads, unsigned int launching_blocks,
                                                 unsigned int launching_threads_per_block);
template void CudaLauncher<2u, ITERS_PER_KERNEL>(ThreadContext<2u, ITERS_PER_KERNEL>* dev_thread_sandbox_arr,
                                                 CudaCase* dev_problems, size_t promblems_vector_size, real dt,
                                                 real timeout, real* dev_timestamps, size_t n_timestamps,
                                                 real* dev_functional_args, Record* dev_records,
                                                 unsigned int* dev_active_threads, unsigned int launching_blocks,
                                                 unsigned int launching_threads_per_block);
template void CudaLauncher<3u, ITERS_PER_KERNEL>(ThreadContext<3u, ITERS_PER_KERNEL>* dev_thread_sandbox_arr,
                                                 CudaCase* dev_problems, size_t promblems_vector_size, real dt,
                                                 real timeout, real* dev_timestamps, size_t n_timestamps,
                                                 real* dev_functional_args, Record* dev_records,
                                                 unsigned int* dev_active_threads, unsigned int launching_blocks,
                                                 unsigned int launching_threads_per_block);
