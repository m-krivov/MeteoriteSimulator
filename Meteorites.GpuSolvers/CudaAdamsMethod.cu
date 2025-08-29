#include "CudaAdamsMethod.h"

#include "Meteorites.Core/Adams.h"

#include "GPUParameters.h"
#include "CudaStructures.h"

template <unsigned int STEPS>
__device__ void InitContext(ThreadContext<STEPS>& c, Case& curr_case, real dt, size_t grid_idx,
                            Record*& curr_record)
{
  Adams::Unchangeable params(curr_case);
  Adams::SetLayer(c.steps[STEPS], params,
                  curr_case.V0, curr_case.Gamma0, curr_case.h0, curr_case.l0, curr_case.M0);
  *curr_record = {
      0.0, c.steps[STEPS].M, c.steps[STEPS].V, c.steps[STEPS].h, c.steps[STEPS].l, c.steps[STEPS].Gamma};
  curr_record++;

  Adams::OneStepIteration(c.steps[STEPS - 1], c.steps[STEPS], params, dt);
  *curr_record = {dt,
                  c.steps[STEPS - 1].M,
                  c.steps[STEPS - 1].V,
                  c.steps[STEPS - 1].h,
                  c.steps[STEPS - 1].l,
                  c.steps[STEPS - 1].Gamma};
  curr_record++;
  if constexpr (STEPS >= 2) {
    Adams::TwoStepIteration(c.steps[STEPS - 2], c.steps[STEPS - 1], c.steps[STEPS], params, dt);
    *curr_record = {dt * 2,
                    c.steps[STEPS - 2].M,
                    c.steps[STEPS - 2].V,
                    c.steps[STEPS - 2].h,
                    c.steps[STEPS - 2].l,
                    c.steps[STEPS - 2].Gamma};
    curr_record++;
  }

  if constexpr (STEPS >= 3) {
    Adams::ThreeStepIteration(c.steps[STEPS - 3], c.steps[STEPS - 2], c.steps[STEPS - 1], c.steps[STEPS], params, dt);
    *curr_record = {dt * 3,
                    c.steps[STEPS - 3].M,
                    c.steps[STEPS - 3].V,
                    c.steps[STEPS - 3].h,
                    c.steps[STEPS - 3].l,
                    c.steps[STEPS - 3].Gamma};
    curr_record++;
  }
}

template <unsigned int STEPS, unsigned int ITERS>
__global__ void AdamsKernel(ThreadContext<STEPS>* contexts, Case* problems, real dt, real timeout,
                            real* timestamps, size_t n_timestamps, real* functional_args, Record* records,
                            unsigned int* active_threads)
{
  size_t grid_idx = blockIdx.x * blockDim.x + threadIdx.x;

  if (grid_idx * CASES_PER_THREAD >= CASE_NUM) return;

  ThreadContext<STEPS>& c = contexts[grid_idx];

  if (c.ended) {
    return;
  }
  if (c.curr_case_num == CASES_PER_THREAD) {
    c.ended = true;
    atomicSub(active_threads, 1);
    return;
  }

  Case& current_case = problems[grid_idx * CASES_PER_THREAD + c.curr_case_num];
  Record* curr_record = records + ITERS_PER_KERNEL * grid_idx;
  real t;
  real *V_arg, *h_arg;
  size_t nxt, timestamp;
  unsigned int iters_count;

  // restore context
  Adams::Unchangeable params(current_case);
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
      V_arg[timestamp] = step.V;
      h_arg[timestamp] = step.h;
      timestamp++;
    }

    // Compute values for the next step, store them
    Adams::Iteration<STEPS>(c.steps, params, nxt, dt);
    auto M = c.steps[nxt].M;
    auto h = c.steps[nxt].h;

    t += dt;
    *curr_record = {t, M, c.steps[nxt].V, h, c.steps[nxt].l, c.steps[nxt].Gamma};

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
void CudaLauncher(ThreadContext<STEPS>* dev_thread_sandbox_arr, Case* dev_problems,
                  size_t promblems_vector_size, real dt, real timeout, real* dev_timestamps, size_t n_timestamps,
                  real* dev_functional_args, Record* dev_records, unsigned int* dev_active_threads,
                  unsigned int launching_blocks, unsigned int launching_threads_per_block)
{

  AdamsKernel<STEPS, ITERS><<<launching_blocks, launching_threads_per_block>>>(
      dev_thread_sandbox_arr, dev_problems, dt, timeout, dev_timestamps, n_timestamps, dev_functional_args, dev_records,
      dev_active_threads);
}

// template-specified functions compiles only this way
template void CudaLauncher<1u, ITERS_PER_KERNEL>(ThreadContext<1u>* dev_thread_sandbox_arr,
                                                 Case* dev_problems, size_t promblems_vector_size, real dt,
                                                 real timeout, real* dev_timestamps, size_t n_timestamps,
                                                 real* dev_functional_args, Record* dev_records,
                                                 unsigned int* dev_active_threads, unsigned int launching_blocks,
                                                 unsigned int launching_threads_per_block);
template void CudaLauncher<2u, ITERS_PER_KERNEL>(ThreadContext<2u>* dev_thread_sandbox_arr,
                                                 Case* dev_problems, size_t promblems_vector_size, real dt,
                                                 real timeout, real* dev_timestamps, size_t n_timestamps,
                                                 real* dev_functional_args, Record* dev_records,
                                                 unsigned int* dev_active_threads, unsigned int launching_blocks,
                                                 unsigned int launching_threads_per_block);
template void CudaLauncher<3u, ITERS_PER_KERNEL>(ThreadContext<3u>* dev_thread_sandbox_arr,
                                                 Case* dev_problems, size_t promblems_vector_size, real dt,
                                                 real timeout, real* dev_timestamps, size_t n_timestamps,
                                                 real* dev_functional_args, Record* dev_records,
                                                 unsigned int* dev_active_threads, unsigned int launching_blocks,
                                                 unsigned int launching_threads_per_block);
