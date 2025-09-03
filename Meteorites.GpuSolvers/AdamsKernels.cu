#include "AdamsKernels.h"

#include "Meteorites.Core/Constants.h"

#include "GPUParameters.h"

template <unsigned int STEPS>
__device__ void InitContext(ThreadContext<STEPS> &ctx, const Case &meteoroid, real dt, size_t idx,
                            Record *&record)
{
  Adams::Unchangeable params(meteoroid);
  Adams::SetLayer(ctx.steps[STEPS], params,
                  meteoroid.V0, meteoroid.Gamma0, meteoroid.h0, meteoroid.l0, meteoroid.M0);
  *record = { 0.0,
              ctx.steps[STEPS].M,
              ctx.steps[STEPS].V,
              ctx.steps[STEPS].h,
              ctx.steps[STEPS].l,
              ctx.steps[STEPS].Gamma };
  record++;

  Adams::OneStepIteration(ctx.steps[STEPS - 1], ctx.steps[STEPS], params, dt);
  *record = { dt,
              ctx.steps[STEPS - 1].M,
              ctx.steps[STEPS - 1].V,
              ctx.steps[STEPS - 1].h,
              ctx.steps[STEPS - 1].l,
              ctx.steps[STEPS - 1].Gamma };
  record++;
  if constexpr (STEPS >= 2) {
    Adams::TwoStepIteration(ctx.steps[STEPS - 2], ctx.steps[STEPS - 1],
                            ctx.steps[STEPS], params, dt);
    *record = { dt * 2,
                ctx.steps[STEPS - 2].M,
                ctx.steps[STEPS - 2].V,
                ctx.steps[STEPS - 2].h,
                ctx.steps[STEPS - 2].l,
                ctx.steps[STEPS - 2].Gamma };
    record++;
  }

  if constexpr (STEPS >= 3) {
    Adams::ThreeStepIteration(ctx.steps[STEPS - 3], ctx.steps[STEPS - 2],
                              ctx.steps[STEPS - 1], ctx.steps[STEPS], params, dt);
    *record = { dt * 3,
                ctx.steps[STEPS - 3].M,
                ctx.steps[STEPS - 3].V,
                ctx.steps[STEPS - 3].h,
                ctx.steps[STEPS - 3].l,
                ctx.steps[STEPS - 3].Gamma };
    record++;
  }
}

template <unsigned int STEPS>
__global__ void AdamsKernel(ThreadContext<STEPS> *contexts, uint32_t *active_threads,
                            const Case *problems, real dt, real timeout,
                            const real *timestamps, size_t n_timestamps,
                            real *functional_args, Record *records,
                            size_t iterations)
{
  size_t idx = blockIdx.x * blockDim.x + threadIdx.x;

  if (idx * CASES_PER_THREAD >= CASE_NUM)
  { return; }

  ThreadContext<STEPS>& ctx = contexts[idx];
  if (ctx.ended)
  { return; }

  if (ctx.curr_case_num == CASES_PER_THREAD)
  {
    ctx.ended = true;
    atomicSub(active_threads, 1);
    return;
  }

  const Case &meteoroid = problems[idx * CASES_PER_THREAD + ctx.curr_case_num];
  Record *record = records + iterations * idx;
  real t;
  real *V_arg, *h_arg;
  size_t nxt, timestamp;
  uint32_t iters_count;

  // Restore context
  Adams::Unchangeable params(meteoroid);
  if (ctx.t == 0.0) // new case
  {
    InitContext(ctx, meteoroid, dt, idx, record);
    t = dt * (real)STEPS;
    nxt = STEPS;
    timestamp = 0;
    iters_count = STEPS + 1;
  }
  else // old case
  {
    t = ctx.t;
    nxt = ctx.nxt;
    timestamp = ctx.timestamp;
    iters_count = 0;
  }
  V_arg = functional_args + (idx * n_timestamps * 2 * CASES_PER_THREAD) + (n_timestamps * 2 * ctx.curr_case_num);
  h_arg = functional_args + (idx * n_timestamps * 2 * CASES_PER_THREAD) + (n_timestamps * 2 * ctx.curr_case_num) +
          n_timestamps;

  // The main loop
  while (iters_count < iterations)
  {
    // If necessery, update the functional's arguments
    if (timestamp < n_timestamps && t >= timestamps[timestamp])
    {
      const auto& step = ctx.steps[(nxt + 1) % (STEPS + 1)];
      V_arg[timestamp] = step.V;
      h_arg[timestamp] = step.h;
      timestamp++;
    }

    // Compute values for the next step, store them
    Adams::Iteration<STEPS>(ctx.steps, params, nxt, dt);
    auto M = ctx.steps[nxt].M;
    auto h = ctx.steps[nxt].h;

    t += dt;
    *record = {t, M, ctx.steps[nxt].V, h, ctx.steps[nxt].l, ctx.steps[nxt].Gamma};

    // Should we stop the simulation?
    if (M <= (real)0.01 || h <= (real)0.0 || t >= timeout)
    {
      record->t = 0.0; // stop marker
      ctx.t = 0.0;
      ctx.curr_case_num++;
      return;
    }
    iters_count++;
    record++;
    nxt = (nxt + STEPS) % (STEPS + 1);
  }

  // Update context
  ctx.t = t;
  ctx.nxt = nxt;
  ctx.timestamp = timestamp;
};

template <unsigned int STEPS>
void BatchedAdamsKernel(ThreadContext<STEPS> *contexts, uint32_t *active_threads,
                        const Case *problems, size_t n_problems,
                        real dt, real timeout, const real *timestamps, size_t n_timestamps,
                        real *functional_args, Record *records,
                        size_t iterations)
{
  assert(n_problems > 0);
  dim3 threads_per_block{ THREADS_PER_BLOCK };
  dim3 blocks(((n_problems - 1) / CASES_PER_THREAD / threads_per_block.x) + 1);
  AdamsKernel<STEPS><<<blocks, threads_per_block>>>(
      contexts, active_threads, problems, dt, timeout, timestamps, n_timestamps, functional_args, records,
      iterations);
}

template
void BatchedAdamsKernel<1u>(ThreadContext<1u> *contexts, uint32_t *active_threads,
                            const Case *problems, size_t n_problems, real dt, real timeout,
                            const real *timestamps, size_t n_timestamps,
                            real *functional_args, Record *records,
                            size_t iterations);
template
void BatchedAdamsKernel<2u>(ThreadContext<2u> *contexts, uint32_t *active_threads,
                            const Case *problems, size_t n_problems, real dt, real timeout,
                            const real *timestamps, size_t n_timestamps,
                            real *functional_args, Record *records,
                            size_t iterations);
template
void BatchedAdamsKernel<3u>(ThreadContext<3u> *contexts, uint32_t *active_threads,
                            const Case *problems, size_t n_problems, real dt, real timeout,
                            const real *timestamps, size_t n_timestamps,
                            real *functional_args, Record *records,
                            size_t iterations);
