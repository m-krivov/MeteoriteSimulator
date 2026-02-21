#include "BatchedAdamsKernel.h"

#include "Meteorites.Core/Constants.h"


namespace
{

template <unsigned int STEPS>
__device__ void InitContext(ThreadContext<STEPS> &ctx, const VirtualMeteoroid &meteoroid, real dt, size_t idx,
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
__global__ void AdamsKernel(ThreadContext<STEPS> *contexts, int32_t *active_threads,
                            const VirtualMeteoroid *problems, size_t n_problems,
                            real dt, real timeout,
                            Record *records,
                            size_t iterations)
{
  size_t idx = blockIdx.x * blockDim.x + threadIdx.x;

  if (idx >= n_problems)
  { return; }

  ThreadContext<STEPS> &ctx = contexts[idx];
  if (ctx.ended)
  { return; }

  const VirtualMeteoroid &meteoroid = problems[idx];
  Record *record = records + iterations * idx;
  real t;
  size_t nxt;
  uint32_t iters_count;

  // Restore context
  Adams::Unchangeable params(meteoroid);
  if (ctx.t == 0.0) // new meteoroid
  {
    InitContext(ctx, meteoroid, dt, idx, record);
    t = dt * (real)STEPS;
    nxt = STEPS;
    iters_count = STEPS + 1;
  }
  else // previous meteoroid
  {
    t = ctx.t;
    nxt = ctx.nxt;
    iters_count = 0;
  }

  // The main loop
  while (iters_count < iterations)
  {
    // Compute values for the next step, store them
    Adams::Iteration<STEPS>(ctx.steps, params, nxt, dt);
    auto M = ctx.steps[nxt].M;
    auto h = ctx.steps[nxt].h;

    t += dt;
    *record = {t, M, ctx.steps[nxt].V, h, ctx.steps[nxt].l, ctx.steps[nxt].Gamma};

    // Should we stop the simulation?
    if (M <= (real)0.01 || h <= (real)0.0 || t >= timeout)
    {
      record->t = -1.0; // stop marker
      ctx.t = 0.0;
      ctx.ended = true;
      atomicSub(active_threads, 1);
      return;
    }
    iters_count++;
    record++;
    nxt = (nxt + STEPS) % (STEPS + 1);
  }

  // Update context
  ctx.t = t;
  ctx.nxt = nxt;
};

} // unnamed namespace


template <unsigned int STEPS>
void BatchedAdamsKernel(ThreadContext<STEPS> *contexts, int32_t *active_threads,
                        const VirtualMeteoroid *problems, size_t n_problems,
                        real dt, real timeout,
                        Record *records,
                        size_t iterations, size_t threads_per_block, cudaStream_t stream)
{
  assert(n_problems > 0);
  assert(iterations > 0);
  assert(threads_per_block > 0);

  dim3 threads{ (uint32_t)threads_per_block };
  dim3 blocks(((n_problems - 1) / threads.x) + 1);
  AdamsKernel<STEPS><<<blocks, threads, 0, stream>>>(contexts, active_threads,
                                                     problems, n_problems, dt, timeout,
                                                     records, iterations);
  HANDLE_ERROR(cudaGetLastError());
}

template
void BatchedAdamsKernel<1u>(ThreadContext<1u> *contexts, int32_t *active_threads,
                            const VirtualMeteoroid *problems, size_t n_problems, real dt, real timeout,
                            Record *records,
                            size_t iterations, size_t threads_per_block, cudaStream_t stream);
template
void BatchedAdamsKernel<2u>(ThreadContext<2u> *contexts, int32_t *active_threads,
                            const VirtualMeteoroid *problems, size_t n_problems, real dt, real timeout,
                            Record *records,
                            size_t iterations, size_t threads_per_block, cudaStream_t stream);
template
void BatchedAdamsKernel<3u>(ThreadContext<3u> *contexts, int32_t *active_threads,
                            const VirtualMeteoroid *problems, size_t n_problems, real dt, real timeout,
                            Record *records,
                            size_t iterations, size_t threads_per_block, cudaStream_t stream);
