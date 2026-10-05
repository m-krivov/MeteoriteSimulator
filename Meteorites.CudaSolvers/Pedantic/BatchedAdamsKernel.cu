#include "BatchedAdamsKernel.h"

#include "Meteorites.Core/Constants.h"


namespace
{

template <unsigned int STEPS>
__device__ void InitContext(ThreadContext<STEPS> &ctx, const VirtualMeteoroid &meteoroid, real dt,
                            Record *&record)
{
  Adams::Unchangeable params(meteoroid);
  ctx.curr_layer = { meteoroid.V0, meteoroid.Gamma0, meteoroid.h0, meteoroid.l0, meteoroid.M0 };
  *record = Record(0.0, ctx.curr_layer);
  record++;

  Adams::OneStepIteration(ctx.curr_layer, ctx.steps[0], params, dt);
  *record = Record(dt, ctx.curr_layer);
  record++;

  if constexpr (STEPS >= 2) {
    Adams::TwoStepIteration(ctx.curr_layer, ctx.steps[1], ctx.steps[0], params, dt);
    *record = Record(dt * 2, ctx.curr_layer);
    record++;
  }

  if constexpr (STEPS >= 3) {
    Adams::ThreeStepIteration(ctx.curr_layer, ctx.steps[2], ctx.steps[1], ctx.steps[0], params, dt);
    *record = Record(dt * 3, ctx.curr_layer);
    record++;
  }
}

template <unsigned int STEPS>
__global__ void AdamsKernel(ThreadContext<STEPS> *contexts, int32_t *active_threads,
                            const VirtualMeteoroid *problems, size_t n_problems,
                            real dt, real timeout,
                            const real *timestamps, size_t n_timestamps,
                            real *functional_args, Record *records,
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
  real *V_arg, *h_arg;
  size_t nxt, timestamp;
  uint32_t iters_count;

  // Restore context
  Adams::Unchangeable params(meteoroid);
  if (ctx.t == 0.0) // new meteoroid
  {
    InitContext(ctx, meteoroid, dt, record);
    t = dt * (real)STEPS;
    nxt = 0;
    timestamp = 0;
    iters_count = STEPS + 1;
  }
  else // previous meteoroid
  {
    t = ctx.t;
    nxt = ctx.nxt;
    timestamp = ctx.timestamp;
    iters_count = 0;
  }
  V_arg = functional_args + (idx * n_timestamps * 2);
  h_arg = functional_args + (idx * n_timestamps * 2) + n_timestamps;

  // The main loop
  while (iters_count < iterations)
  {
    // If necessery, update the functional's arguments
    if (timestamp < n_timestamps && t >= timestamps[timestamp])
    {
      V_arg[timestamp] = ctx.curr_layer.V;
      h_arg[timestamp] = ctx.curr_layer.h;
      timestamp++;
    }

    // Compute values for the next step, store them
    Adams::Iteration<STEPS>(ctx.curr_layer, ctx.steps, params, nxt, dt);

    t += dt;
    *record = Record(t, ctx.curr_layer);

    // Should we stop the simulation?
    if (ctx.curr_layer.M <= (real)0.01 || ctx.curr_layer.h <= (real)0.0 || t >= timeout)
    {
      record->t = 0.0; // stop marker
      ctx.t = 0.0;
      ctx.ended = true;
      atomicSub(active_threads, 1);
      return;
    }
    iters_count++;
    record++;
    nxt = (nxt + 1) % STEPS;
  }

  // Update context
  ctx.t = t;
  ctx.nxt = nxt;
  ctx.timestamp = timestamp;
};

} // unnamed namespace


template <unsigned int STEPS>
void BatchedAdamsKernel(ThreadContext<STEPS> *contexts, int32_t *active_threads,
                        const VirtualMeteoroid *problems, size_t n_problems,
                        real dt, real timeout, const real *timestamps, size_t n_timestamps,
                        real *functional_args, Record *records,
                        size_t iterations, size_t threads_per_block, cudaStream_t stream)
{
  assert(n_problems > 0);
  assert(iterations > 0);
  assert(threads_per_block > 0);

  dim3 threads{ (uint32_t)threads_per_block };
  dim3 blocks(((n_problems - 1) / threads.x) + 1);
  AdamsKernel<STEPS><<<blocks, threads, 0, stream>>>(contexts, active_threads,
                                                     problems, n_problems, dt, timeout,
                                                     timestamps, n_timestamps,
                                                     functional_args, records, iterations);
  HANDLE_ERROR(cudaGetLastError());
}

template
void BatchedAdamsKernel<1u>(ThreadContext<1u> *contexts, int32_t *active_threads,
                            const VirtualMeteoroid *problems, size_t n_problems, real dt, real timeout,
                            const real *timestamps, size_t n_timestamps,
                            real *functional_args, Record *records,
                            size_t iterations, size_t threads_per_block, cudaStream_t stream);
template
void BatchedAdamsKernel<2u>(ThreadContext<2u> *contexts, int32_t *active_threads,
                            const VirtualMeteoroid *problems, size_t n_problems, real dt, real timeout,
                            const real *timestamps, size_t n_timestamps,
                            real *functional_args, Record *records,
                            size_t iterations, size_t threads_per_block, cudaStream_t stream);
template
void BatchedAdamsKernel<3u>(ThreadContext<3u> *contexts, int32_t *active_threads,
                            const VirtualMeteoroid *problems, size_t n_problems, real dt, real timeout,
                            const real *timestamps, size_t n_timestamps,
                            real *functional_args, Record *records,
                            size_t iterations, size_t threads_per_block, cudaStream_t stream);
