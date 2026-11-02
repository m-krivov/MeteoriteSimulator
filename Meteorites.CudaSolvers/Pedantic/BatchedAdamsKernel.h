#pragma once
#include "Meteorites.CudaSolvers/CudaDefs.h"

#include "Meteorites.Core/Meteoroids/VirtualMeteoroid.h"
#include "Meteorites.Core/Solvers/Adams.h"


// Contains a few values extracted from Adams::State
// They could be passed to formatters
struct Record
{
  real t{}, M{}, V{}, h{}, l{}, Gamma{};

  Record() = default;
  DEVICE Record(real t, const Adams::State &state)
    : t(t), M(state.M), V(state.V), h(state.h), l(state.l), Gamma(state.Gamma)
  { }
  Record(const Record &) = delete;
  Record &operator =(const Record &) = default;

#if !defined(NDEBUG)
  void Print() const
  {
    printf("Record: {t=%f, M=%f, V=%f, h=%f, l=%f, Gamma=%f\n",
           (float)t, (float)M, (float)V, (float)h, (float)l, Gamma);
  }
#endif
};

// Context allows a thread to suspend and resume computations for the same meteorite
// Its size depends on the number of steps in Adams' method
// Use 'SizeOfThreadContext()' to determine the actual size
struct ThreadContext
{
  VirtualMeteoroid params{};
  size_t nxt{};
  real t{};
  bool ended{};
  Adams::State curr_layer;
  Adams::State steps[0];
};

constexpr size_t SizeOfThreadContext(size_t steps)
{
  assert(steps >= 1);
  assert(steps <= 3);
  return sizeof(ThreadContext) + steps * sizeof(Adams::State);
}


void BatchedAdamsKernel(ThreadContext *contexts, int32_t *active_threads,
                        const VirtualMeteoroid *problems, size_t n_problems,
                        size_t adams_steps, real dt, real timeout,
                        Record *records,
                        size_t iterations, size_t threads_per_block, cudaStream_t stream);
