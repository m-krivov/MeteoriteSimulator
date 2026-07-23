#pragma once
#include "Meteorites.CudaSolvers/CudaDefs.h"

#include "Meteorites.Core/Meteoroids/VirtualMeteoroid.h"
#include "Meteorites.Core/Solvers/Adams.h"


// Contains a few values extracted from Adams::Layer
// They could be passed to formatters
struct alignas(32) Record
{
  real t{}, M{}, V{}, h{}, l{}, Gamma{};

  Record() = default;
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
template <uint32_t STEPS>
struct ThreadContext
{
  Adams::Layer steps[STEPS + 1];
  VirtualMeteoroid params{};
  size_t nxt{};
  real t{};
  size_t timestamp{};
  bool ended{};

  ThreadContext() = default;
  ThreadContext(const ThreadContext<STEPS> &) = delete;
  ThreadContext<STEPS> &operator =(const ThreadContext<STEPS> &) = delete;
};


template <unsigned int STEPS>
void BatchedAdamsKernel(ThreadContext<STEPS> *contexts, int32_t *active_threads,
                        const VirtualMeteoroid *problems, size_t n_problems, real dt, real timeout,
                        const real *timestamps, size_t n_timestamps,
                        real *functional_args, Record *records,
                        size_t iterations, size_t threads_per_block, cudaStream_t stream);