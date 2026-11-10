#pragma once

#include "Meteorites.CudaSolvers/CudaDefs.h"

class FlopsEstimator
{
  public:
    FlopsEstimator() = delete;
    FlopsEstimator(const FlopsEstimator &) = delete;
    FlopsEstimator &operator =(const FlopsEstimator &) = delete;

    FlopsEstimator(size_t flops_per_step, size_t total_threads);
    ~FlopsEstimator();

    uint64_t *GetStepsCounters();
    double    CountTflops(double seconds, size_t absoluteTflopsTolerance);
  
  private:
    CudaPtr<uint64_t> steps_counters_;
    size_t total_threads_;
    size_t flops_per_step_;
};