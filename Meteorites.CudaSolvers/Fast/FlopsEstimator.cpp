#include "FlopsEstimator.h"

#include <iostream>

FlopsEstimator::FlopsEstimator(size_t flops_per_step, size_t total_threads)
  : flops_per_step_(flops_per_step), total_threads_(total_threads)
{
  HANDLE_ERROR(CudaAlloc(steps_counters_, total_threads));
}

FlopsEstimator::~FlopsEstimator()
{
  try
  {
    steps_counters_.reset();
  }
  catch(std::exception &)
  {
  }
}

uint64_t *FlopsEstimator::GetStepsCounters()
{
  return steps_counters_.get();
}

double FlopsEstimator::CountTflops(double seconds, size_t absoluteTflopsTolerance)
{
  std::vector<uint64_t> host_steps_counters(total_threads_);
  HANDLE_ERROR(cudaMemcpy(host_steps_counters.data(), steps_counters_.get(),
                          sizeof(uint64_t) * total_threads_, cudaMemcpyDeviceToHost));

  constexpr size_t teraPow = 12;
  int safeDropDigits = teraPow - absoluteTflopsTolerance - (std::to_string(total_threads_).length() + 1);
  if (safeDropDigits < 0)
  {
    std::cerr << "Very high precision. Reducing it by " << std::abs(safeDropDigits) << "..." << std::endl;
    safeDropDigits = 0;
  }

  uint64_t steps_sum = 0;
  for (size_t i = 0; i < host_steps_counters.size(); i++)
  {
    if (host_steps_counters[i] / std::pow(10, safeDropDigits) > std::numeric_limits<uint64_t>::max() - steps_sum) {
      std::cerr << "TFLOPS counter has overflowed. Try to reduce precision." << std::endl;
    }
    steps_sum += host_steps_counters[i] / (uint64_t)std::pow(10, safeDropDigits);
  }

  return (double)steps_sum / std::pow(10, teraPow - safeDropDigits) * (double)flops_per_step_ / seconds;
}