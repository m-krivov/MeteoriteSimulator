#pragma once

#include "Meteorites.CudaSolvers/CudaDefs.h"
#include "Meteorites.Core/Meteoroids/VirtualMeteoroid.h"

// Simple hardcoded generator for use in CUDA kernels
class CudaMeteoroidsGenerator
{
  public:
    __device__ CudaMeteoroidsGenerator(uint64_t seed, uint64_t curand_offset = 0)
      : current_offset_(curand_offset)
    {
      curand_init(seed, 0, curand_offset, &current_state_);
    }
    CudaMeteoroidsGenerator() = delete;
    CudaMeteoroidsGenerator(const CudaMeteoroidsGenerator &) = delete;
    CudaMeteoroidsGenerator &operator =(const CudaMeteoroidsGenerator &) = delete;

    __device__ void Next(VirtualMeteoroid &meteoroid, const TrajectoryPoint &ref_poin0)
    {
      meteoroid.H = 1e5 + curand_uniform(&current_state_) * (5e6 - 1e5);
      meteoroid.Ch = 0.1 + curand_uniform(&current_state_) * (0.9 - 0.1);
      meteoroid.Rho = 2000.0 + curand_uniform(&current_state_) * (5000.0 - 2000.0);
      meteoroid.Cd = 0.5 + curand_uniform(&current_state_) * (2.5 - 0.5);
      meteoroid.Cl = 0.0 + curand_uniform(&current_state_) * (0.25 - 0.0);
      meteoroid.M0 = 10.0 + curand_uniform(&current_state_) * (500.0 - 10.0);
      meteoroid.V0 = ref_poin0.v;
      meteoroid.h0 = ref_poin0.h;
      meteoroid.Gamma0 = curand_uniform(&current_state_) * (M_PI / 2);

      current_offset_ += number_of_parameters;
    }

    // Returns the actual curand_offset for the last generated meteoroid
    __device__ uint64_t GetOffset()
    {
      assert(current_offset_ != 0); // No one meteoroid had been generated yet

      return current_offset_ - number_of_parameters;
    }

  private:
    curandState current_state_;
    // Raw offset for curand_init. This does NOT represent the number of generated meteoroids
    uint64_t    current_offset_;

    constexpr static uint64_t number_of_parameters = 7;
};