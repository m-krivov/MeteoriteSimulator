#pragma once

#include "Meteorites.CudaSolvers/CudaDefs.h"
#include "Meteorites.Core/IMeteorite.h"
#include "Meteorites.Core/Solvers/ISolver.h"

// Settings for FastCudaSolver
struct FastCudaSolverConfig
{
  constexpr static size_t meteoroids_per_thread = 64;
  
  // How many CUDA threads must be spawned per each block
  // Will be used to configure the grid for a CUDA kernel
  constexpr static size_t threads_per_block = 32;

  // How many blocks must be spawned per each CUDA streaming multiprocessor
  // Large numbers can lead to high memory usage
  constexpr static size_t blocks_per_sm = 4;

  // Size of local arrays of the best meteorites for each thread
  // MUST be constexpr static for using as template parameter
  constexpr static size_t best_meteoroids_buffer_size = 4;

  FastCudaSolverConfig() = default;
  FastCudaSolverConfig(const FastCudaSolverConfig &) = default;
  FastCudaSolverConfig &operator =(const FastCudaSolverConfig &) = default;
};

// Version performing calculations
class FastCudaSolver
{
  public:
    FastCudaSolver(FastCudaSolverConfig config = FastCudaSolverConfig());
    ~FastCudaSolver();

    size_t BatchSize() const;
    NumericalAlgorithm Algorithm() const { return algorithm_; }

    void Configure(NumericalAlgorithm alg, real dt, real timeout);

    std::vector<std::pair<VirtualMeteoroid, real> >
    Solve(const IMeteorite &meteorite, size_t n_meteoroids, size_t m_bests);

  private:
    NumericalAlgorithm algorithm_ = NumericalAlgorithm::ONE_STEP_ADAMS;
    size_t steps_ = 1;
    real dt_ = (real)0.001;
    real timeout_ = (real)1000.0;

    cudaDeviceProp props_{};
    const FastCudaSolverConfig config_{};

    thrust::device_vector<uint64_t> seeds_;
    thrust::device_vector<real> functional_args_;
    thrust::device_vector<real> functional_args_references_;
    thrust::device_vector<real> best_meteoroids_deviations_;
    thrust::device_vector<curandState> best_meteoroids_curand_states_;
    thrust::device_vector<real> timestamps_;
    size_t n_timestamps_{0};
};