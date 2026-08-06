#pragma once

#include "Meteorites.CudaSolvers/CudaDefs.h"
#include "Meteorites.Core/IMeteorite.h"
#include "Meteorites.Core/Solvers/ISolver.h"

// Settings for FastCudaSolver
struct FastCudaSolverConfig
{
  uint32_t meteoroids_per_thread = 4096 * 4;
  
  // How many CUDA threads must be spawned per each block
  // Will be used to configure the grid for the CUDA kernel
  uint32_t threads_per_block = 64;

  // How many blocks must be spawned per each CUDA streaming multiprocessor
  // Large numbers can lead to high memory usage
  uint32_t blocks_per_sm = 8;

  // Size of the local array of best meteoroids per thread
  // MUST be constexpr static to be used as template parameter
  constexpr static uint32_t best_meteoroids_per_thread = 4;

  FastCudaSolverConfig() = default;
  FastCudaSolverConfig(const FastCudaSolverConfig &) = default;
  FastCudaSolverConfig &operator =(const FastCudaSolverConfig &) = default;
};

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
    NumericalAlgorithm algorithm_ = NumericalAlgorithm::TWO_STEP_ADAMS;
    size_t             steps_ = 2;
    real               dt_ = (real)0.001;
    real               timeout_ = (real)1000.0;

    cudaDeviceProp             props_{};
    const FastCudaSolverConfig config_{};

    CudaPtr<TrajectoryPoint>    context_points_;
    CudaPtr<TrajectoryPoint>    reference_points_;
    CudaPtr<real>               timestamps_;
    uint32_t                    n_timestamps_{0};
};