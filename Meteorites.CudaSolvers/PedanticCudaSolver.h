#pragma once
#include "Meteorites.CudaSolvers/CudaDefs.h"

#include "Meteorites.Core/Solvers/BasicSolver.h"


// Settings for PedanticCudaSolver
struct PedanticCudaConfig
{
  // How many CUDA threads must be spawned per each block
  // Will be used to configure the grid for a CUDA kernel
  size_t threads_per_block = 32;

  // Each block performs no more than the requested number of Adams' iterations
  // After that, it sends the intermediate simulation results to host
  size_t iterations_per_block = 200;

  // How many blocks must be spawned per each CUDA streaming multiprocessor
  // Large numbers can lead to high memory usage
  size_t blocks_per_sm = 4;
};

// Version for debugging and testing
// Provides identical results but works faster than GoldSolver
class PedanticCudaSolver : public BasicSolver
{
  public:
    // Definition for pImpl
    struct DeviceContext;

    PedanticCudaSolver(PedanticCudaConfig config = PedanticCudaConfig());
    virtual ~PedanticCudaSolver();

  private:
    // BasicSolver method
    // How many meteorites must be simulated at one time
    virtual size_t BatchSize() const override final;

    virtual void SolveAny(MeteoroidEnumerator &enumerator,
                          const IFunctional &functional,
                          ISimulationRecorder &results) override final;

    cudaDeviceProp props_{};
    const PedanticCudaConfig config_{};
    std::shared_ptr<DeviceContext> context_;
};
