#pragma once
#include "Meteorites.CudaSolvers/CudaDefs.h"

#include "Meteorites.Core/Solvers/BasicSolver.h"


// Settings for CudaSolver
struct PedanticCudaSolverConfig
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

  PedanticCudaSolverConfig() = default;
  PedanticCudaSolverConfig(const PedanticCudaSolverConfig &) = default;
  PedanticCudaSolverConfig &operator =(const PedanticCudaSolverConfig &) = default;
};

// Version for CUDA kernel debugging
class PedanticCudaSolver : public BasicSolver
{
  public:
    PedanticCudaSolver(PedanticCudaSolverConfig config = PedanticCudaSolverConfig());
    virtual ~PedanticCudaSolver();

  private:
    // BasicSolver method
    // How many meteorites must be simulated at one time
    virtual size_t BatchSize() const override final;

    virtual void SolveAny(MeteoroidEnumerator &enumerator,
                          const IFunctional &functional,
                          ISimulationRecorder &results) override final;

    cudaDeviceProp props_{};
    const PedanticCudaSolverConfig config_{};

    std::array<cudaStream_t, 2> streams_;
    std::array<cudaEvent_t, 2> copy_events_;
    std::array<cudaEvent_t, 2> kernel_events_;
    std::array<uint8_t *, 2> pinned_buffers_;

    CudaPtr<uint8_t> buffer_counter_;
    CudaPtr<uint8_t> buffer_problems_;
    CudaPtr<uint8_t> buffer_contexts_;
    CudaPtr<uint8_t> buffer_records_;
    size_t n_timestamps_{0};    // implicitly defines the size of 'buffer_functional_' and 'buffer_timestamps_'
};
