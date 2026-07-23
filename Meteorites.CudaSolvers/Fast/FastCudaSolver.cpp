#include "FastCudaSolver.h"

#include "FastAdamsKernel.h"
#include "MeteoroidsManager.h"

FastCudaSolver::FastCudaSolver(FastCudaSolverConfig config)
  : config_(config)
{
  int device = 0;
  HANDLE_ERROR(cudaGetDevice(&device));
  HANDLE_ERROR(cudaGetDeviceProperties(&props_, device));
}

FastCudaSolver::~FastCudaSolver()
{
  try
  {
    context_points_.reset();
    reference_points_.reset();
    timestamps_.reset();
  }
  catch(std::exception &)
  {
  }
}

// How many meteorites must be simulated at one kernel call
size_t FastCudaSolver::BatchSize() const
{
  return config_.meteoroids_per_thread * config_.threads_per_block * config_.blocks_per_sm * props_.multiProcessorCount;
}

void FastCudaSolver::Configure(NumericalAlgorithm alg, real dt, real timeout)
{
  size_t steps{};
  switch (alg) {
    case NumericalAlgorithm::ONE_STEP_ADAMS:
      steps = 1;
    break;

    case NumericalAlgorithm::TWO_STEP_ADAMS:
      steps = 2;
    break;

    case NumericalAlgorithm::THREE_STEP_ADAMS:
      steps = 3;
    break;

    default:
      throw std::runtime_error("unsupported numerical algorithm");
  }

  if (dt <= (real)0.0 || timeout < dt * (steps + 1))
  {
    throw std::runtime_error("wrong time step and/or timeout");
  }
  algorithm_   = alg;
  steps_       = steps;
  dt_          = dt;
  timeout_     = timeout;
}

std::vector<std::pair<VirtualMeteoroid, real> >
FastCudaSolver::Solve(const IMeteorite &meteorite, size_t n_meteoroids, size_t m_bests)
{
  assert(n_meteoroids >= m_bests);

  // Make n_meteoroids a multiple of BatchSize to prevent additional checks in kernels
  if (n_meteoroids % BatchSize())
  { n_meteoroids += BatchSize() - n_meteoroids % BatchSize(); }
  printf("\nTotal meteoroids: %.2lf mln\t\tBatch size: %.2lf mln\n\n",
      (double)n_meteoroids / 1000000, (double)BatchSize() / 1000000);

  uint32_t total_threads = config_.threads_per_block * config_.blocks_per_sm * props_.multiProcessorCount;

  // Prepare device buffers
  MeteoroidsManager meteoroids_manager(m_bests, total_threads, config_.best_meteoroids_per_thread);
  // Looks not clear :(
  // But if we move all buffers to MeteoroidsManager, it will become FastCudaSolver
  TrajectoryPoint ref_point0;
  {
  size_t n_timestamps;
  const real *timestamps, *v, *h;
  meteorite.Trajectory(n_timestamps, timestamps, v, h);
  ref_point0 = {v[0], h[0]};
  TrajectoryPoint reference_points[n_timestamps];
  for (size_t i = 0; i < n_timestamps; i++) { reference_points[i] = {v[i], h[i]}; }

  n_timestamps_ = n_timestamps;

  HANDLE_ERROR(CudaAlloc(context_points_, n_timestamps_ * total_threads));

  HANDLE_ERROR(CudaAlloc(reference_points_, n_timestamps_));
  HANDLE_ERROR(cudaMemcpy(reference_points_.get(), reference_points,
                          sizeof(TrajectoryPoint) * n_timestamps_, cudaMemcpyHostToDevice));

  HANDLE_ERROR(CudaAlloc(timestamps_, n_timestamps_));
  HANDLE_ERROR(cudaMemcpy(timestamps_.get(), timestamps,
                          sizeof(real) * n_timestamps_, cudaMemcpyHostToDevice));
  }

  // Batchs cicle
  while (n_meteoroids > 0)
  {
    meteoroids_manager.GenerateSeeds();

    // Launch kernel
    cudaEvent_t start, stop;
    float elapsedTime;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start, 0);    

    switch (Algorithm())
    {
      case NumericalAlgorithm::ONE_STEP_ADAMS:
        FastAdamsKernel<1u, FastCudaSolverConfig::best_meteoroids_per_thread>
            (meteoroids_manager.GetDeviceSeeds(), timestamps_.get(), n_timestamps_,
             reference_points_.get(), context_points_.get(), dt_, timeout_,
             meteoroids_manager.GetDeviceDeviationsBuffer(),
             meteoroids_manager.GetActualBorderDeviation(),
             config_.meteoroids_per_thread,
             config_.blocks_per_sm * props_.multiProcessorCount, config_.threads_per_block);
        break;

      case NumericalAlgorithm::TWO_STEP_ADAMS:
        FastAdamsKernel<2u, FastCudaSolverConfig::best_meteoroids_per_thread>
            (meteoroids_manager.GetDeviceSeeds(), timestamps_.get(), n_timestamps_,
             reference_points_.get(), context_points_.get(), dt_, timeout_,
             meteoroids_manager.GetDeviceDeviationsBuffer(),
             meteoroids_manager.GetActualBorderDeviation(),
             config_.meteoroids_per_thread,
             config_.blocks_per_sm * props_.multiProcessorCount, config_.threads_per_block);
        break;

      case NumericalAlgorithm::THREE_STEP_ADAMS:
        FastAdamsKernel<3u, FastCudaSolverConfig::best_meteoroids_per_thread>
            (meteoroids_manager.GetDeviceSeeds(), timestamps_.get(), n_timestamps_,
             reference_points_.get(), context_points_.get(), dt_, timeout_,
             meteoroids_manager.GetDeviceDeviationsBuffer(),
             meteoroids_manager.GetActualBorderDeviation(),
             config_.meteoroids_per_thread,
             config_.blocks_per_sm * props_.multiProcessorCount, config_.threads_per_block);
        break;

      default:
        assert(false);
    }

    cudaEventRecord(stop,0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&elapsedTime, start, stop);
    printf("Batch time by GPU: %.2lf sec\nSpeed:%.2lf mln/sec\n",
        (double)elapsedTime / 1000, ((double)BatchSize() / 1000000) / ((double)elapsedTime / 1000));
    cudaEventDestroy(start);
    cudaEventDestroy(stop);

    meteoroids_manager.UpdateTopMeteoroids();

    n_meteoroids -= BatchSize();

    printf("%.2lf mln meteoroids left. Achieved border deviation: %f\n\n",
        (double)n_meteoroids / 1000000, meteoroids_manager.GetActualBorderDeviation());
  }

  return meteoroids_manager.GetTopMeteoroids(ref_point0);
}