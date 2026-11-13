#include "FastCudaSolver.h"

#include "FastAdamsKernel.h"
#include "MeteoroidsManager.h"
#include "Meteorites.Core/Solvers/Adams.h"

FastCudaSolver::FastCudaSolver(FastCudaSolverConfig config)
  : config_(config)
{
  int device = 0;
  HANDLE_ERROR(cudaGetDevice(&device));
  HANDLE_ERROR(cudaGetDeviceProperties(&props_, device));
}

// Number of meteoroids simulated per kernel launch
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

namespace
{

// Estimates the performance of FastCudaSolver in floating point operations per second (FLOPS).
// To make these estimates more accurate, configure the project with the '-DMETEORITES_FAST_MATH=OFF' option.
double EstimateFlops(NumericalAlgorithm algorithm, const std::vector<uint64_t> &iterations, double seconds)
{
  if (seconds <= 0.0)
  { return 0.0; }

  OperationCounter operations_per_iteration;
  switch (algorithm)
  {
    case NumericalAlgorithm::ONE_STEP_ADAMS:
      operations_per_iteration = Adams::OneStepIteration();
      break;

    case NumericalAlgorithm::TWO_STEP_ADAMS:
      operations_per_iteration = Adams::TwoStepIteration();
      break;

    case NumericalAlgorithm::THREE_STEP_ADAMS:
      operations_per_iteration = Adams::ThreeStepIteration();
      break;

    default:
      throw std::runtime_error("unsupported numerical algorithm");
  }

  long double flops_per_iteration = operations_per_iteration.arithmetic +
                                    operations_per_iteration.division +
                                    operations_per_iteration.special;

  // Compensated summation reduces precision loss when adding many large values
  long double total_iterations = 0.0L;
  long double compensation = 0.0L;
  for (const uint64_t steps : iterations)
  {
    const long double value = (long double)steps - compensation;
    const long double next_total = total_iterations + value;
    compensation = (next_total - total_iterations) - value;
    total_iterations = next_total;
  }

  const long double total_flops = total_iterations * flops_per_iteration;
  return (double)(total_flops / (long double)seconds);
}

} // unnamed namespace

std::vector<std::pair<VirtualMeteoroid, real> >
FastCudaSolver::Solve(const IMeteorite &meteorite, size_t n_meteoroids, size_t m_bests, double *flops)
{
  assert(n_meteoroids >= m_bests);

  // Make n_meteoroids a multiple of BatchSize to prevent additional checks in kernels
  if (n_meteoroids % BatchSize())
  { n_meteoroids += BatchSize() - n_meteoroids % BatchSize(); }
  printf("\nTotal meteoroids: %.2lf mln\t\tBatch size: %.2lf mln\n\n",
      (double)n_meteoroids / 1000000, (double)BatchSize() / 1000000);

  size_t total_threads = config_.threads_per_block * config_.blocks_per_sm * props_.multiProcessorCount;

  // Prepare device buffers
  MeteoroidsManager meteoroids_manager(m_bests, total_threads, config_.best_meteoroids_per_thread);
  // Looks unclear :(
  // However, moving all buffers into MeteoroidsManager would turn it into FastCudaSolver
  TrajectoryPoint ref_point0;
  size_t n_timestamps;
  CudaPtr<uint64_t>           iteration_counters_dev;
  CudaPtr<TrajectoryPoint>    context_points_dev;
  CudaPtr<TrajectoryPoint>    reference_points_dev;
  CudaPtr<real>               timestamps_dev;
  {
    const real *timestamps, *v, *h;
    meteorite.Trajectory(n_timestamps, timestamps, v, h);
    ref_point0 = {v[0], h[0]};
    std::vector<TrajectoryPoint> reference_points(n_timestamps);
    for (size_t i = 0; i < n_timestamps; i++) { reference_points[i] = {v[i], h[i]}; }

    if (flops != nullptr)
    { HANDLE_ERROR(CudaAlloc(iteration_counters_dev, total_threads)); }
    HANDLE_ERROR(CudaAlloc(context_points_dev, n_timestamps * total_threads));

    HANDLE_ERROR(CudaAlloc(reference_points_dev, n_timestamps));
    HANDLE_ERROR(cudaMemcpy(reference_points_dev.get(), reference_points.data(),
                            sizeof(TrajectoryPoint) * n_timestamps, cudaMemcpyHostToDevice));

    HANDLE_ERROR(CudaAlloc(timestamps_dev, n_timestamps));
    HANDLE_ERROR(cudaMemcpy(timestamps_dev.get(), timestamps,
                            sizeof(real) * n_timestamps, cudaMemcpyHostToDevice));
  }

  // Batches loop
  while (n_meteoroids > 0)
  {
    meteoroids_manager.GenerateSeeds();

    // Launch kernel
    cudaEvent_t start, stop;
    float elapsedTime;
    HANDLE_ERROR(cudaEventCreate(&start));
    HANDLE_ERROR(cudaEventCreate(&stop));
    HANDLE_ERROR(cudaEventRecord(start, 0));    

    switch (Algorithm())
    {
      case NumericalAlgorithm::ONE_STEP_ADAMS:
        FastAdamsKernel<1u, FastCudaSolverConfig::best_meteoroids_per_thread>
            (meteoroids_manager.GetDeviceSeeds(), timestamps_dev.get(), n_timestamps,
             reference_points_dev.get(), context_points_dev.get(), dt_, timeout_,
             meteoroids_manager.GetDeviceDeviationsBuffer(),
             meteoroids_manager.GetActualThresholdDeviation(),
             config_.meteoroids_per_thread,
             config_.blocks_per_sm * props_.multiProcessorCount, config_.threads_per_block,
             iteration_counters_dev.get());
        break;

      case NumericalAlgorithm::TWO_STEP_ADAMS:
        FastAdamsKernel<2u, FastCudaSolverConfig::best_meteoroids_per_thread>
            (meteoroids_manager.GetDeviceSeeds(), timestamps_dev.get(), n_timestamps,
             reference_points_dev.get(), context_points_dev.get(), dt_, timeout_,
             meteoroids_manager.GetDeviceDeviationsBuffer(),
             meteoroids_manager.GetActualThresholdDeviation(),
             config_.meteoroids_per_thread,
             config_.blocks_per_sm * props_.multiProcessorCount, config_.threads_per_block,
             iteration_counters_dev.get());
        break;

      case NumericalAlgorithm::THREE_STEP_ADAMS:
        FastAdamsKernel<3u, FastCudaSolverConfig::best_meteoroids_per_thread>
            (meteoroids_manager.GetDeviceSeeds(), timestamps_dev.get(), n_timestamps,
             reference_points_dev.get(), context_points_dev.get(), dt_, timeout_,
             meteoroids_manager.GetDeviceDeviationsBuffer(),
             meteoroids_manager.GetActualThresholdDeviation(),
             config_.meteoroids_per_thread,
             config_.blocks_per_sm * props_.multiProcessorCount, config_.threads_per_block,
             iteration_counters_dev.get());
        break;

      default:
        assert(false);
    }

    // TODO: pretty statistical output
    HANDLE_ERROR(cudaEventRecord(stop,0));
    HANDLE_ERROR(cudaEventSynchronize(stop));
    HANDLE_ERROR(cudaEventElapsedTime(&elapsedTime, start, stop));
    printf("Batch time by GPU: %.2lf sec\nSpeed: %.2lf mln/sec\n",
        (double)elapsedTime / 1000, ((double)BatchSize() / 1000000) / ((double)elapsedTime / 1000));

    if (iteration_counters_dev != nullptr && flops != nullptr)
    {
      std::vector<uint64_t> iteration_counters_host(total_threads);
      HANDLE_ERROR(cudaMemcpy(iteration_counters_host.data(), iteration_counters_dev.get(),
                              sizeof(uint64_t) * total_threads, cudaMemcpyDeviceToHost));
      *flops = EstimateFlops(algorithm_, iteration_counters_host, elapsedTime * 1e-3);
      printf("Achieved TFLOPS: %.3lf\n", std::round(*flops * 1e-9) * 1e-3);
    }

    HANDLE_ERROR(cudaEventDestroy(start));
    HANDLE_ERROR(cudaEventDestroy(stop));

    meteoroids_manager.UpdateTopMeteoroids();

    n_meteoroids -= BatchSize();

    printf("%.2lf mln meteoroids left. Achieved threshold deviation: %f\n\n",
        (double)n_meteoroids / 1000000, meteoroids_manager.GetActualThresholdDeviation());
  }

  return meteoroids_manager.GetTopMeteoroids(ref_point0);
}