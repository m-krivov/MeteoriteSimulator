#include "FastCudaSolver.h"

#include "FastAdamsKernel.h"
#include "Random.h"

void GenerateSeeds(std::vector<uint64_t> &seeds) {

  Random::State rootState;
  Random::Initialize(rootState);

  std::vector<Random::State> tempStates(seeds.size());

  Random::Multiply(rootState, tempStates);

  for (size_t i = 0; i < seeds.size(); i++)
  {
    uint32_t part1 = Random::Next(tempStates[i]);
    uint32_t part2 = Random::Next(tempStates[i]);
    seeds[i] = (static_cast<uint64_t>(part1) << 32) | part2;
  }
}

FastCudaSolver::FastCudaSolver(FastCudaSolverConfig config)
  : config_(config)
{
  int device = 0;
  HANDLE_ERROR(cudaGetDevice(&device));
  HANDLE_ERROR(cudaGetDeviceProperties(&props_, device));

  seeds_.resize(BatchSize());
}

FastCudaSolver::~FastCudaSolver()
{
  try
  {
    seeds_.clear();
    functional_args_.clear();
    functional_args_references_.clear();
    best_meteoroids_deviations_.clear();
    best_meteoroids_curand_states_.clear();
    timestamps_.clear();
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

  // Make n_meteoroids a multiple of meteoroids_per_thread to prevent additional checks in kernels
  if (n_meteoroids % config_.meteoroids_per_thread)
  { n_meteoroids += config_.meteoroids_per_thread - n_meteoroids % config_.meteoroids_per_thread; }

  size_t threads_num = config_.threads_per_block * config_.blocks_per_sm * props_.multiProcessorCount;

  {
    size_t n_timestamps;
    const real *timestamps, *v, *h;
    meteorite.Trajectory(n_timestamps, timestamps, v, h);
    if (n_timestamps_ < n_timestamps)
    {
      n_timestamps_ = n_timestamps;
      timestamps_.resize(n_timestamps_);
      functional_args_.resize(n_timestamps_ * 2 * threads_num);
      functional_args_references_.resize(n_timestamps_ * 2);

      thrust::copy(timestamps, timestamps + n_timestamps, timestamps_.begin());
      thrust::copy(v, v + n_timestamps, functional_args_references_.begin());
      thrust::copy(h, h + n_timestamps, functional_args_references_.begin() + n_timestamps_);
    }
    // Top m_bests meteoroids are placed in memory before and next to local thread buffers
    // to simplify sorting between butch runs
    best_meteoroids_deviations_.resize(m_bests + config_.best_meteoroids_buffer_size * threads_num);
    best_meteoroids_curand_states_.resize(m_bests + config_.best_meteoroids_buffer_size * threads_num);
    thrust::fill(thrust::cuda::par,
                 best_meteoroids_deviations_.begin(),
                 best_meteoroids_deviations_.begin() + m_bests,
                 std::numeric_limits<real>::max());
  }

  std::vector<uint64_t> seeds(BatchSize());
  real border_deviation = std::numeric_limits<real>::max();

  while (n_meteoroids > 0)
  {
    if (n_meteoroids < BatchSize())
    {
      seeds.resize(n_meteoroids);
      seeds_.resize(n_meteoroids);
    }
    GenerateSeeds(seeds);
    thrust::copy(seeds.begin(), seeds.end(), seeds_.begin());

    switch (Algorithm())
    {
      case NumericalAlgorithm::ONE_STEP_ADAMS:
        FastAdamsKernel<1u, FastCudaSolverConfig::best_meteoroids_buffer_size>
            (thrust::raw_pointer_cast(seeds_.data()), seeds_.size(),
             dt_, timeout_,
             thrust::raw_pointer_cast(timestamps_.data()), n_timestamps_,
             thrust::raw_pointer_cast(functional_args_.data()),
             thrust::raw_pointer_cast(functional_args_references_.data()),
             thrust::raw_pointer_cast(best_meteoroids_deviations_.data()) + m_bests,
             thrust::raw_pointer_cast(best_meteoroids_curand_states_.data()) + m_bests,
             border_deviation,
             config_.meteoroids_per_thread,
             config_.blocks_per_sm * props_.multiProcessorCount, config_.threads_per_block);
        break;

      case NumericalAlgorithm::TWO_STEP_ADAMS:
        FastAdamsKernel<2u, FastCudaSolverConfig::best_meteoroids_buffer_size>
            (thrust::raw_pointer_cast(seeds_.data()), seeds_.size(),
             dt_, timeout_,
             thrust::raw_pointer_cast(timestamps_.data()), n_timestamps_,
             thrust::raw_pointer_cast(functional_args_.data()),
             thrust::raw_pointer_cast(functional_args_references_.data()),
             thrust::raw_pointer_cast(best_meteoroids_deviations_.data()) + m_bests,
             thrust::raw_pointer_cast(best_meteoroids_curand_states_.data()) + m_bests,
             border_deviation,
             config_.meteoroids_per_thread,
             config_.blocks_per_sm * props_.multiProcessorCount, config_.threads_per_block);
        break;

      case NumericalAlgorithm::THREE_STEP_ADAMS:
        FastAdamsKernel<3u, FastCudaSolverConfig::best_meteoroids_buffer_size>
            (thrust::raw_pointer_cast(seeds_.data()), seeds_.size(),
             dt_, timeout_,
             thrust::raw_pointer_cast(timestamps_.data()), n_timestamps_,
             thrust::raw_pointer_cast(functional_args_.data()),
             thrust::raw_pointer_cast(functional_args_references_.data()),
             thrust::raw_pointer_cast(best_meteoroids_deviations_.data()) + m_bests,
             thrust::raw_pointer_cast(best_meteoroids_curand_states_.data()) + m_bests,
             border_deviation,
             config_.meteoroids_per_thread,
             config_.blocks_per_sm * props_.multiProcessorCount, config_.threads_per_block);
        break;

      default:
        assert(false);
    }

    HANDLE_ERROR(cudaDeviceSynchronize());

    thrust::sort_by_key(thrust::cuda::par,
                        best_meteoroids_deviations_.begin(),
                        best_meteoroids_deviations_.end(),
                        best_meteoroids_curand_states_.begin());

    border_deviation = best_meteoroids_deviations_[m_bests - 1];
    n_meteoroids -= seeds.size();
    /*test*/printf("%zu meteoroids left. Current border_deviation: %f\n", n_meteoroids, border_deviation);
  }
  best_meteoroids_deviations_.resize(m_bests);
  best_meteoroids_curand_states_.resize(m_bests);

  thrust::device_vector<VirtualMeteoroid> best_meteoroids_(m_bests);
  RestoreMeteoroidsKernel(thrust::raw_pointer_cast(best_meteoroids_.data()),
                          thrust::raw_pointer_cast(best_meteoroids_curand_states_.data()),
                          m_bests,
                          functional_args_references_[0],
                          functional_args_references_[n_timestamps_]);
  HANDLE_ERROR(cudaDeviceSynchronize());

  std::vector<VirtualMeteoroid> best_meteoroids(m_bests);
  std::vector<real> best_deviations(m_bests);
  thrust::copy(best_meteoroids_.begin(), best_meteoroids_.end(), best_meteoroids.begin());
  thrust::copy(best_meteoroids_deviations_.begin(), best_meteoroids_deviations_.end(), best_deviations.begin());
  std::vector<std::pair<VirtualMeteoroid, real> > result(m_bests);
  for (size_t i = 0; i < m_bests; i++)
  {
    result[i] = {best_meteoroids[i], best_deviations[i]};
  }
  return result;
}