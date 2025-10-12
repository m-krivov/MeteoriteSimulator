#include "CudaSolver.h"

#include "BatchedAdamsKernel.h"

//--------------------
//--- BatchedAdams ---
//--------------------

namespace
{

// Simulation results are stored in some slice-based format
// This helper allows us to traverse them
class Enumerator
{
  public:
    Enumerator() = delete;
    Enumerator(const Enumerator &) = delete;
    Enumerator &operator =(const Enumerator &) = delete;

    Enumerator(const std::vector<std::unique_ptr<Record[]>> &blocks,
               size_t meteorites_per_block, size_t iterations_per_block,
               size_t meteorite_idx)
      : iterations_per_block_(iterations_per_block), meteorite_idx_(meteorite_idx),
        block_(0), iteration_(0), finished_(false), blocks_(blocks)
    {
      assert(!blocks.empty());
      assert(meteorites_per_block > 0);
      assert(iterations_per_block > 0);
      assert(meteorite_idx < meteorites_per_block);
    }

    // Moves to the next record, assign it to 'value' and returns true
    // If no records left, does not update 'value' and returns false
    bool MoveNext(const Record *&value)
    {
      if (finished_)
      { return false; }

      if (++iteration_ > iterations_per_block_)
      {
        iteration_ = 0;
        block_ += 1;
        assert(block_ < blocks_.size());
      }

      const Record *res = blocks_[block_].get() + iteration_ + meteorite_idx_ * iterations_per_block_;
      if (res->t == 0.0)
      {
        finished_ = true;
        return false;
      }

      value = res;
      return true;
    }

  private:
    const size_t iterations_per_block_{}, meteorite_idx_{};
    size_t block_{}, iteration_{};
    bool finished_{};
    const std::vector<std::unique_ptr<Record[]>> &blocks_;
};


// Performs simulation for a batch of all cases
// Expects that all buffers points to device-accessible memory and have valid sizes
template <uint32_t STEPS>
void BatchedAdams(const std::vector<Case> &problems, real dt, real timeout,
                  const IFunctional &functional, IResultFormatter &results,

                  size_t batch_size, size_t iterations_per_batch, size_t threads_per_block,

                  uint32_t *dev_active_meteorites, Case *dev_problems,
                  ThreadContext<STEPS> *dev_contexts, Record *dev_records,
                  real *dev_timestamps, real *dev_functional_args)
{
  assert(dt > (real)0.0);
  assert(timeout > dt * STEPS);
  assert(threads_per_block > 0);
  assert(iterations_per_batch > 0);

  // Cache timestamps, they are the same for all meteorites
  const real *timestamps = nullptr;
  size_t n_timestamps = 0;
  functional.GetTimeStamps(n_timestamps, timestamps);
  assert(n_timestamps > 0);
  HANDLE_ERROR(cudaMemcpy(dev_timestamps, timestamps,
                          sizeof(real) * n_timestamps, cudaMemcpyHostToDevice));

  size_t functional_args_size = n_timestamps * 2 * batch_size;
  auto functional_args = std::make_unique<real[]>(functional_args_size);

  // Split all problems into a set of batches
  size_t n_batches = (problems.size() + batch_size - 1) / batch_size;
  for (size_t batch = 0; batch < n_batches; batch++)
  {
    // Re-initialize buffers
    size_t n_meteorites = std::min(batch_size, problems.size() - batch * batch_size);
    HANDLE_ERROR(cudaMemset(dev_contexts, 0, sizeof(ThreadContext<STEPS>) * batch_size));
    HANDLE_ERROR(cudaMemset(dev_functional_args, 0, sizeof(real) * functional_args_size));

    uint32_t active_meteorites = n_meteorites;
    HANDLE_ERROR(cudaMemcpy(dev_active_meteorites, &active_meteorites,
                            sizeof(uint32_t), cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMemcpy(dev_problems, problems.data() + batch * batch_size,
                            sizeof(Case) * n_meteorites, cudaMemcpyHostToDevice));

    // Perform Adams' iterations while at least one thread of the batch is active
    std::vector<std::unique_ptr<Record[]>> records;
    while (active_meteorites > 0)
    {
      BatchedAdamsKernel<STEPS>(dev_contexts, dev_active_meteorites,
                                dev_problems, n_meteorites, dt, timeout,
                                dev_timestamps, n_timestamps,
                                dev_functional_args, dev_records,
                                iterations_per_batch, threads_per_block);

      records.emplace_back(std::make_unique<Record[]>(iterations_per_batch * batch_size));
      HANDLE_ERROR(cudaMemcpy(records[records.size() - 1].get(), dev_records,
                              sizeof(Record) * iterations_per_batch * batch_size, cudaMemcpyDeviceToHost));

      HANDLE_ERROR(cudaMemcpy(&active_meteorites, dev_active_meteorites,
                              sizeof(active_meteorites), cudaMemcpyDeviceToHost));
    }

    // Register the simulation results
    HANDLE_ERROR(cudaMemcpy(functional_args.get(), dev_functional_args,
                            sizeof(real) * functional_args_size, cudaMemcpyDeviceToHost));

    for (size_t meteorite = 0; meteorite < n_meteorites; ++meteorite)
    {
      // Submit trajectory to the formatter
      auto t_next = results.Started(problems[batch * batch_size + meteorite]);
    
      Enumerator en(records, n_meteorites, iterations_per_batch, meteorite);
      const Record *record{};
      while (en.MoveNext(record))
      {
        if (record->t >= t_next)
        {
          t_next = results.Store(record->t, record->M, record->V, record->h,
                                  record->l, record->Gamma);
        }
      }

      // Compute value of the functional for this meteorite
      real *V_args = functional_args.get() + (meteorite * n_timestamps * 2);
      real *h_args = V_args + n_timestamps;
      size_t timestamp = 0;
      while (timestamp < n_timestamps && V_args[timestamp] != (real)0.0)
      { timestamp++; };
      auto f_val = functional.Compute(timestamp, V_args, h_args);

      // Finalize the meteorite
      // We can continue using 'record' as it is still valid
      results.Finished(IResultFormatter::Classify(record->t, record->M, record->h), f_val);
    }
  }
}

} // unnamed namespace


//------------------
//--- CudaSolver ---
//------------------

CudaSolver::CudaSolver(Config config)
  : config_(config)
{
  int device = 0;
  HANDLE_ERROR(cudaGetDevice(&device));
  HANDLE_ERROR(cudaGetDeviceProperties(&props_, device));

  size_t batch_size = BatchSize();
  HANDLE_ERROR(CudaAlloc(buffer_counter_, sizeof(uint32_t)));
  HANDLE_ERROR(CudaAlloc(buffer_problems_, batch_size * sizeof(Case)));
  HANDLE_ERROR(CudaAlloc(buffer_contexts_, batch_size * sizeof(ThreadContext<3>)));
  HANDLE_ERROR(CudaAlloc(buffer_records_,  batch_size * config_.iterations_per_block * sizeof(Record)));
}

CudaSolver::~CudaSolver()
{
  try
  {
    buffer_counter_.reset();
    buffer_problems_.reset();
    buffer_contexts_.reset();
    buffer_records_.reset();
    buffer_timestamps_.reset();
    buffer_functional_.reset();
  }
  catch (std::exception &)
  {
    // Sadly, but we can do nothing
  }
}

size_t CudaSolver::BatchSize() const
{
  return config_.threads_per_block * config_.blocks_per_sm * props_.multiProcessorCount;
}

void CudaSolver::Solve(ICaseGenerator &generator, const IFunctional &functional, IResultFormatter &results)
{
  Case problem;
  std::vector<Case> problems;
  size_t batch_size = BatchSize();

  problems.reserve(batch_size);
  do
  {
    problems.clear();
    while (problems.size() < batch_size && generator.Next(problem))
    { problems.emplace_back(std::move(problem)); }

    if (!problems.empty())
    { Solve(problems, functional, results); }
  }
  while (!problems.empty());
}

void CudaSolver::Solve(const std::vector<Case> &problems, const IFunctional &functional, IResultFormatter &results)
{
  // Resize buffers for functional arguments and timestamps (if needed)
  size_t n_timestamps{};
  {
    const real *_{};
    functional.GetTimeStamps(n_timestamps, _);
  }

  if (n_timestamps_ < n_timestamps)
  {
    n_timestamps_ = n_timestamps;
    HANDLE_ERROR(CudaAlloc(buffer_timestamps_, n_timestamps_ * sizeof(real)));
    HANDLE_ERROR(CudaAlloc(buffer_functional_, n_timestamps_ * 2 * BatchSize() * sizeof(real)));
  }

  // Perform simulations
  switch (Algorithm())
  {
    case NumericalAlgorithm::ONE_STEP_ADAMS:
      BatchedAdams<1>(problems, Dt(), Timeout(), functional, results,
                      BatchSize(), config_.iterations_per_block, config_.threads_per_block,
                      (uint32_t *)buffer_counter_.get(), (Case *)buffer_problems_.get(),
                      (ThreadContext<1> *)buffer_contexts_.get(), (Record *)buffer_records_.get(),
                      (real *)buffer_timestamps_.get(), (real *)buffer_functional_.get());
      break;

    case NumericalAlgorithm::TWO_STEP_ADAMS:
      BatchedAdams<2>(problems, Dt(), Timeout(), functional, results,
                      BatchSize(), config_.iterations_per_block, config_.threads_per_block,
                      (uint32_t *)buffer_counter_.get(), (Case *)buffer_problems_.get(),
                      (ThreadContext<2> *)buffer_contexts_.get(), (Record *)buffer_records_.get(),
                      (real *)buffer_timestamps_.get(), (real *)buffer_functional_.get());
      break;

    case NumericalAlgorithm::THREE_STEP_ADAMS:
      BatchedAdams<3>(problems, Dt(), Timeout(), functional, results,
                      BatchSize(), config_.iterations_per_block, config_.threads_per_block,
                      (uint32_t *)buffer_counter_.get(), (Case *)buffer_problems_.get(),
                      (ThreadContext<3> *)buffer_contexts_.get(), (Record *)buffer_records_.get(),
                      (real *)buffer_timestamps_.get(), (real *)buffer_functional_.get());
      break;

    default:
      assert(false);
  }
}

void CudaSolver::Solve(const Case &problem, const IFunctional &functional, IResultFormatter &results)
{
  std::vector<Case> problems = { problem };
  Solve(problems, functional, results);
}
