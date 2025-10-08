#include "CudaManager.h"

#include "AdamsKernels.h"
#include "GPUParameters.h"


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


template <uint32_t STEPS>
void CudaManager(const std::vector<Case> &problems, real dt, real timeout, size_t block_size,
                 const IFunctional &functional, IResultFormatter &results)
{
  assert(block_size > 0);
  assert(dt > (real)0.0);
  assert(timeout > dt * 3);

  // Allocate GPU buffers
  CudaPtr<ThreadContext<STEPS>> dev_contexts;
  HANDLE_ERROR(CudaAlloc(dev_contexts, block_size));
  CudaPtr<Case> dev_problems;
  HANDLE_ERROR(CudaAlloc(dev_problems, block_size));

  const real *timestamps = nullptr;
  size_t n_timestamps = 0;
  functional.GetTimeStamps(n_timestamps, timestamps);
  CudaPtr<real> dev_timestamps;
  HANDLE_ERROR(CudaAlloc(dev_timestamps, n_timestamps));
  HANDLE_ERROR(cudaMemcpy(dev_timestamps.get(), timestamps,
                          sizeof(real) * n_timestamps, cudaMemcpyHostToDevice));

  size_t functional_args_size = n_timestamps * 2 * block_size;
  auto functional_args = std::make_unique<real[]>(functional_args_size);
  CudaPtr<real> dev_functional_args;
  HANDLE_ERROR(CudaAlloc(dev_functional_args, functional_args_size));

  CudaPtr<Record> dev_records;
  CudaAlloc(dev_records, ITERS_PER_KERNEL * block_size);

  CudaPtr<uint32_t> dev_active_meteorites;
  HANDLE_ERROR(CudaAlloc(dev_active_meteorites, 1));

  // Split all problems into a set of blocks
  size_t n_blocks = (problems.size() + block_size - 1) / block_size;
  for (size_t block = 0; block < n_blocks; block++)
  {
    size_t n_meteorites = std::min(block_size, problems.size() - block * block_size);
    HANDLE_ERROR(cudaMemset(dev_contexts.get(), 0, sizeof(ThreadContext<STEPS>) * block_size));
    HANDLE_ERROR(cudaMemset(dev_functional_args.get(), 0, functional_args_size));

    uint32_t active_meteorites = n_meteorites;
    HANDLE_ERROR(cudaMemcpy(dev_active_meteorites.get(), &active_meteorites,
                            sizeof(uint32_t), cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMemcpy(dev_problems.get(), problems.data() + block * block_size,
                            sizeof(Case) * n_meteorites, cudaMemcpyHostToDevice));

    // Perform Adams' iterations while at least one thread of the block is active
    std::vector<std::unique_ptr<Record[]>> records;
    while (active_meteorites > 0)
    {
      BatchedAdamsKernel<STEPS>(dev_contexts.get(), dev_active_meteorites.get(),
                                dev_problems.get(), n_meteorites, dt, timeout,
                                dev_timestamps.get(), n_timestamps,
                                dev_functional_args.get(), dev_records.get(), ITERS_PER_KERNEL);

      records.emplace_back(std::make_unique<Record[]>(ITERS_PER_KERNEL * block_size));
      HANDLE_ERROR(cudaMemcpy(records[records.size() - 1].get(), dev_records.get(),
                              sizeof(Record) * ITERS_PER_KERNEL * block_size, cudaMemcpyDeviceToHost));

      HANDLE_ERROR(cudaMemcpy(&active_meteorites, dev_active_meteorites.get(),
                              sizeof(active_meteorites), cudaMemcpyDeviceToHost));
    }

    // Register the simulation results
    HANDLE_ERROR(cudaMemcpy(functional_args.get(), dev_functional_args.get(),
                            sizeof(real) * functional_args_size, cudaMemcpyDeviceToHost));

    for (size_t meteorite = 0; meteorite < n_meteorites; ++meteorite)
    {
      // Submit trajectory to the formatter
      auto t_next = results.Started(problems[block * block_size + meteorite]);
    
      Enumerator en(records, problems.size(), ITERS_PER_KERNEL, meteorite);
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

// template-specified functions compiles only this way
template void CudaManager<1u>(const std::vector<Case>& problems, real dt, real timeout, size_t block_size,
                              const IFunctional& functional, IResultFormatter& results);
template void CudaManager<2u>(const std::vector<Case>& problems, real dt, real timeout, size_t block_size,
                              const IFunctional& functional, IResultFormatter& results);
template void CudaManager<3u>(const std::vector<Case>& problems, real dt, real timeout, size_t block_size,
                              const IFunctional& functional, IResultFormatter& results);
