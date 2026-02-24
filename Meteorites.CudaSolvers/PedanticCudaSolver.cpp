#include "IterationAllocator.h"
#include "PedanticCudaSolver.h"

#include "BatchedAdamsKernel.h"

//--------------------
//--- BatchedAdams ---
//--------------------

namespace
{

using UniquePtr = std::unique_ptr<Record, IterationAllocator::Deleter>;

// Simulation results (trajectories) are stored in some internal slice-based format
// This helper allows us to traverse them
class TrajectoryEnumerator
{
  public:
    TrajectoryEnumerator() = delete;
    TrajectoryEnumerator(const TrajectoryEnumerator &) = delete;
    TrajectoryEnumerator &operator =(const TrajectoryEnumerator &) = delete;

    TrajectoryEnumerator(const std::vector<UniquePtr> &blocks,
                         size_t meteorites_per_block, size_t iterations_per_block,
                         size_t meteorite_idx)
      : blocks_(blocks), iterations_per_block_(iterations_per_block), meteorite_idx_(meteorite_idx),
        iteration_{-1}
    {
      assert(!blocks.empty());
      assert(meteorites_per_block > 0);
      assert(iterations_per_block > 0);
      assert(meteorite_idx < meteorites_per_block);

      dt_ = (blocks_[0].get() + 1)->t;
      assert(blocks_[0]->t == (real)0.0);
      assert(dt_ > (real)0.0);
    }

    // Returns a pointer to the current trajectory point if the enumerator is valid
    // Or nullptr otherwise
    const Record *Current() const
    { return current_; }

    // Moves to the next record and returns true (or false, if no records left)
    // Must be called before the first use of 'Current()'
    bool MoveNext()
    {
      if (finished_)
      { return false; }
      
      iteration_ += 1;
      assert(iteration_ >= 0);
      auto block = iteration_ / iterations_per_block_;
      auto offset = iteration_ % iterations_per_block_;
      if (block >= blocks_.size() ||
          (current_ = blocks_[block].get() + offset + meteorite_idx_ * iterations_per_block_)->t < 0.0)
      {
        finished_ = true;
        current_ = nullptr;
        return false;
      }

      return true;
    }

    // Tries to move to the first record made after time 't'
    // If no such record exists, returns false
    bool MoveTo(real t)
    {
      t = std::max((real)0, t);
      iteration_ = (int64_t)(t / dt_) - 1;
      current_ = nullptr;
      finished_ = false;

      // To avoid accuracy issues, let's traverse one excessive record
      iteration_ = std::max((int64_t)-1, iteration_ - 1);
      while (MoveNext()) {
        if (Current()->t >= t) {
          return true;
        }
      }
      return false;
    }

    // Moves to the last record
    // Since at least one record is always exist, this method does not return status
    void MoveLast()
    {
      // We need to traverse the last block to detect the stop marker
      iteration_ = (blocks_.size() - 1) * iterations_per_block_ - 1;
      current_ = nullptr;
      finished_ = false;
      while (MoveNext()) {}

      iteration_ = std::max((int64_t)-1, iteration_ - 2);
      current_ = nullptr;
      finished_ = false;
      bool status = MoveNext();
      assert(status);
    }

  private:
    const std::vector<UniquePtr> &blocks_;
    const size_t iterations_per_block_{}, meteorite_idx_{};
    real dt_{};
    int64_t iteration_{};
    const Record *current_{};
    bool finished_{};
};


// Performs simulation for a batch of all virtual meteoroids
// Expects that all buffers points to device-accessible memory and have valid sizes
template <uint32_t STEPS>
void BatchedAdams(BasicSolver::MeteoroidEnumerator &problems, real dt, real timeout,
                  const IFunctional &functional, ISimulationRecorder &results,

                  size_t iterations_per_batch, size_t threads_per_block,

                  int32_t *dev_active_meteorites, VirtualMeteoroid *dev_problems,
                  ThreadContext<STEPS> *dev_contexts, Record *dev_records,
                
                  const std::array<cudaStream_t, 2> &streams,
                  const std::array<cudaEvent_t, 2> &copy_events,
                  const std::array<cudaEvent_t, 2> &kernel_events,
                  const std::array<uint8_t *, 2> &pinned_buffers)
{
  assert(dt > (real)0.0);
  assert(timeout > dt * STEPS);
  assert(threads_per_block > 0);
  assert(iterations_per_batch > 0);

  IterationAllocator allocator;
  std::future<void> futures[2];
  
  auto memcpy_async = [](void *dst, const void *src, size_t size) {
    return std::async(std::launch::async, [=]() {
      std::memcpy(dst, src, size);
    });
  };

  // Process all virtual meteoroids by batches
  const VirtualMeteoroid *meteoroids = nullptr;
  size_t n_meteoroids = 0;
  while (problems.MoveNext(meteoroids, n_meteoroids))
  {
    size_t half_size = n_meteoroids * iterations_per_batch;
    std::vector<UniquePtr> records;

    // Re-initialize buffers
    HANDLE_ERROR(cudaMemset(dev_contexts, 0, sizeof(ThreadContext<STEPS>) * n_meteoroids));

    int32_t active_meteorites = n_meteoroids;
    HANDLE_ERROR(cudaMemcpy(dev_active_meteorites, &active_meteorites,
                            sizeof(int32_t), cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMemcpy(dev_problems, meteoroids,
                            sizeof(VirtualMeteoroid) * n_meteoroids, cudaMemcpyHostToDevice));

    // The first step of pipeline
    int iter = 0;

    BatchedAdamsKernel<STEPS>(dev_contexts, dev_active_meteorites,
                              dev_problems, n_meteoroids, dt, timeout,
                              dev_records + iter * half_size,
                              iterations_per_batch, threads_per_block, streams[iter]);
    HANDLE_ERROR(cudaEventRecord(kernel_events[iter], streams[iter]));
    HANDLE_ERROR(cudaMemcpyAsync(pinned_buffers[iter],
                                 dev_records + iter * half_size,
                                 half_size * sizeof(Record),
                                 cudaMemcpyDeviceToHost, streams[iter]));
    HANDLE_ERROR(cudaEventRecord(copy_events[iter], streams[iter]));

    HANDLE_ERROR(cudaEventRecord(copy_events[!iter], streams[!iter]));
    futures[!iter] = std::async([](){});

    while (active_meteorites > 0)
    {
      iter = !iter;

      HANDLE_ERROR(cudaStreamWaitEvent(streams[iter], copy_events[iter]));
      HANDLE_ERROR(cudaStreamWaitEvent(streams[iter], kernel_events[!iter]));
      BatchedAdamsKernel<STEPS>(dev_contexts, dev_active_meteorites,
                                dev_problems, n_meteoroids, dt, timeout,
                                dev_records + iter * half_size,
                                iterations_per_batch, threads_per_block, streams[iter]);
      HANDLE_ERROR(cudaEventRecord(kernel_events[iter], streams[iter]));

      futures[iter].wait();
      HANDLE_ERROR(cudaEventSynchronize(copy_events[!iter]));
      HANDLE_ERROR(cudaEventSynchronize(kernel_events[iter]));

      records.emplace_back(allocator.Alloc<Record>(half_size));
      futures[!iter] = memcpy_async(records.back().get(), pinned_buffers[!iter], half_size * sizeof(Record));
      HANDLE_ERROR(cudaMemcpy(&active_meteorites, dev_active_meteorites,
                              sizeof(int32_t), cudaMemcpyDeviceToHost));

      HANDLE_ERROR(cudaMemcpyAsync(pinned_buffers[iter],
                                   dev_records + iter * half_size,
                                   half_size * sizeof(Record),
                                   cudaMemcpyDeviceToHost, streams[iter]));
      HANDLE_ERROR(cudaEventRecord(copy_events[iter], streams[iter]));
    }

    // The final step of pipeline
    iter = !iter;

    futures[iter].wait();
    HANDLE_ERROR(cudaEventSynchronize(copy_events[!iter]));

    records.emplace_back(allocator.Alloc<Record>(half_size));
    futures[!iter] = memcpy_async(records.back().get(), pinned_buffers[!iter], half_size * sizeof(Record));
    futures[!iter].wait();

    HANDLE_ERROR(cudaDeviceSynchronize());

    // Extract timestamps, they are the same for all meteorites
    const real *timestamps = nullptr;
    size_t n_timestamps = 0;
    functional.GetTimeStamps(n_timestamps, timestamps);
    assert(n_timestamps > 0);

    std::vector<real> v, h;
    v.reserve(n_timestamps);
    h.reserve(n_timestamps);

    for (size_t meteoroid = 0; meteoroid < n_meteoroids; ++meteoroid)
    {
      // Submit trajectory to the formatter
      auto t_next = results.Started(meteoroids[meteoroid]);
    
      TrajectoryEnumerator en(records, n_meteoroids, iterations_per_batch, meteoroid);
      while (en.MoveTo(t_next))
      {
        auto record = en.Current();
        assert(record != nullptr);
        t_next = results.Store(record->t, record->M, record->V, record->h,
                               record->l, record->Gamma);
      }

      // Collect values for the functional and compute it
      v.clear();
      h.clear();
      for (size_t i = 0; i < n_timestamps; i++)
      {
        if (!en.MoveTo(timestamps[i]))
        { break; }

        v.emplace_back(en.Current()->V);
        h.emplace_back(en.Current()->h);
      }
      assert(v.size() <= n_timestamps);
      assert(h.size() == v.size());
      auto f_val = functional.Compute(v.size(), v.data(), h.data());
      
      // Finalize the meteoroid
      en.MoveLast();
      auto reason = ISimulationRecorder::Classify(en.Current()->t, en.Current()->M, en.Current()->h);
      results.Finished(reason, f_val);
    }
  }
}

} // unnamed namespace


//------------------
//--- CudaSolver ---
//------------------

PedanticCudaSolver::PedanticCudaSolver(PedanticCudaSolverConfig config)
  : config_(config)
{
  int device = 0;
  HANDLE_ERROR(cudaGetDevice(&device));
  HANDLE_ERROR(cudaGetDeviceProperties(&props_, device));
  size_t batch_size = BatchSize();

  HANDLE_ERROR(cudaStreamCreate(&streams_[0]));
  HANDLE_ERROR(cudaStreamCreate(&streams_[1]));
  HANDLE_ERROR(cudaEventCreate(&copy_events_[0]));
  HANDLE_ERROR(cudaEventCreate(&copy_events_[1]));
  HANDLE_ERROR(cudaEventCreate(&kernel_events_[0]));
  HANDLE_ERROR(cudaEventCreate(&kernel_events_[1]));

  size_t half_size = batch_size * config.iterations_per_block * sizeof(Record);
  HANDLE_ERROR(cudaMallocHost(&pinned_buffers_[0], half_size));
  HANDLE_ERROR(cudaMallocHost(&pinned_buffers_[1], half_size));

  HANDLE_ERROR(CudaAlloc(buffer_counter_, sizeof(int32_t)));
  HANDLE_ERROR(CudaAlloc(buffer_problems_, batch_size * sizeof(VirtualMeteoroid)));
  HANDLE_ERROR(CudaAlloc(buffer_contexts_, batch_size * sizeof(ThreadContext<3>)));

  // Double record buffer for 2 CUDA streams
  HANDLE_ERROR(CudaAlloc(buffer_records_, half_size * 2));
}

PedanticCudaSolver::~PedanticCudaSolver()
{
  try
  {
    buffer_counter_.reset();
    buffer_problems_.reset();
    buffer_contexts_.reset();
    buffer_records_.reset();

    cudaFreeHost(pinned_buffers_[0]);
    cudaFreeHost(pinned_buffers_[1]);

    cudaStreamDestroy(streams_[0]);
    cudaStreamDestroy(streams_[1]);
    cudaEventDestroy(copy_events_[0]);
    cudaEventDestroy(copy_events_[1]);
    cudaEventDestroy(kernel_events_[0]);
    cudaEventDestroy(kernel_events_[1]);
  }
  catch (std::exception &)
  {
    // Sadly, but we can do nothing
  }
}

size_t PedanticCudaSolver::BatchSize() const
{
  return config_.threads_per_block * config_.blocks_per_sm * props_.multiProcessorCount;
}

void PedanticCudaSolver::SolveAny(MeteoroidEnumerator &problems,
                                  const IFunctional &functional,
                                  ISimulationRecorder &results)
{
  // Perform simulations
  auto method = [&](auto contexts)
  {
    BatchedAdams(problems, Dt(), Timeout(), functional, results,
                 config_.iterations_per_block, config_.threads_per_block,
                 (int32_t *)buffer_counter_.get(), (VirtualMeteoroid *)buffer_problems_.get(),
                 contexts, (Record *)buffer_records_.get(),
                 streams_, copy_events_, kernel_events_, pinned_buffers_);
  };

  switch (Algorithm())
  {
    case NumericalAlgorithm::ONE_STEP_ADAMS:
      method((ThreadContext<1> *)buffer_contexts_.get());
      break;

    case NumericalAlgorithm::TWO_STEP_ADAMS:
      method((ThreadContext<2> *)buffer_contexts_.get());
      break;

    case NumericalAlgorithm::THREE_STEP_ADAMS:
      method((ThreadContext<3> *)buffer_contexts_.get());
      break;

    default:
      assert(false);
  }
}
