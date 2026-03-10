#include "IterationAllocator.h"
#include "PedanticCudaSolver.h"

#include "BatchedAdamsKernel.h"


//-----------------------------------------
//--- PedanticCudaSolver::DeviceContext ---
//-----------------------------------------

// All allocated CUDA resources are stored in a single place
struct PedanticCudaSolver::DeviceContext
{
  // Buffers were allocated for batches of this size
  size_t batch_size{};

  std::array<cudaStream_t, 2> streams;
  std::array<cudaEvent_t, 2> copy_events;
  std::array<cudaEvent_t, 2> kernel_events;

  std::array<uint8_t *, 2> pinned_buffers;
  
  CudaPtr<int32_t> counter;
  CudaPtr<VirtualMeteoroid> meteoroids;
  CudaPtr<uint8_t> thread_contexts;
  CudaPtr<Record> records;
};


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
void BatchedAdams(BasicSolver::MeteoroidEnumerator &problems, real dt, real timeout,
                  const IFunctional &functional, ISimulationRecorder &results,

                  const PedanticCudaSolver::DeviceContext &context,
                  size_t adams_steps, size_t iterations_per_batch, size_t threads_per_block)
{
  assert(dt > (real)0.0);
  assert(timeout > dt * adams_steps);
  assert(threads_per_block > 0);
  assert(iterations_per_batch > 0);

  IterationAllocator allocator;
  std::future<void> futures[2];
  
  auto memcpy_async = [](void *dst, const void *src, size_t size) {
    return std::async(std::launch::async, [=]() {
      std::memcpy(dst, src, size);
    });
  };

  const auto dev_active_meteorites = context.counter.get();
  const auto dev_meteoroids = context.meteoroids.get();
  const auto dev_thread_contexts = (ThreadContext *)context.thread_contexts.get();
  const auto dev_records = context.records.get();
  const auto pinned_records = context.pinned_buffers;

  // Process all virtual meteoroids by batches
  const VirtualMeteoroid *meteoroids = nullptr;
  size_t n_meteoroids = 0;
  while (problems.MoveNext(meteoroids, n_meteoroids))
  {
    assert(n_meteoroids <= context.batch_size);
    size_t half_size = n_meteoroids * iterations_per_batch;
    std::vector<UniquePtr> records;

    // Re-initialize buffers
    HANDLE_ERROR(cudaMemset(dev_thread_contexts, 0,
                            SizeOfThreadContext(adams_steps) * n_meteoroids));

    int32_t active_meteorites = n_meteoroids;
    HANDLE_ERROR(cudaMemcpy(dev_active_meteorites, &active_meteorites,
                            sizeof(int32_t), cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMemcpy(dev_meteoroids, meteoroids,
                            sizeof(VirtualMeteoroid) * n_meteoroids, cudaMemcpyHostToDevice));

    // The first step of pipeline
    int iter = 0;

    BatchedAdamsKernel(dev_thread_contexts, dev_active_meteorites,
                       dev_meteoroids, n_meteoroids, adams_steps, dt, timeout,
                       dev_records + iter * half_size,
                       iterations_per_batch, threads_per_block, context.streams[iter]);
    HANDLE_ERROR(cudaEventRecord(context.kernel_events[iter], context.streams[iter]));
    HANDLE_ERROR(cudaMemcpyAsync(pinned_records[iter],
                                 dev_records + iter * half_size,
                                 half_size * sizeof(Record),
                                 cudaMemcpyDeviceToHost, context.streams[iter]));
    HANDLE_ERROR(cudaEventRecord(context.copy_events[iter], context.streams[iter]));

    HANDLE_ERROR(cudaEventRecord(context.copy_events[!iter], context.streams[!iter]));
    futures[!iter] = std::async([](){});

    while (active_meteorites > 0)
    {
      iter = !iter;

      HANDLE_ERROR(cudaStreamWaitEvent(context.streams[iter], context.copy_events[iter]));
      HANDLE_ERROR(cudaStreamWaitEvent(context.streams[iter], context.kernel_events[!iter]));
      BatchedAdamsKernel(dev_thread_contexts, dev_active_meteorites,
                         dev_meteoroids, n_meteoroids, adams_steps, dt, timeout,
                         dev_records + iter * half_size,
                         iterations_per_batch, threads_per_block, context.streams[iter]);
      HANDLE_ERROR(cudaEventRecord(context.kernel_events[iter], context.streams[iter]));

      futures[iter].wait();
      HANDLE_ERROR(cudaEventSynchronize(context.copy_events[!iter]));
      HANDLE_ERROR(cudaEventSynchronize(context.kernel_events[iter]));

      records.emplace_back(allocator.Alloc<Record>(half_size));
      futures[!iter] = memcpy_async(records.back().get(), pinned_records[!iter],
                                    half_size * sizeof(Record));
      HANDLE_ERROR(cudaMemcpy(&active_meteorites, dev_active_meteorites,
                              sizeof(int32_t), cudaMemcpyDeviceToHost));

      HANDLE_ERROR(cudaMemcpyAsync(pinned_records[iter],
                                   dev_records + iter * half_size,
                                   half_size * sizeof(Record),
                                   cudaMemcpyDeviceToHost, context.streams[iter]));
      HANDLE_ERROR(cudaEventRecord(context.copy_events[iter], context.streams[iter]));
    }

    // The final step of pipeline
    iter = !iter;

    futures[iter].wait();
    HANDLE_ERROR(cudaEventSynchronize(context.copy_events[!iter]));

    records.emplace_back(allocator.Alloc<Record>(half_size));
    futures[!iter] = memcpy_async(records.back().get(), pinned_records[!iter],
                                  half_size * sizeof(Record));
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

PedanticCudaSolver::PedanticCudaSolver(PedanticCudaConfig config)
  : config_(config), context_{ new DeviceContext() }
{
  int device = 0;
  HANDLE_ERROR(cudaGetDevice(&device));
  HANDLE_ERROR(cudaGetDeviceProperties(&props_, device));
  context_->batch_size = BatchSize();

  HANDLE_ERROR(cudaStreamCreate(&context_->streams[0]));
  HANDLE_ERROR(cudaStreamCreate(&context_->streams[1]));
  HANDLE_ERROR(cudaEventCreate(&context_->copy_events[0]));
  HANDLE_ERROR(cudaEventCreate(&context_->copy_events[1]));
  HANDLE_ERROR(cudaEventCreate(&context_->kernel_events[0]));
  HANDLE_ERROR(cudaEventCreate(&context_->kernel_events[1]));

  size_t half_size = context_->batch_size * config.iterations_per_block * sizeof(Record);
  HANDLE_ERROR(cudaMallocHost(&context_->pinned_buffers[0], half_size));
  HANDLE_ERROR(cudaMallocHost(&context_->pinned_buffers[1], half_size));

  HANDLE_ERROR(CudaAlloc(context_->counter, 1));
  HANDLE_ERROR(CudaAlloc(context_->meteoroids, context_->batch_size));
  HANDLE_ERROR(CudaAlloc(context_->thread_contexts,
                         context_->batch_size * SizeOfThreadContext(3 /* the actual N is unknown */)));

  // Double record buffer for two CUDA streams
  HANDLE_ERROR(CudaAlloc(context_->records, context_->batch_size * config.iterations_per_block * 2));
}

PedanticCudaSolver::~PedanticCudaSolver()
{
  try
  {
    context_->counter.reset();
    context_->meteoroids.reset();
    context_->thread_contexts.reset();
    context_->records.reset();

    cudaFreeHost(context_->pinned_buffers[0]);
    cudaFreeHost(context_->pinned_buffers[1]);

    cudaStreamDestroy(context_->streams[0]);
    cudaStreamDestroy(context_->streams[1]);
    cudaEventDestroy(context_->copy_events[0]);
    cudaEventDestroy(context_->copy_events[1]);
    cudaEventDestroy(context_->kernel_events[0]);
    cudaEventDestroy(context_->kernel_events[1]);
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
  BatchedAdams(problems, Dt(), Timeout(), functional, results,
               *context_, Steps(), config_.iterations_per_block, config_.threads_per_block);
}
