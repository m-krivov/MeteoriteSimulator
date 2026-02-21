#include "IterationAllocator.h"
#include "PedanticCudaSolver.h"

#include "BatchedAdamsKernel.h"

//--------------------
//--- BatchedAdams ---
//--------------------

namespace
{

using UniquePtr = std::unique_ptr<Record, IterationAllocator::Deleter>;

// Simulation results (trajectories) are stored in some slice-based format
// This helper allows us to traverse them and distinguish records that corresponds to timestamps
class TrajectoryEnumerator
{
  public:
    TrajectoryEnumerator() = delete;
    TrajectoryEnumerator(const TrajectoryEnumerator &) = delete;
    TrajectoryEnumerator &operator =(const TrajectoryEnumerator &) = delete;

    TrajectoryEnumerator(const real *timestamps, size_t n_timestamps,
                         const std::vector<UniquePtr> &blocks,
                         size_t meteorites_per_block, size_t iterations_per_block,
                         size_t meteorite_idx)
      : timestamps_(timestamps), n_timestamps_(n_timestamps), next_timestamp_{0},
        iterations_per_block_(iterations_per_block), meteorite_idx_(meteorite_idx),
        block_(0), iteration_(0), finished_(false), blocks_(blocks)
    {
      assert(!blocks.empty());
      assert(meteorites_per_block > 0);
      assert(iterations_per_block > 0);
      assert(meteorite_idx < meteorites_per_block);
    }

    // Moves to the next record, assign it to 'value' and returns true
    // If no records left, does not update 'value' and returns false
    // The 'timestamp' flag is set if valus of this record must be passed to the functional
    bool MoveNext(const Record *&value, bool &is_timestamp)
    {
      if (finished_)
      { return false; }

      if (iteration_ >= iterations_per_block_)
      {
        iteration_ = 0;
        block_ += 1;
        assert(block_ < blocks_.size());
      }

      is_timestamp = false;
      value = blocks_[block_].get() + iteration_ + meteorite_idx_ * iterations_per_block_;
      if (value->t < 0.0)
      {
        finished_ = true;
        return false;
      }

      if (next_timestamp_ < n_timestamps_ &&
          value->t >= timestamps_[next_timestamp_])
      {
        next_timestamp_ += 1;
        is_timestamp = true;
      }

      iteration_ += 1;
      return true;
    }

  private:
    const real *timestamps_{};
    size_t n_timestamps_{}, next_timestamp_{};
    const size_t iterations_per_block_{}, meteorite_idx_{};
    size_t block_{}, iteration_{};
    bool finished_{};
    const std::vector<UniquePtr> &blocks_;
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
      v.clear();
      h.clear();

      // Submit trajectory to the formatter, collect values for functional
      auto t_next = results.Started(meteoroids[meteoroid]);
    
      TrajectoryEnumerator en(timestamps, n_timestamps, records,
                              n_meteoroids, iterations_per_batch, meteoroid);
      const Record *record{};
      bool is_timestamp{};
      while (en.MoveNext(record, is_timestamp))
      {
        if (record->t >= t_next)
        {
          t_next = results.Store(record->t, record->M, record->V, record->h,
                                 record->l, record->Gamma);
        }
        if (is_timestamp)
        {
          v.emplace_back(record->V);
          h.emplace_back(record->h);
        }
      }

      // Compute value of the functional and finalize the meteoroid
      // We can continue using 'record' as it is still valid
      assert(v.size() <= n_timestamps);
      assert(h.size() == v.size());
      auto f_val = functional.Compute(v.size(), v.data(), h.data());
      auto reason = ISimulationRecorder::Classify(record->t, record->M, record->h);
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
