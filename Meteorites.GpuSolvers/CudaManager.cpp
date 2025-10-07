#include "CudaManager.h"

#include "AdamsKernels.h"
#include "GPUParameters.h"


template <uint32_t STEPS>
void CudaManager(const std::vector<Case> &problems, const IFunctional &functional,
                 real dt, real timeout, IResultFormatter &results)
{
  size_t n_threads = problems.size() / CASES_PER_THREAD;

  // Prepare data for GPU
  CudaPtr<ThreadContext<STEPS>> contexts;
  HANDLE_ERROR(CudaAlloc(contexts, n_threads));
  CudaPtr<Case> dev_problems;
  HANDLE_ERROR(CudaAlloc(dev_problems, problems.size()));
  HANDLE_ERROR(cudaMemcpy(dev_problems.get(), problems.data(),
                          sizeof(Case) * problems.size(), cudaMemcpyHostToDevice));

  const real *timestamps = nullptr;
  size_t n_timestamps = 0;
  functional.GetTimeStamps(n_timestamps, timestamps);
  CudaPtr<real> dev_timestamps;
  HANDLE_ERROR(CudaAlloc(dev_timestamps, n_timestamps));
  HANDLE_ERROR(cudaMemcpy(dev_timestamps.get(), timestamps,
                          sizeof(real) * n_timestamps, cudaMemcpyHostToDevice));

  size_t functional_args_size = n_timestamps * 2 * n_threads * CASES_PER_THREAD;
  auto functional_args = std::make_unique<real[]>(functional_args_size);
  CudaPtr<real> dev_functional_args;
  HANDLE_ERROR(CudaAlloc(dev_functional_args, functional_args_size));

  std::vector<std::unique_ptr<Record[]>> records;
  CudaPtr<Record> dev_records;
  CudaAlloc(dev_records, ITERS_PER_KERNEL * n_threads);

  uint32_t active_threads = n_threads;
  CudaPtr<uint32_t> dev_active_threads;
  HANDLE_ERROR(CudaAlloc(dev_active_threads, 1));
  HANDLE_ERROR(cudaMemcpy(dev_active_threads.get(), &active_threads,
                          sizeof(uint32_t), cudaMemcpyHostToDevice));

  // Perform Adams' iterations while at least one thread is active
  int global_iter = 0;
  while (active_threads) {
    global_iter++;

    BatchedAdamsKernel<STEPS>(contexts.get(), dev_active_threads.get(),
                              dev_problems.get(), problems.size(), dt, timeout,
                              dev_timestamps.get(), n_timestamps,
                              dev_functional_args.get(), dev_records.get(), ITERS_PER_KERNEL);

    HANDLE_ERROR(cudaGetLastError());

    records.push_back(std::make_unique<Record[]>(ITERS_PER_KERNEL * n_threads));
    HANDLE_ERROR(cudaMemcpy(records[records.size() - 1].get(), dev_records.get(),
                 sizeof(Record) * ITERS_PER_KERNEL * n_threads, cudaMemcpyDeviceToHost));

    HANDLE_ERROR(cudaMemcpy(&active_threads, dev_active_threads.get(),
                            sizeof(active_threads), cudaMemcpyDeviceToHost));
  }

  // Register the simulation results
  HANDLE_ERROR(cudaMemcpy(functional_args.get(), dev_functional_args.get(),
                          sizeof(real) * functional_args_size, cudaMemcpyDeviceToHost));

  size_t problem_idx = 0;
  for (size_t thread_idx = 0; thread_idx < n_threads; ++thread_idx) {
    size_t records_idx = 0, rec_i = 0; // global and local indices
    for (size_t case_idx = 0; case_idx < CASES_PER_THREAD; ++case_idx) {

      // Prepare new case before writing
      auto t_next = results.Started(problems[problem_idx]);
      Record* rec_block = records[records_idx].get() + thread_idx * ITERS_PER_KERNEL;
      for (uint32_t i = 0; i < STEPS + 1; i++) {
        t_next = results.Store(rec_block[rec_i].t, rec_block[rec_i].M, rec_block[rec_i].V, rec_block[rec_i].h,
                               rec_block[rec_i].l, rec_block[rec_i].Gamma);
        rec_i++;
      }

      // Write case's necessary records to formatter and functional
      while (true) {

        // Go to next records block if the current one is ended
        if (rec_i == ITERS_PER_KERNEL) {
          rec_i = 0;
          records_idx++;
          assert(records_idx != records.size());
          rec_block = records[records_idx].get() + thread_idx * ITERS_PER_KERNEL;
        }

        // End of case's records
        if (rec_block[rec_i].t == (real)0.0) {
          real* V_args =
              functional_args.get() + (thread_idx * n_timestamps * 2 * CASES_PER_THREAD) + (n_timestamps * 2 * case_idx);
          real* h_args = V_args + n_timestamps;
          size_t timestamp = 0;
          while (timestamp < n_timestamps && V_args[timestamp] != (real)0.0) {
            timestamp++;
          };
          if (rec_block[rec_i].M <= (real)0.01) {
            results.Finished(IResultFormatter::Reason::Burnt, functional.Compute(timestamp, V_args, h_args));
          } else if (rec_block[rec_i].h <= (real)0.0) {
            results.Finished(IResultFormatter::Reason::Collided, functional.Compute(timestamp, V_args, h_args));
          } else {
            results.Finished(IResultFormatter::Reason::Timeouted, functional.Compute(timestamp, V_args, h_args));
          }
          rec_i = 0; // new case's records starts from next kernel launch
          records_idx++;
          break;
        }

        // Store necessary record
        if (rec_block[rec_i].t >= t_next) {
          t_next = results.Store(rec_block[rec_i].t, rec_block[rec_i].M, rec_block[rec_i].V, rec_block[rec_i].h,
                                 rec_block[rec_i].l, rec_block[rec_i].Gamma);
        }

        rec_i++;
      }
      problem_idx++;
    }
  }
}

// template-specified functions compiles only this way
template void CudaManager<1u>(const std::vector<Case>& problems, const IFunctional& functional,
                              real dt, real timeout, IResultFormatter& results);
template void CudaManager<2u>(const std::vector<Case>& problems, const IFunctional& functional,
                              real dt, real timeout, IResultFormatter& results);
template void CudaManager<3u>(const std::vector<Case>& problems, const IFunctional& functional,
                              real dt, real timeout, IResultFormatter& results);
