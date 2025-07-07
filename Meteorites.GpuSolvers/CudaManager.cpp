#include "CudaManager.h"

#include "CudaAdamsMethod.h"
#include "CudaStructures.h"
#include "GPUParameters.h"
#include "CudaHandleError.h"
#include "CudaUniquePtr.h"

#include <memory>

template <unsigned int STEPS, unsigned int ITERS>
void CudaManager(const std::vector<Case>& problems_vector, const IFunctional& functional, real dt, real timeout,
                 IResultFormatter& results)
{
  // configure grids (make number of cases the main parameter)

  cudaDeviceProp prop;
  HANDLE_ERROR(cudaGetDeviceProperties(&prop, 0));
  unsigned int launching_threads_per_block = THREADS_PER_BLOCK;
  unsigned int cases_per_thread = CASES_PER_THREAD;
  unsigned int launching_blocks = ((problems_vector.size() - 1) / cases_per_thread / launching_threads_per_block) + 1;
  unsigned int count_of_threads = problems_vector.size() / CASES_PER_THREAD; // effective threads
  printf("\nGPU: %s\nConfig:\n>Blocks: %u\n>Threads per block: %u\n>Cases per thread: %u\nTotal threads:%u\nTotal "
         "cases:%lu\n",
         prop.name, launching_blocks, launching_threads_per_block, cases_per_thread,
         launching_threads_per_block * launching_blocks, problems_vector.size());

  // preparing buffers

  // thread contexts
  auto dev_thread_sandbox_arr = cuda_make_unique<ThreadContext<STEPS, ITERS>>(count_of_threads);

  // problems (read-only)
  auto dev_problems = cuda_make_unique<CudaCase>(problems_vector.size());
  HANDLE_ERROR(
      cudaMemcpy(dev_problems.get(), problems_vector.data(), sizeof(CudaCase) * problems_vector.size(), cudaMemcpyHostToDevice)); // damn reinterpret cast

  // meteorite's timestamps (read-only)
  const real* timestamps = nullptr;
  size_t n_timestamps = 0;
  functional.GetTimeStamps(n_timestamps, timestamps);
  auto dev_timestamps = cuda_make_unique<real>(n_timestamps);
  HANDLE_ERROR(cudaMemcpy(dev_timestamps.get(), timestamps, sizeof(real) * n_timestamps, cudaMemcpyHostToDevice));

  // (V_arg[] and h_arg[] for functional)'s for each case
  size_t functional_args_size = n_timestamps * 2 * count_of_threads * CASES_PER_THREAD;
  auto functional_args = std::make_unique<real[]>(functional_args_size);
  auto dev_functional_args = cuda_make_unique<real>(functional_args_size);

  // vector of Record's arrays
  // load to CPU vector every iteration
  std::vector<std::unique_ptr<Record[]>> records;
  auto dev_records = cuda_make_unique<Record>(ITERS_PER_KERNEL * count_of_threads);

  // active threads atomic counter (make reduction (or no))
  unsigned int active_threads = count_of_threads;
  auto dev_active_threads = cuda_make_unique<unsigned int>(1);
  HANDLE_ERROR(cudaMemcpy(dev_active_threads.get(), &active_threads, sizeof(unsigned int), cudaMemcpyHostToDevice));

  int global_iter = 0;
  while (active_threads) {
    global_iter++;

    CudaLauncher<STEPS, ITERS>(dev_thread_sandbox_arr.get(), dev_problems.get(), problems_vector.size(), dt, timeout,
                               dev_timestamps.get(), n_timestamps, dev_functional_args.get(), dev_records.get(), dev_active_threads.get(),
                               launching_blocks, launching_threads_per_block);

    cudaDeviceSynchronize();

    records.push_back(std::make_unique<Record[]>(ITERS_PER_KERNEL * count_of_threads));
    HANDLE_ERROR(cudaMemcpy(records[records.size() - 1].get(), dev_records.get(), sizeof(Record) * ITERS_PER_KERNEL * count_of_threads, cudaMemcpyDeviceToHost));

    HANDLE_ERROR(cudaMemcpy(&active_threads, dev_active_threads.get(), sizeof(active_threads), cudaMemcpyDeviceToHost));
  }

  // process records
  printf("process records...\n");

  HANDLE_ERROR(cudaMemcpy(functional_args.get(), dev_functional_args.get(), sizeof(real) * functional_args_size, cudaMemcpyDeviceToHost));

  size_t problem_idx = 0;
  for (size_t thread_idx = 0; thread_idx < count_of_threads; ++thread_idx) {
    size_t records_idx = 0, rec_i = 0; // global and local idxs
    for (size_t case_idx = 0; case_idx < CASES_PER_THREAD; ++case_idx) {

      // prepare new case before writing
      auto t_next = results.Started(problems_vector[problem_idx]);
      Record* rec_block = records[records_idx].get() + thread_idx * ITERS_PER_KERNEL;
      for (unsigned int i = 0; i < STEPS + 1; i++) {
        t_next = results.Store(rec_block[rec_i].t_, rec_block[rec_i].m_, rec_block[rec_i].v_, rec_block[rec_i].h_,
                               rec_block[rec_i].l_, rec_block[rec_i].gamma_);
        rec_i++;
      }

      // write case's necessary records to formatter and functional
      while (true) {

        // go to next records block if current ended
        if (rec_i == ITERS_PER_KERNEL) {
          rec_i = 0;
          records_idx++;
          if (records_idx == records.size()) {
            std::cout << "logically impossible error: uncompleted case\n";
            break;
          }
          rec_block = records[records_idx].get() + thread_idx * ITERS_PER_KERNEL;
        }

        // end of case's records
        if (rec_block[rec_i].t_ == (real)0.0) {
          real* V_args =
              functional_args.get() + (thread_idx * n_timestamps * 2 * CASES_PER_THREAD) + (n_timestamps * 2 * case_idx);
          real* h_args = V_args + n_timestamps;
          size_t timestamp = 0;
          while (timestamp < n_timestamps && V_args[timestamp] != (real)0.0) {
            timestamp++;
          };
          if (rec_block[rec_i].m_ <= (real)0.01) {
            results.Finished(IResultFormatter::Reason::Burnt, functional.Compute(timestamp, V_args, h_args));
          } else if (rec_block[rec_i].h_ <= (real)0.0) {
            results.Finished(IResultFormatter::Reason::Collided, functional.Compute(timestamp, V_args, h_args));
          } else {
            results.Finished(IResultFormatter::Reason::Timeouted, functional.Compute(timestamp, V_args, h_args));
          }
          rec_i = 0; // new case's records starts from next kernel launch
          records_idx++;
          break;
        }

        // store necessary record
        if (rec_block[rec_i].t_ >= t_next) {
          t_next = results.Store(rec_block[rec_i].t_, rec_block[rec_i].m_, rec_block[rec_i].v_, rec_block[rec_i].h_,
                                 rec_block[rec_i].l_, rec_block[rec_i].gamma_);
        }

        rec_i++;
      }
      problem_idx++;
    }
  }
}

// template-specified functions compiles only this way
template void CudaManager<1u, ITERS_PER_KERNEL>(const std::vector<Case>& problems_vector, const IFunctional& functional,
                                                real dt, real timeout, IResultFormatter& results);
template void CudaManager<2u, ITERS_PER_KERNEL>(const std::vector<Case>& problems_vector, const IFunctional& functional,
                                                real dt, real timeout, IResultFormatter& results);
template void CudaManager<3u, ITERS_PER_KERNEL>(const std::vector<Case>& problems_vector, const IFunctional& functional,
                                                real dt, real timeout, IResultFormatter& results);