#include "CudaManager.h"
#include "CudaAdamsMethod.h"
#include "CudaStructures.cu"
#include "GPUParameters.h"

static void HandleError(cudaError_t err, const char* file, int line)
{
  if (err != cudaSuccess) {
    printf("%s in %s at line %d\n", cudaGetErrorString(err), file, line);
    exit(EXIT_FAILURE);
  }
}
#define HANDLE_ERROR(err) (HandleError(err, __FILE__, __LINE__))

template <unsigned int STEPS, unsigned int ITERS>
void CudaManager(const std::vector<Case>& problems_vector, const IFunctional& functional, real dt, real timeout,
                 IResultFormatter& results)
{
  // configure grids (make number of cases the main parameter)

  cudaDeviceProp prop;
  HANDLE_ERROR(cudaGetDeviceProperties(&prop, 0));
  unsigned int launching_blocks = prop.multiProcessorCount * BLOCKS_PER_SM;
  unsigned int launching_threads_per_block = THREADS_PER_BLOCK;
  // unsigned int cases_per_thread = problems_vector.size() / launching_blocks / launching_threads_per_block;
  unsigned int cases_per_thread = CASES_PER_THREAD;
  unsigned int count_of_threads = launching_blocks * launching_threads_per_block;
  printf("\nGPU: %s\nConfig:\n>Blocks: %u\n>Threads per block: %u\n>Cases per thread: %u\nTotal threads:%u\nTotal "
         "cases:%u\n",
         prop.name, launching_blocks, launching_threads_per_block, cases_per_thread, count_of_threads,
         problems_vector.size());
  printf("\nDid you update parameters in GPUParameters.h? If yes tap any key to continue...\n");
  getchar();

  // preparing buffers

  // thread contexts
  ThreadContext<STEPS, ITERS>* dev_thread_sandbox_arr;
  size_t thread_sandbox_arr_size = sizeof(ThreadContext<STEPS, ITERS>) * count_of_threads;
  HANDLE_ERROR(cudaMalloc(&dev_thread_sandbox_arr, thread_sandbox_arr_size));

  // problems (read-only)
  CudaCase* dev_problems;
  size_t problems_size = sizeof(CudaCase) * problems_vector.size();
  HANDLE_ERROR(cudaMalloc(&dev_problems, problems_size));
  HANDLE_ERROR(
      cudaMemcpy(dev_problems, problems_vector.data(), problems_size, cudaMemcpyHostToDevice)); // damn reinterpret cast

  // meteorite's timestamps (read-only)
  const real* timestamps = nullptr;
  size_t n_timestamps = 0;
  functional.GetTimeStamps(n_timestamps, timestamps);
  real* dev_timestamps;
  HANDLE_ERROR(cudaMalloc(&dev_timestamps, sizeof(real) * n_timestamps));
  HANDLE_ERROR(cudaMemcpy(dev_timestamps, timestamps, sizeof(real) * n_timestamps, cudaMemcpyHostToDevice));

  // (V_arg[] and h_arg[] for functional)'s for each case
  size_t functional_args_size = n_timestamps * sizeof(real) * 2 * count_of_threads * CASES_PER_THREAD;
  real* functional_args = (real*)malloc(functional_args_size);
  real* dev_functional_args;
  HANDLE_ERROR(cudaMalloc(&dev_functional_args, functional_args_size));

  // vector of Record's arrays
  // load to CPU vector every iteration
  std::vector<Record*> records;
  size_t records_size = sizeof(Record) * ITERS_PER_KERNEL * count_of_threads;
  Record* dev_records;
  HANDLE_ERROR(cudaMalloc(&dev_records, records_size));

  // active threads atomic counter (make reduction (or no))
  unsigned int active_threads = count_of_threads, *dev_active_threads;
  HANDLE_ERROR(cudaMalloc(&dev_active_threads, sizeof(unsigned int)));
  HANDLE_ERROR(cudaMemcpy(dev_active_threads, &active_threads, sizeof(unsigned int), cudaMemcpyHostToDevice));

  int global_iter = 0;
  while (active_threads) {
    global_iter++;

    CudaLauncher<STEPS, ITERS>(dev_thread_sandbox_arr, dev_problems, problems_vector.size(), dt, timeout,
                               dev_timestamps, n_timestamps, dev_functional_args, dev_records, dev_active_threads,
                               launching_blocks, launching_threads_per_block);

    cudaDeviceSynchronize();

    records.push_back((Record*)malloc(records_size));
    HANDLE_ERROR(cudaMemcpy(records[records.size() - 1], dev_records, records_size, cudaMemcpyDeviceToHost));

    HANDLE_ERROR(cudaMemcpy(&active_threads, dev_active_threads, sizeof(active_threads), cudaMemcpyDeviceToHost));
  }

  cudaFree(dev_thread_sandbox_arr);
  cudaFree(dev_problems);
  cudaFree(dev_timestamps);
  cudaFree(dev_records);
  cudaFree(dev_active_threads);

  // process records

  HANDLE_ERROR(cudaMemcpy(functional_args, dev_functional_args, functional_args_size, cudaMemcpyDeviceToHost));

  size_t problem_idx = 0;
  for (size_t thread_idx = 0; thread_idx < count_of_threads; ++thread_idx) {
    size_t records_idx = 0, rec_i = 0; // global and local idxs
    for (size_t case_idx = 0; case_idx < CASES_PER_THREAD; ++case_idx) {

      // prepare new case before writing
      auto t_next = results.Started(problems_vector[problem_idx]);
      Record* rec_block = records[records_idx] + thread_idx * ITERS_PER_KERNEL;
      for (unsigned int i = 0; i < STEPS + 1; i++) {
        /*test*/ if (problem_idx == 0) {
          rec_block[rec_i].print();
        }
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
          rec_block = records[records_idx] + thread_idx * ITERS_PER_KERNEL;
        }

        // end of case's records
        if (rec_block[rec_i].t_ == (real)0.0) {
          real* V_args =
              functional_args + (thread_idx * n_timestamps * 2 * CASES_PER_THREAD) + (n_timestamps * 2 * case_idx);
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

  cudaFree(dev_functional_args);

  free(functional_args);
  for (auto i = records.begin(); i != records.end(); ++i) {
    free(*i);
  }
}

// template-specified functions compiles only this way
template void CudaManager<1u, ITERS_PER_KERNEL>(const std::vector<Case>& problems_vector, const IFunctional& functional,
                                                real dt, real timeout, IResultFormatter& results);
template void CudaManager<2u, ITERS_PER_KERNEL>(const std::vector<Case>& problems_vector, const IFunctional& functional,
                                                real dt, real timeout, IResultFormatter& results);
template void CudaManager<3u, ITERS_PER_KERNEL>(const std::vector<Case>& problems_vector, const IFunctional& functional,
                                                real dt, real timeout, IResultFormatter& results);