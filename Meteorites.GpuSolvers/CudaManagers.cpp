#include "CudaManagers.h"

#include "Meteorites.Core/IMeteorite.h"
#include "Meteorites.Core/Case.h"

#include "CudaHandleError.h"
#include "CudaUniquePtr.h"
#include "AdamsMethod.h"
#include "CudaStructures.h"
#include "RestoreCases.h"
#include "Random.h"

void GenerateSeeds(uint64_t* seeds, int N)
{
  Random::State rootState;
  Random::Initialize(rootState);
        
  for (int i = 0; i < N; i++)
  {
    Random::State tempState;
    std::vector<Random::State> tempVector(1);
    Random::Multiply(rootState, tempVector);
    tempState = tempVector[0];
            
    uint32_t part1 = Random::Next(tempState);
    uint32_t part2 = Random::Next(tempState);
    seeds[i] = (static_cast<uint64_t>(part1) << 32) | part2;
            
    Random::Next(rootState);
  }
}

template <uint32_t STEPS>
void CudaAdamsMethodManager(const IMeteorite &meteorite, uint64_t *seeds, int N,
                            real dt, real timeout, real border_dev,
                            uint32_t BLOCKS_NUM, uint32_t THREADS_PER_BLOCK,
                            uint32_t BEST_CASES_BUFFER_SIZE, uint32_t CASES_PER_THREAD,
                            curandState *result_curand_states, real *result_deviations)
{
  size_t n_timestamps = 0;
  const real *timestamps, *v, *h;
  meteorite.Trajectory(n_timestamps, timestamps, v, h);

  // prepare seeds for GPU
  /*test*/printf("Generate seeds on CPU...\n");
  GenerateSeeds(seeds, N);

  auto dev_seeds = cuda_make_unique<uint64_t>(N);
  HANDLE_ERROR(cudaMemcpy(dev_seeds.get(), seeds, sizeof(uint64_t) * N, cudaMemcpyHostToDevice));

  // prepare context for GPU
  auto dev_timestamps = cuda_make_unique<real>(n_timestamps);
  HANDLE_ERROR(cudaMemcpy(dev_timestamps.get(), timestamps, sizeof(real) * n_timestamps, cudaMemcpyHostToDevice));
  auto dev_ref_v = cuda_make_unique<real>(n_timestamps);
  HANDLE_ERROR(cudaMemcpy(dev_ref_v.get(), v, sizeof(real) * n_timestamps, cudaMemcpyHostToDevice));
  auto dev_ref_h = cuda_make_unique<real>(n_timestamps);
  HANDLE_ERROR(cudaMemcpy(dev_ref_h.get(), h, sizeof(real) * n_timestamps, cudaMemcpyHostToDevice));
  uint32_t GLOBAL_BEST_BUFFER_SIZE = BLOCKS_NUM * THREADS_PER_BLOCK * BEST_CASES_BUFFER_SIZE;
  auto dev_best_cases_devs = cuda_make_unique<real>(GLOBAL_BEST_BUFFER_SIZE);
  auto dev_best_cases_states = cuda_make_unique<curandState>(GLOBAL_BEST_BUFFER_SIZE);
  auto dev_v_args_arr = cuda_make_unique<real>(BLOCKS_NUM * THREADS_PER_BLOCK * n_timestamps);
  auto dev_h_args_arr = cuda_make_unique<real>(BLOCKS_NUM * THREADS_PER_BLOCK * n_timestamps);

  cudaEvent_t start, stop;
  float elapsedTime;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);

  /*test*/printf("Start calculations...\n");
  CudaAdamsMethodLauncher<STEPS>(dev_seeds.get(), dev_timestamps.get(), n_timestamps,
                                 dev_ref_v.get(), dev_ref_h.get(),
                                 dt, timeout, dev_best_cases_devs.get(), dev_best_cases_states.get(),
                                 border_dev, dev_v_args_arr.get(), dev_h_args_arr.get(),
                                 CASES_PER_THREAD, BLOCKS_NUM, THREADS_PER_BLOCK);
  
  cudaEventRecord(stop,0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&elapsedTime, start, stop);
  /*test*/printf("Time spent executing by the GPU: %.2f milliseconds\n", elapsedTime);
  cudaEventDestroy(start);
  cudaEventDestroy(stop);

  // load best deviations and indices from GPU, find seeds, return
  /*test*/printf("Load top results to CPU...\n");
  HANDLE_ERROR(cudaMemcpy(result_curand_states, dev_best_cases_states.get(),
                         sizeof(curandState) * GLOBAL_BEST_BUFFER_SIZE, cudaMemcpyDeviceToHost));
  HANDLE_ERROR(cudaMemcpy(result_deviations, dev_best_cases_devs.get(),
                         sizeof(real) * GLOBAL_BEST_BUFFER_SIZE, cudaMemcpyDeviceToHost));
}

template void CudaAdamsMethodManager<1u>(const IMeteorite &meteorite, uint64_t *seeds, int N,
                                         real dt, real timeout, real border_dev,
                                         uint32_t BLOCKS_NUM, uint32_t THREADS_PER_BLOCK,
                                         uint32_t BEST_CASES_BUFFER_SIZE, uint32_t CASES_PER_THREAD,
                                         curandState *result_curand_states, real *result_deviations);
template void CudaAdamsMethodManager<2u>(const IMeteorite &meteorite, uint64_t *seeds, int N,
                                         real dt, real timeout, real border_dev,
                                         uint32_t BLOCKS_NUM, uint32_t THREADS_PER_BLOCK,
                                         uint32_t BEST_CASES_BUFFER_SIZE, uint32_t CASES_PER_THREAD,
                                         curandState *result_curand_states, real *result_deviations);
template void CudaAdamsMethodManager<3u>(const IMeteorite &meteorite, uint64_t *seeds, int N,
                                         real dt, real timeout, real border_dev,
                                         uint32_t BLOCKS_NUM, uint32_t THREADS_PER_BLOCK,
                                         uint32_t BEST_CASES_BUFFER_SIZE, uint32_t CASES_PER_THREAD,
                                         curandState *result_curand_states, real *result_deviations);

std::vector<Case> CudaRestoreCasesManager(const curandState *states, int N, real v0, real h0)
{
  auto dev_cases = cuda_make_unique<Case>(N);
  auto dev_states = cuda_make_unique<curandState>(N);
  HANDLE_ERROR(cudaMemcpy(dev_states.get(), states, sizeof(curandState) * N, cudaMemcpyHostToDevice));
  CudaRestoreCasesLauncher(dev_cases.get(), dev_states.get(), N, v0, h0);
  auto cases = std::make_unique<Case[]>(N);
  HANDLE_ERROR(cudaMemcpy(cases.get(), dev_cases.get(), sizeof(Case) * N, cudaMemcpyDeviceToHost));
  std::vector<Case> good_cases;
  good_cases.reserve(N);
  std::move(cases.get(), cases.get() + N, std::back_inserter(good_cases));
  return good_cases;
}