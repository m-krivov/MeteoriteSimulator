#include <chrono>
#include <iostream>

#include "Meteorites.Core/Functionals.h"
#include "Meteorites.Core/ResultFormatters.h"
#include "Meteorites.Core/ISolver.h"
#include "Meteorites.KnowledgeBase/KnownMeteorites.h"
#include "Meteorites.KnowledgeBase/PossibleParameters.h"

#include "Meteorites.CpuSolvers/GoldSolver.h"
#include "Meteorites.KnowledgeBase/MonteCarloGenerator.h"

#include "Meteorites.GpuSolvers/CudaManagers.h"
#include "Meteorites.GpuSolvers/AdamsMethod.h" // BEST_CASES_BUFFER_SIZE defined here
#include "Meteorites.GpuSolvers/CudaStructures.h"

#include <curand_kernel.h>

constexpr auto METEORITE = KnownMeteorites::ID::LOST_CITY;

constexpr uint32_t STAGE1_CASES_PER_THREAD  = 1024 * 8;
constexpr uint32_t STAGE1_THREADS_PER_BLOCK = 32 * 2;
constexpr uint32_t STAGE1_BLOCKS_NUM        = 128 * 2;
constexpr uint64_t STAGE1_N_TOTAL = STAGE1_BLOCKS_NUM * STAGE1_THREADS_PER_BLOCK * STAGE1_CASES_PER_THREAD;
constexpr uint32_t STAGE1_N_TOP   = STAGE1_BLOCKS_NUM * STAGE1_THREADS_PER_BLOCK * BEST_CASES_BUFFER_SIZE;
constexpr auto     STAGE1_METHOD  = NumericalAlgorithm::TWO_STEP_ADAMS;
constexpr real     STAGE1_DT      = (real)1e-3;

constexpr auto     STAGE2_METHOD  = NumericalAlgorithm::THREE_STEP_ADAMS;
constexpr real     STAGE2_DT      = (real)1e-4;
constexpr uint32_t STAGE2_N_TOTAL = 100;

int main()
{
  using clock = std::chrono::high_resolution_clock;

  const IMeteorite *meteorite  = KnownMeteorites::Get(METEORITE);

  real STAGE1_TIMEOUT = (real)0.0;

  size_t records;
  const real *time, *v, *h;
  meteorite->Trajectory(records, time, v, h);
  STAGE1_TIMEOUT = time[records - 1] + STAGE1_DT + 30 * 0;

  //GPU

  // Stage 1.
  std::cout << "Stage 1. Computing huge amount of virtual meteorites with low precision. Fixed small number of AdamsSteps";
  std::cout << std::endl;
  std::cout << "\tCases:   " << STAGE1_N_TOTAL   << " pcs" << std::endl;
  std::cout << "\tMethod:  " << (uint32_t)STAGE1_METHOD << "-step Adams" << std::endl;
  std::cout << "\tTimeout: " << STAGE1_TIMEOUT << " seconds" << std::endl;
  std::cout << "\tdt:      " << STAGE1_DT << " seconds" << std::endl;

  printf("\nGPU config:\n\tBlocks:%u\n\tThreads per block:%u\n\tCases per thread:%u\n",
         STAGE1_BLOCKS_NUM, STAGE1_THREADS_PER_BLOCK, STAGE1_CASES_PER_THREAD);

  auto started = clock::now();

  auto seeds = std::make_unique<uint64_t[]>(STAGE1_BLOCKS_NUM * STAGE1_THREADS_PER_BLOCK);
  auto result_curand_states = std::make_unique<curandState[]>(STAGE1_N_TOP);
  auto result_deviations = std::make_unique<real[]>(STAGE1_N_TOP);
  CudaAdamsMethodManager<(uint32_t)STAGE1_METHOD>
      (*meteorite, seeds.get(), STAGE1_BLOCKS_NUM * STAGE1_THREADS_PER_BLOCK,
       STAGE1_DT, STAGE1_TIMEOUT, 0.1,
       STAGE1_BLOCKS_NUM, STAGE1_THREADS_PER_BLOCK,
       BEST_CASES_BUFFER_SIZE, STAGE1_CASES_PER_THREAD,
       result_curand_states.get(), result_deviations.get());

  std::vector<std::pair<real, curandState>> pairs;
  pairs.reserve(STAGE1_N_TOP);

  printf("Sort best cases...\n");

  for (uint32_t i = 0; i < STAGE1_N_TOP; i++)
  {
      pairs.emplace_back(result_deviations[i], result_curand_states[i]);
  }
  std::sort(pairs.begin(), pairs.end(), [](const auto& a, const auto& b) { return a.first < b.first; });

#define PRINT_TOP
#ifdef PRINT_TOP
  int display_top_conter = 0;
  printf("\nSTAGE1 results:\n");
  for (uint32_t i = 0; i < STAGE1_N_TOP; i++) {
    result_deviations[i] = pairs[i].first;
    result_curand_states[i] = pairs[i].second;
    if (display_top_conter < 10)
    {
      printf("\ttop%d: dev:%f\n", i, result_deviations[i]);
      display_top_conter++;
    }
  }
  printf("\t...\n\ttop%d: dev:%f\n", STAGE1_N_TOP - 2,  result_deviations[STAGE1_N_TOP - 2]);
  printf("\ttop%d: dev:%f\n", STAGE1_N_TOP - 1,  result_deviations[STAGE1_N_TOP - 1]);
#endif

  auto ended = clock::now();
  std::cout << std::endl;
  std::cout << "Done in " << std::chrono::duration_cast<std::chrono::minutes>(ended - started).count()
            << " minutes" << std::endl << std::endl;


  // CPU

  // Stage 2.
  std::cout << "Stage 2. Recomputing trajectories of the best virtual meteorites with high precision";
  std::cout << std::endl;
  std::cout << "\tCases:  " << STAGE2_N_TOTAL << " pcs" << std::endl;
  std::cout << "\tMethod: " << (uint32_t)STAGE2_METHOD << "-step Adams" << std::endl;
  std::cout << "\tdt:     " << STAGE2_DT << " seconds" << std::endl;

  started = clock::now();
  {
    std::vector<Case> good_cases = CudaRestoreCasesManager(result_curand_states.get(), STAGE2_N_TOTAL, v[0], h[0]);

    GoldSolver solver;
    solver.Configure(STAGE2_METHOD, STAGE2_DT, STAGE1_TIMEOUT);

    L2Functional functional(*meteorite);
    CsvFromatter csv_fmt(meteorite->Name(), 0.01f);
    for (size_t i = 0; i < good_cases.size(); i++)
    {
      solver.Solve(good_cases[i], functional, csv_fmt);
    }
  }
  ended = clock::now();
  std::cout << std::endl;
  std::cout << "Done in " << std::chrono::duration_cast<std::chrono::minutes>(ended - started).count()
            << " minutes" << std::endl << std::endl;
}
