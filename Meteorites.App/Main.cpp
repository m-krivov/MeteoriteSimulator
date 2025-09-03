#include <chrono>
#include <iostream>

#include "Meteorites.Core/Functionals.h"
#include "Meteorites.Core/ResultFormatters.h"
#include "Meteorites.CpuSolvers/GoldSolver.h"
#include "Meteorites.KnowledgeBase/KnownMeteorites.h"
#include "Meteorites.KnowledgeBase/PossibleParameters.h"
#include "Meteorites.KnowledgeBase/MonteCarloGenerator.h"

//
#include "Meteorites.GpuSolvers/CudaSolver.h"
//

constexpr auto METEORITE         = KnownMeteorites::ID::LOST_CITY;
constexpr auto PARAMETERS        = Distribution::UNIFORM_ANY;
constexpr uint32_t SEED          = 25102018;
constexpr real TIMEOUT           = (real)60.0 * 30;

constexpr size_t STAGE1_N_TOTAL  = 1000;
constexpr size_t STAGE1_N_TOP    = 5;
constexpr auto   STAGE1_METHOD   = NumericalAlgorithm::TWO_STEP_ADAMS;
constexpr real   STAGE1_DT       = (real)1e-3;

constexpr auto   STAGE2_METHOD   = NumericalAlgorithm::THREE_STEP_ADAMS;
constexpr real   STAGE2_DT       = (real)1e-4;


int main()
{
  using clock = std::chrono::high_resolution_clock;

  const IMeteorite* meteorite = KnownMeteorites::Get(METEORITE);
  const ParameterSet& params    = PossibleParameters::Get(PARAMETERS);
  auto progress_callback = [](float) -> void {
    std::cout << '#';
    std::cout.flush();
  };

  // For the first stage, we don't want to simulate meteorite flight till the end
  // So we may use the last record as timeout
  real t_end = (real)0.0;
  {
    size_t records;
    const real *time, *v, *h;
    meteorite->Trajectory(records, time, v, h);
    t_end = time[records - 1];
  }


  //GPU

  // Stage 1.
  std::cout << "Stage 1. Computing huge amount of virtual meteorites with low precision";
  std::cout << std::endl;
  std::cout << "     Cases:  " << STAGE1_N_TOTAL   << " pcs" << std::endl;
  std::cout << "     Method: " << (uint32_t)STAGE1_METHOD << "-step Adams" << std::endl;
  std::cout << "     dt:     " << STAGE1_DT << " seconds" << std::endl;
  std::cout << "Progress: ";
  std::vector<std::pair<Case, double> > good_cases;
  auto started = clock::now();
  {
    MonteCarloGenerator cases(*meteorite, params, STAGE1_N_TOTAL, SEED);
    cases.OnProgress(progress_callback, 0.01f);
    L2Functional functional(*meteorite);
    MetaFormatter meta_fmt(STAGE1_N_TOP, STAGE1_N_TOP * 10);

    CudaSolver cuda_solver;
    cuda_solver.Configure(STAGE1_METHOD, STAGE1_DT, t_end + (real)0.1);
    cuda_solver.Solve(cases, functional, meta_fmt);

    meta_fmt.ExportAndReset(good_cases);
    assert(good_cases.size() == STAGE1_N_TOP);
  }
  auto ended = clock::now();
  std::cout << std::endl;
  std::cout << "Done in " << std::chrono::duration_cast<std::chrono::minutes>(ended - started).count()
            << " minutes" << std::endl << std::endl;

  //CPU

  // Stage 1.
  // Compute trajectories for 'STAGE1_N_TOTAL' virtual meteorites, select 'STAGE1_N_TOTAL' best of them
  /*std::cout << "Stage 1. Computing huge amount of virtual meteorites with low precision";
  std::cout << std::endl;
  std::cout << "     Cases:  " << STAGE1_N_TOTAL   << " pcs" << std::endl;
  std::cout << "     Method: " << (uint32_t)STAGE1_METHOD << "-step Adams" << std::endl;
  std::cout << "     dt:     " << STAGE1_DT << " seconds" << std::endl;
  std::cout << "Progress: ";
  std::vector<std::pair<Case, double> > good_cases;
  auto started = clock::now();
  {
    MonteCarloGenerator cases(*meteorite, params, STAGE1_N_TOTAL, SEED);
    cases.OnProgress(progress_callback, 0.01f);
    L2Functional functional(*meteorite);
    MetaFormatter meta_fmt(STAGE1_N_TOP, STAGE1_N_TOP * 10);

    std::unique_ptr<ISolver> solver(new GoldSolver());
    solver->Configure(STAGE1_METHOD, STAGE1_DT, t_end + (real)0.1);
    solver->Solve(cases, functional, meta_fmt);
    
    meta_fmt.ExportAndReset(good_cases);
    assert(good_cases.size() == STAGE1_N_TOP);
  }
  auto ended = clock::now();
  std::cout << std::endl;
  std::cout << "Done in " << std::chrono::duration_cast<std::chrono::minutes>(ended - started).count()
            << " minutes" << std::endl << std::endl;*/
            
  // Stage 2.
  // Recompute them with better precision, store as tables
  std::cout << "Stage 2. Recomputing trajectories of the best virtual meteorites with high precision";
  std::cout << std::endl;
  std::cout << "     Cases:  " << STAGE1_N_TOP << " pcs" << std::endl;
  std::cout << "     Method: " << (uint32_t)STAGE2_METHOD << "-step Adams" << std::endl;
  std::cout << "     dt:     " << STAGE2_DT << " seconds" << std::endl;
  std::cout << "Progress: ";
  started = clock::now();
  {
    GoldSolver solver;
    solver.Configure(STAGE2_METHOD, STAGE2_DT, TIMEOUT);

    L2Functional functional(*meteorite);
    CsvFromatter csv_fmt(meteorite->Name(), 0.01f);
    size_t on_progress = good_cases.size() / 100;
    for (size_t i = 0; i < good_cases.size(); i++)
    {
      if (i >= on_progress)
      {
        on_progress += good_cases.size() / 100;
        std::cout << "#";
        std::cout.flush();
      }
      solver.Solve(good_cases[i].first, functional, csv_fmt);
    }
  }
  std::cout << "#";
  ended = clock::now();
  std::cout << std::endl;
  std::cout << "Done in " << std::chrono::duration_cast<std::chrono::minutes>(ended - started).count()
            << " minutes" << std::endl << std::endl;

  return 0;
}
