#include <chrono>
#include <iostream>
#include <indicators/cursor_control.hpp>
#include <indicators/progress_bar.hpp>

#include "Meteorites.Core/Functionals.h"
#include "Meteorites.Core/ResultFormatters.h"
#include "Meteorites.Core/Meteoroids/CollectionMeteoroidGenerator.h"
#include "Meteorites.CpuSolvers/GoldSolver.h"
#include "Meteorites.KnowledgeBase/KnownMeteorites.h"
#include "Meteorites.KnowledgeBase/PossibleParameters.h"
#include "Meteorites.KnowledgeBase/MonteCarloGenerator.h"

#if defined(METEORITES_CUDA)
  #include "Meteorites.CudaSolvers/PedanticCudaSolver.h"
#endif


constexpr auto   METEORITE       = KnownMeteorites::ID::INNISFREE;
constexpr auto   PARAMETERS      = Distribution::UNIFORM_ANY;
constexpr uint32_t SEED          = 25102018;
constexpr real   TIMEOUT         = (real)60.0 * 30;

constexpr size_t STAGE1_N_TOTAL  = 1000000;
constexpr size_t STAGE1_N_TOP    = 100;
constexpr auto   STAGE1_METHOD   = NumericalAlgorithm::TWO_STEP_ADAMS;
constexpr real   STAGE1_DT       = (real)1e-3;

constexpr auto   STAGE2_METHOD   = NumericalAlgorithm::THREE_STEP_ADAMS;
constexpr real   STAGE2_DT       = (real)1e-4;

#if defined(METEORITES_CUDA)
  constexpr bool   USE_GPU = true;
#endif


// Configure progress bar from a 3rd-party library
// Just an old-school bars without colors and animations
class MyProgressBar : public indicators::ProgressBar
{
  public:
    MyProgressBar()
      : ProgressBar(indicators::option::BarWidth{ 60 },
                    indicators::option::Fill{ "#" },
                    indicators::option::Lead{ "#" },
                    indicators::option::Remainder{ "-" },
                    indicators::option::PrefixText{ "     " },
                    indicators::option::ShowPercentage{ true },
                    indicators::option::ShowElapsedTime{ true },
                    indicators::option::ShowRemainingTime{ true })
    { }
};


int main()
{
  using clock = std::chrono::high_resolution_clock;

  decltype(auto) meteorite = KnownMeteorites::Ref().Get(METEORITE);
  decltype(auto) params    = PossibleParameters::Get(PARAMETERS);
  auto progress_callback = [](float) -> void {
    std::cout << '#';
    std::cout.flush();
  };

  // Print all available information about the simulated meteorite
  std::cout << "Meteorite: " << meteorite.Name() << std::endl;
  std::cout << "     Date: " << meteorite.Date() << std::endl;
  std::cout << "     Fall: " << meteorite.FallLocation() << std::endl;
  std::cout << "     DOI:  " << meteorite.DOI() << std::endl;
  std::cout << std::endl;

  // For the first stage, we don't want to simulate meteorite flight till the end
  // So we may use the last record as timeout
  real t_end = (real)0.0;
  {
    size_t records;
    const real *time, *v, *h;
    meteorite.Trajectory(records, time, v, h);
    t_end = time[records - 1];
  }

  indicators::show_console_cursor(false);
  try
  {
    // Stage 1.
    // Compute trajectories for 'STAGE1_N_TOTAL' virtual meteoroids, select 'STAGE1_N_TOTAL' best of them
    std::cout << "Stage 1. Computing huge amount of virtual meteoroids with low precision";
    std::cout << std::endl;
    std::cout << "     Meteoroids: " << STAGE1_N_TOTAL   << " pcs" << std::endl;
    std::cout << "     Method:     " << (uint32_t)STAGE1_METHOD << "-step Adams" << std::endl;
    std::cout << "     dt:         " << STAGE1_DT << " seconds" << std::endl;
    std::vector<std::pair<VirtualMeteoroid, double> > good_meteoroids;
    {
      MonteCarloGenerator generator(meteorite, params, STAGE1_N_TOTAL, SEED);
      generator.OnProgress
      (
        [bar = std::make_shared<MyProgressBar>()](float progress) mutable -> void
        { bar->set_progress((size_t)(100 * progress)); }, 0.01f
      );
      L2Functional functional(meteorite);
      MetaFormatter meta_fmt(STAGE1_N_TOP, STAGE1_N_TOP * 10);

      std::unique_ptr<ISolver> solver;
    #if defined(METEORITES_CUDA)
      if constexpr (USE_GPU)
      { solver.reset(new PedanticCudaSolver()); }
      else
    #endif
      { solver.reset(new GoldSolver()); }

      solver->Configure(STAGE1_METHOD, STAGE1_DT, t_end + (real)0.1);
      solver->Solve(generator, functional, meta_fmt);
    
      meta_fmt.ExportAndReset(good_meteoroids);
      assert(good_meteoroids.size() == STAGE1_N_TOP);
    }
    std::cout << std::endl;

    // Stage 2.
    // Recompute them with better precision, store as tables
    std::cout << "Stage 2. Recomputing trajectories of the best virtual meteoroids with high precision";
    std::cout << std::endl;
    std::cout << "     Meteoroids: " << STAGE1_N_TOP << " pcs" << std::endl;
    std::cout << "     Method:     " << (uint32_t)STAGE2_METHOD << "-step Adams" << std::endl;
    std::cout << "     dt:         " << STAGE2_DT << " seconds" << std::endl;
    {
      GoldSolver solver;
      solver.Configure(STAGE2_METHOD, STAGE2_DT, TIMEOUT);

      L2Functional functional(meteorite);
      CsvFromatter csv_fmt(meteorite.Name(), 0.01f);
      CollectionMeteoroidGenerator generator(good_meteoroids);
      generator.OnProgress
      (
        [bar = std::make_shared<MyProgressBar>()](float progress) mutable -> void
        { bar->set_progress((size_t)(100 * progress)); }, 0.01f
      );
      ((ISolver &)solver).Solve(generator, functional, csv_fmt);
    }
    std::cout << std::endl;
    std::cout << "Success!" << std::endl;
  }
  catch (std::exception &ex)
  { std::cerr << "Error: " << ex.what(); }
  indicators::show_console_cursor(true);
  
  return 0;
}
