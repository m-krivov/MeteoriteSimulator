#include "Meteorites.Core/Defs.h"

#include <iostream>
#include <indicators/cursor_control.hpp>
#include <indicators/progress_bar.hpp>
#include <indicators/progress_spinner.hpp>

#include "Meteorites.Core/Functionals/CFunctional.h"
#include "Meteorites.Core/Functionals/L1Functional.h"
#include "Meteorites.Core/Functionals/L2Functional.h"
#include "Meteorites.Core/Functionals/BurnTimePenaltyFunctional.h"
#include "Meteorites.Core/Recorders/BufferingRecorder.h"
#include "Meteorites.Core/Recorders/MetaRecorder.h"
#include "Meteorites.Core/Meteoroids/CollectionGenerator.h"
#include "Meteorites.Core/Meteoroids/MonteCarloGenerator.h"
#include "Meteorites.CpuSolvers/GoldSolver.h"
#include "Meteorites.KnowledgeBase/KnownMeteorites.h"
#include "Meteorites.KnowledgeBase/PossibleParameters.h"

#if defined(METEORITES_CUDA)
  #include "Meteorites.CudaSolvers/PedanticCudaSolver.h"
#endif

#include "Exporters/MeanStdevExporter.h"
#include "Exporters/TrajectoryExporter.h"
#include "Exporters/TrajectoryVisualizer.h"
#include "Exporters/HistogramBuilder.h"


constexpr auto     METEORITE         = KnownMeteorites::ID::INNISFREE;
constexpr auto     PARAMETERS        = Distribution::UNIFORM_ANY;
constexpr uint32_t SEED              = 25102018;
constexpr float    PROGRESS_BAR_STEP = 0.01f;
constexpr real     TIMEOUT           = (real)60.0 * 30;

// Stage 0: Quick scan to optimize entry angle range
constexpr size_t   STAGE0_N_TOTAL    = 10000;
constexpr size_t   STAGE0_N_TOP      = 100;
constexpr auto     STAGE0_METHOD     = NumericalAlgorithm::ONE_STEP_ADAMS;
constexpr real     STAGE0_DT         = (real)1e-2;
using              STAGE0_FUNC       = L2Functional;

constexpr size_t   STAGE1_N_TOTAL    = 1000000;
constexpr size_t   STAGE1_N_TOP      = 1000;
constexpr auto     STAGE1_METHOD     = NumericalAlgorithm::TWO_STEP_ADAMS;
constexpr real     STAGE1_DT         = (real)1e-3;
using              STAGE1_FUNC       = L2Functional;
// Note: To use weighted functionals (for meteorites where later measurements are less accurate):
// 1. Get measurement count from meteorite trajectory
// 2. Generate decaying weights: auto weights = BasicFunctional::GenerateDecayingWeights(n_records, 0.5);
// 3. Create functional: L2Functional functional(meteorite, 1.0, 1.0, weights);
// This gives exponentially decaying weights to later measurements (as recommended for Halliday1996 data)

constexpr size_t   STAGE2_N_TOP      = 100;
constexpr auto     STAGE2_METHOD     = NumericalAlgorithm::THREE_STEP_ADAMS;
constexpr real     STAGE2_DT         = (real)1e-4;
using              STAGE2_FUNC       = L2Functional;

constexpr size_t   EXPORT_N_GROUPS   = 4;
constexpr real     EXPORT_DT         = (real)0.1;

#if defined(METEORITES_CUDA)
  constexpr bool   USE_GPU = true;
#endif


// Configure progress bar from a 3rd-party library
// Just an old-school bars without colors and animations
class MyProgressBar : public indicators::ProgressBar
{
  public:
    MyProgressBar()
      : ProgressBar
      (
        indicators::option::BarWidth{ 60 },
        indicators::option::Fill{ "#" },
        indicators::option::Lead{ "#" },
        indicators::option::Remainder{ "-" },
        indicators::option::PrefixText{ "     " },
        indicators::option::ShowPercentage{ true },
        indicators::option::ShowElapsedTime{ true },
        indicators::option::ShowRemainingTime{ true }
      )
    { }

    // Saves a few lines of code
    static std::function<void(float)> Create()
    {
      return [bar = std::make_shared<MyProgressBar>()](float progress) mutable -> void
      { bar->set_progress((size_t)(100 * progress)); };
    }
};

// The same for spinner: simple animation, percentage and nothing more
class MyProgressSpinner : public indicators::ProgressSpinner
{
  public:
    MyProgressSpinner(const std::string &title)
      : ProgressSpinner
      (
        indicators::option::PrefixText{title},
        indicators::option::SpinnerStates{std::vector<std::string>{" ", " " }},
        indicators::option::ShowSpinner(true),
        indicators::option::ShowPercentage(true)
      )
    {}

    static std::function<void(float)> Create(const std::string &title)
    {
      return [spinner = std::make_shared<MyProgressSpinner>(title)](float progress) mutable -> void
      { spinner->set_progress((size_t)(100 * progress)); };
    }

};

std::string ToString(NumericalAlgorithm alg)
{
  switch (alg)
  {
    case NumericalAlgorithm::ONE_STEP_ADAMS:
      return "One-step Adams";

    case NumericalAlgorithm::TWO_STEP_ADAMS:
      return "Two-step Adams";

    case NumericalAlgorithm::THREE_STEP_ADAMS:
      return "Three-step Adams";

    default:
      return "Unknown algorithm";
  }
}

std::unique_ptr<ISolver> CreateSolver()
{
#if defined(METEORITES_CUDA)
  if constexpr (USE_GPU)
  { return std::make_unique<PedanticCudaSolver>(); }
  else
#endif
  { return std::make_unique<GoldSolver>(); }
}


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
  std::cout << "Meteorite: " << meteorite->Name() << std::endl;
  std::cout << "     Date: " << meteorite->Date() << std::endl;
  std::cout << "     Fall: " << meteorite->FallLocation() << std::endl;
  std::cout << "     DOI:  " << meteorite->DOI() << std::endl;
  std::cout << std::endl;

  // For the first stage, we don't want to simulate meteorite flight till the end
  // So we may use the last record as timeout
  real t_end = (real)0.0;
  {
    size_t records;
    const real *time, *v, *h;
    meteorite->Trajectory(records, time, v, h);
    t_end = time[records - 1];
  }

  indicators::show_console_cursor(false);
  try
  {
    // Stage 0.
    // Quick scan to find optimal entry angle range
    std::cout << "Stage 0. Quick scan to optimize entry angle range";
    std::cout << std::endl;
    std::cout << "     Meteoroids: " << STAGE0_N_TOTAL   << " pcs" << std::endl;
    std::cout << "     Method:     " << ToString(STAGE0_METHOD) << std::endl;
    std::cout << "     dt:         " << STAGE0_DT << " seconds" << std::endl;
    std::vector<std::pair<VirtualMeteoroid, double> > stage0_meteoroids;
    real gamma_min = (real)90.0, gamma_max = (real)0.0;
    {
      MonteCarloGenerator generator(meteorite, params, STAGE0_N_TOTAL, SEED);
      generator.OnProgress(MyProgressBar::Create(), PROGRESS_BAR_STEP);
      STAGE0_FUNC functional(meteorite);
      MetaRecorder recorder(STAGE0_N_TOP, STAGE0_N_TOP * 10);

      std::unique_ptr<ISolver> solver = CreateSolver();
      solver->Configure(STAGE0_METHOD, STAGE0_DT, t_end + (real)0.1);
      solver->Solve(generator, functional, recorder);
    
      recorder.ExportAndReset(stage0_meteoroids);
      assert(stage0_meteoroids.size() == STAGE0_N_TOP);
      
      // Analyze the best meteoroids to find optimal gamma range
      for (const auto &meteoroid : stage0_meteoroids)
      {
        real gamma_deg = meteoroid.first.Gamma0 * (real)180.0 / (real)M_PI;
        gamma_min = std::min(gamma_min, gamma_deg);
        gamma_max = std::max(gamma_max, gamma_deg);
      }
      
      // Add some margin (10%) to the range
      real gamma_margin = (gamma_max - gamma_min) * (real)0.1;
      gamma_min = std::max((real)0.0, gamma_min - gamma_margin);
      gamma_max = std::min((real)90.0, gamma_max + gamma_margin);
      
      std::cout << "     Optimal entry angle range: [" << gamma_min << ", " << gamma_max << "] degrees" << std::endl;
    }
    std::cout << std::endl;

    // Update parameters with refined gamma range
    real gamma_min_rad = gamma_min * (real)M_PI / (real)180.0;
    real gamma_max_rad = gamma_max * (real)M_PI / (real)180.0;
    ParameterSet refined_params(params.H(), params.Ch(), params.Rho(), 
                                params.Cd(), params.Cl(), params.M0(),
                                std::make_pair(gamma_min_rad, gamma_max_rad));

    // Stage 1.
    // Compute trajectories for 'STAGE1_N_TOTAL' virtual meteoroids, select 'STAGE1_N_TOTAL' best of them
    std::cout << "Stage 1. Computing huge amount of virtual meteoroids with low precision";
    std::cout << std::endl;
    std::cout << "     Meteoroids: " << STAGE1_N_TOTAL   << " pcs" << std::endl;
    std::cout << "     Method:     " << ToString(STAGE1_METHOD) << std::endl;
    std::cout << "     dt:         " << STAGE1_DT << " seconds" << std::endl;
    std::vector<std::pair<VirtualMeteoroid, double> > stage1_meteoroids;
    {
      MonteCarloGenerator generator(meteorite, refined_params, STAGE1_N_TOTAL, SEED);
      generator.OnProgress(MyProgressBar::Create(), PROGRESS_BAR_STEP);
      STAGE1_FUNC functional(meteorite);
      MetaRecorder recorder(STAGE1_N_TOP, STAGE1_N_TOP * 10);

      std::unique_ptr<ISolver> solver = CreateSolver();
      solver->Configure(STAGE1_METHOD, STAGE1_DT, t_end + (real)0.1);
      solver->Solve(generator, functional, recorder);
    
      recorder.ExportAndReset(stage1_meteoroids);
      assert(stage1_meteoroids.size() == STAGE1_N_TOP);
    }
    std::cout << std::endl;

    // Stage 2.
    // Recompute them with better precision, store their trajectories
    // After that, sort them again and truncate
    std::cout << "Stage 2. Recomputing trajectories of the best virtual meteoroids with high precision";
    std::cout << std::endl;
    std::cout << "     Meteoroids: " << STAGE1_N_TOP << " pcs" << std::endl;
    std::cout << "     Method:     " << ToString(STAGE2_METHOD) << std::endl;
    std::cout << "     dt:         " << STAGE2_DT << " seconds" << std::endl;
    std::vector<MeteoroidTrajectory> stage2_trajectories;
    {
      std::unique_ptr<ISolver> solver = std::make_unique<GoldSolver>();
      solver->Configure(STAGE2_METHOD, STAGE2_DT, TIMEOUT);

      STAGE2_FUNC functional(meteorite);
      BufferingRecorder recorder(EXPORT_DT);
      CollectionGenerator generator(stage1_meteoroids);
      generator.OnProgress(MyProgressBar::Create(), PROGRESS_BAR_STEP);
      stage1_meteoroids.clear();
      solver->Solve(generator, functional, recorder);
      
      recorder.MoveTo(stage2_trajectories);
      assert(stage2_trajectories.size() == STAGE1_N_TOP);
      std::sort(stage2_trajectories.begin(), stage2_trajectories.end(),
                [](const MeteoroidTrajectory &a, const MeteoroidTrajectory &b) -> bool
                { return a.Accuracy() < b.Accuracy(); });
      
      assert(STAGE2_N_TOP <= stage2_trajectories.size());
      stage2_trajectories.erase(stage2_trajectories.begin() + STAGE2_N_TOP,
                                stage2_trajectories.end());
    }
    std::cout << std::endl;

    // Done! Save the results using different exporters
    std::cout << "Saving the results" << std::endl;
    std::cout << "     Meteoroids: " << STAGE2_N_TOP << " pcs" << std::endl;
    std::vector<std::pair<std::string, std::shared_ptr<IExporter>>> exporters
    {
      std::make_pair(std::string("Represent meteoroids' parameters as mean and standard deviation"),
                     std::shared_ptr<IExporter>(new MeanStdevExporter(EXPORT_N_GROUPS))),
      std::make_pair(std::string("Store trajectories as tables"),
                     std::shared_ptr<IExporter>(new TrajectoryExporter(EXPORT_DT))),
    };
#if defined(METEORITES_GNUPLOT)
    exporters.emplace_back
    (
      std::make_pair(std::string("Visualize trajectories as images"),
                     std::shared_ptr<IExporter>(new TrajectoryVisualizer(meteorite)))
    );
    exporters.emplace_back
    (
      std::make_pair(std::string("Save meteoroid parameter distributions as histograms"),
                     std::shared_ptr<IExporter>(new HistogramBuilder(meteorite)))
    );
#endif

    for (const auto &rec : exporters)
    {
      rec.second->OnProgress(MyProgressSpinner::Create(std::string("     ") + rec.first), PROGRESS_BAR_STEP);
      rec.second->SetDirectory(std::filesystem::current_path() / meteorite->Name());
      rec.second->SetMetaData("now", meteorite); // TODO: use C++20 and 'date' to format the actual time
      rec.second->Export(stage2_trajectories);
      std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "Success!" << std::endl;
  }
  catch (std::exception &ex)
  { std::cerr << "Error: " << ex.what(); }
  indicators::show_console_cursor(true);
  
  return 0;
}
