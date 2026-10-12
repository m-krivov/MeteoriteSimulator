#include "OpenmpSolver.h"

#include <omp.h>
#include "GoldAdamsMethod.h"
#include "Meteorites.Core/Recorders/MetaRecorder.h"
#include "Meteorites.Core/Recorders/BufferingRecorder.h"

OpenmpSolver::OpenmpSolver(size_t cores)
{
  cores_ = cores;
  if (cores_ == 0)
  { cores_ = (size_t)std::max(1, omp_get_num_procs() - 1); }
}

void OpenmpSolver::SolveAny(MeteoroidEnumerator &problems,
                            const IFunctional &functional,
                            ISimulationRecorder &results)
{
  auto method = CreateGoldAdamsMethod(Algorithm());

  std::vector<std::unique_ptr<ISimulationRecorder>> recorders(BatchSize());
  std::function<void()> commit_and_reset;

  // If trajectory is needed, we must cache all trajectories of a batch
  if (results.NeedTrajectory())
  {
    std::vector<BufferingRecorder *> refs(recorders.size());
    for (size_t i = 0; i < recorders.size(); i++)
    {
      recorders[i] = std::unique_ptr<ISimulationRecorder>(refs[i] = new BufferingRecorder((real)0.0));
    }
    commit_and_reset = [refs = std::move(refs), &results]() -> void {
      for (auto ref : refs)
      {
        assert(ref->Trajectories().size() <= 1);
        if (ref->Trajectories().size() == 1)
        {
          ref->Trajectories().front().ExportTo(results);
          ref->Reset();
        }
      }
    };
  }
  // Otherwise, we can process meteoroids fully in parallel
  else
  {
    std::vector<MetaRecorder *> refs(recorders.size());
    for (size_t i = 0; i < recorders.size(); i++)
    {
      recorders[i] = std::unique_ptr<ISimulationRecorder>(refs[i] = new MetaRecorder(1, 1));
    }
    commit_and_reset = [refs = std::move(refs), &results]() -> void {
      for (auto ref : refs)
      {
        std::vector<MeteoroidSummary> summary;
        ref->MoveTo(summary);
        assert(summary.size() <= 1);

        if (summary.size() == 1)
        { summary.front().ExportTo(results); }
      }
    };
  }
  
  // Now, let's process all meteoroids in parallel
  size_t size = 0;
  const VirtualMeteoroid *meteoroid = nullptr;
  while (problems.MoveNext(meteoroid, size))
  {
    assert(size <= recorders.size());
    assert(meteoroid != nullptr);

    auto cores = std::min(cores_, size);
    #pragma omp parallel for num_threads(cores)
    for (int i = 0; i < size; i++)
    { method(meteoroid[i], Dt(), Timeout(), functional, *recorders[i]); }
    
    commit_and_reset();
  }
}
