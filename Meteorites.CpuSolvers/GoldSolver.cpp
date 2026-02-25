#include "GoldSolver.h"

#include "GoldAdamsMethod.h"

void GoldSolver::SolveAny(MeteoroidEnumerator &problems,
                          const IFunctional &functional,
                          ISimulationRecorder &results)
{
  auto method = CreateGoldAdamsMethod(Algorithm());

  size_t size = 0;
  const VirtualMeteoroid *meteoroid = nullptr;
  while (problems.MoveNext(meteoroid, size))
  {
    assert(size == 1);
    assert(meteoroid != nullptr);
    method(*meteoroid, Dt(), Timeout(), functional, results);
  }
}
