#include "CudaSolver.h"

#include "CudaManager.h"
#include "GPUParameters.h"


void CudaSolver::Solve(ICaseGenerator &generator, const IFunctional &functional, IResultFormatter &results)
{
  Case problem;
  std::vector<Case> problems;
  problems.reserve(PORTION_SIZE);
  do
  {
    problems.clear();
    while (problems.size() < PORTION_SIZE && generator.Next(problem))
    { problems.emplace_back(std::move(problem)); }

    if (!problems.empty())
    { Solve(problems, functional, results); }
  }
  while (!problems.empty());
}

void CudaSolver::Solve(const std::vector<Case> &problems, const IFunctional &functional, IResultFormatter &results)
{
  switch (Algorithm())
  {
    case NumericalAlgorithm::ONE_STEP_ADAMS:
      CudaManager<1u>(problems, Dt(), Timeout(), 256, functional, results);
      break;

    case NumericalAlgorithm::TWO_STEP_ADAMS:
      CudaManager<2u>(problems, Dt(), Timeout(), 256, functional, results);
      break;

    case NumericalAlgorithm::THREE_STEP_ADAMS:
      CudaManager<3u>(problems, Dt(), Timeout(), 256, functional, results);
      break;

    default:
      assert(false);
  }
}

void CudaSolver::Solve(const Case &problem, const IFunctional &functional, IResultFormatter &results)
{
  std::vector<Case> problems = {problem};
  Solve(problems, functional, results);
}
