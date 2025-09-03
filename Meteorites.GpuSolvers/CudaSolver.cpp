#include "CudaSolver.h"

#include "CudaManager.h"
#include "GPUParameters.h"

/*void CudaSolver::Configure(NumericalAlgorithm alg, real dt, real timeout)
{
  size_t steps = 0;
  switch (alg) {
  case NumericalAlgorithm::ONE_STEP_ADAMS:
    steps = 1;
    break;

  case NumericalAlgorithm::TWO_STEP_ADAMS:
    steps = 2;
    break;

  case NumericalAlgorithm::THREE_STEP_ADAMS:
    steps = 3;
    break;

  default:
    throw std::runtime_error("unsupported numerical algorithm");
  }

  if (dt <= (real)0.0 || timeout < dt * (steps + 1)) {
    throw std::runtime_error("wrong time step and/or timeout");
  }

  adams_steps_ = steps;
  dt_ = dt;
  timeout_ = timeout;
}*/

void CudaSolver::Solve(ICaseGenerator& generator, const IFunctional& functional, IResultFormatter& results)
{
  std::vector<Case> problems;
  Case problem;
  //while (generator.Next(problem)) {
  for (unsigned int i = 0; i < CASE_NUM; i++) {
    generator.Next(problem);
    problems.emplace_back(std::move(problem));
  }
  Solve(problems, functional, results);
}

void CudaSolver::Solve(const std::vector<Case> &problems, const IFunctional &functional, IResultFormatter &results)
{
  switch (Algorithm())
  {
    case NumericalAlgorithm::ONE_STEP_ADAMS:
      CudaManager<1u, ITERS_PER_KERNEL>(problems, functional, Dt(), Timeout(), results);
      break;

    case NumericalAlgorithm::TWO_STEP_ADAMS:
      CudaManager<2u, ITERS_PER_KERNEL>(problems, functional, Dt(), Timeout(), results);
      break;

    case NumericalAlgorithm::THREE_STEP_ADAMS:
      CudaManager<3u, ITERS_PER_KERNEL>(problems, functional, Dt(), Timeout(), results);
      break;

    default:
      assert(false);
  }
}

void CudaSolver::Solve(const Case &problem, const IFunctional& functional, IResultFormatter& results)
{
  std::vector<Case> problems = {problem};
  Solve(problems, functional, results);
}
