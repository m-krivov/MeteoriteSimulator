#include "BasicSolver.h"


void BasicSolver::Configure(NumericalAlgorithm alg, real dt, real timeout)
{
  size_t steps{};
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

  if (dt <= (real)0.0 || timeout < dt * (steps + 1))
  {
    throw std::runtime_error("wrong time step and/or timeout");
  }

  algorithm_   = alg;
  steps_       = steps;
  dt_          = dt;
  timeout_     = timeout;
}

void BasicSolver::Solve(const std::vector<VirtualMeteoroid> &problems,
                        const IFunctional &functional,
                        ISimulationRecorder &results)
{
  for (const auto &problem : problems)
  { ((ISolver *)this)->Solve(problem, functional, results); }
}

void BasicSolver::Solve(IMeteoroidGenerator &generator,
                        const IFunctional &functional,
                        ISimulationRecorder &results)
{
  VirtualMeteoroid problem;
  while (generator.MoveNext())
  { ((ISolver *)this)->Solve(generator.Current(), functional, results); }
}
