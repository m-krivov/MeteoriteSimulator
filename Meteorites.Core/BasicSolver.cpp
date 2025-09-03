#include "BasicSolver.h"

void BasicSolver::Configure(NumericalAlgorithm alg, real dt, real timeout)
{
  if (alg != NumericalAlgorithm::ONE_STEP_ADAMS &&
      alg != NumericalAlgorithm::TWO_STEP_ADAMS &&
      alg != NumericalAlgorithm::THREE_STEP_ADAMS)
  { throw std::runtime_error("unsupported numerical algorithm"); }

  if (dt <= (real)0.0 || timeout < dt * 4) // 4 as an iteration of three-step Adams
  {
    throw std::runtime_error("wrong time step and/or timeout");
  }

  algorithm_   = alg;
  dt_          = dt;
  timeout_     = timeout;
}

void BasicSolver::Solve(const std::vector<Case> &problems,
                        const IFunctional &functional,
                        IResultFormatter &results)
{
  for (const auto &problem : problems)
  { ((ISolver *)this)->Solve(problem, functional, results); }
}

void BasicSolver::Solve(ICaseGenerator &generator,
                        const IFunctional &functional,
                        IResultFormatter &results)
{
  Case problem;
  while (generator.Next(problem))
  { ((ISolver *)this)->Solve(problem, functional, results); }
}
