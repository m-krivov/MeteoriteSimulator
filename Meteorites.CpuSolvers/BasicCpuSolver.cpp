#include "BasicCpuSolver.h"

void BasicCpuSolver::Configure(NumericalAlgorithm alg, real dt, real timeout)
{
  size_t steps = 0;
  switch (alg)
  {
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

  adams_steps_ = steps;
  dt_          = dt;
  timeout_     = timeout;
}
