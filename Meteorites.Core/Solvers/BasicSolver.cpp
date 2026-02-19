#include "BasicSolver.h"

//---------------------------
//--- MeteoroidEnumerator ---
//---------------------------

BasicSolver::MeteoroidEnumerator::MeteoroidEnumerator(const VirtualMeteoroid &meteoroid)
{
  impl_ = [&meteoroid, first = true]
          (const VirtualMeteoroid *&batch, size_t &size) mutable -> bool
  {
    if (first)
    {
      batch = &meteoroid;
      size = 1;
      first = false;
      return true;
    }
    else
    {
      batch = nullptr;
      size = 0;
      return false;
    }
  };
}

BasicSolver::MeteoroidEnumerator::MeteoroidEnumerator(const std::vector<VirtualMeteoroid> &meteoroids,
                                                      size_t batch_size)
{
  assert(batch_size >= 1);
  impl_ = [&meteoroids, batch_size, next = (size_t)0]
          (const VirtualMeteoroid *&batch, size_t &size) mutable -> bool
  {
    if (next < meteoroids.size())
    {
      batch = meteoroids.data() + next;
      size = std::min(batch_size, meteoroids.size() - next);
      next += size;
      return true;
    }
    else
    {
      batch = nullptr;
      size = 0;
      return false;
    }
  };
}

BasicSolver::MeteoroidEnumerator::MeteoroidEnumerator(IMeteoroidGenerator &gen,
                                                      size_t batch_size)
{
  assert(batch_size >= 1);
  impl_ = [&gen, buffer = std::vector<VirtualMeteoroid>(batch_size)]
          (const VirtualMeteoroid *&batch, size_t &size) mutable -> bool
  {
    size = 0;
    while (size < buffer.size() && gen.MoveNext())
    {
      buffer[size++] = gen.Current();
    }
    batch = size > 0 ? buffer.data() : nullptr;
    return size > 0;
  };
}

bool BasicSolver::MeteoroidEnumerator::MoveNext(const VirtualMeteoroid *&batch, size_t &size)
{
  return impl_(batch, size);
}

//-------------------
//--- BasicSolver ---
//-------------------

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

void BasicSolver::Solve(const VirtualMeteoroid &problem,
                        const IFunctional &functional,
                        ISimulationRecorder &results)
{
  MeteoroidEnumerator en(problem);
  SolveAny(en, functional, results);
}

void BasicSolver::Solve(const std::vector<VirtualMeteoroid> &problems,
                        const IFunctional &functional,
                        ISimulationRecorder &results)
{
  MeteoroidEnumerator en(problems, BatchSize());
  SolveAny(en, functional, results);
}

void BasicSolver::Solve(IMeteoroidGenerator &generator,
                        const IFunctional &functional,
                        ISimulationRecorder &results)
{
  MeteoroidEnumerator en(generator, BatchSize());
  SolveAny(en, functional, results);
}
