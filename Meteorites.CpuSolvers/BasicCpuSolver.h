#pragma once
#include "Meteorites.Core/Defs.h"

#include "Meteorites.Core/ISolver.h"

// Implements some common logic that is used by CPU solvers 
class BasicCpuSolver : public ISolver
{
  public:
    // ISolver method
    virtual void Configure(NumericalAlgorithm alg, real dt, real timeout) final;

  protected:
    BasicCpuSolver()
      : adams_steps_(1), dt_((real)0.001), timeout_((real)1000.0)
    { /*nothing*/ }

    size_t AdamsSteps() const { return adams_steps_; }

    real Dt() const { return dt_; }

    real Timeout() const { return timeout_; }

  private:
    size_t adams_steps_;
    real dt_;
    real timeout_;
};
