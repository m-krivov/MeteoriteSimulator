#pragma once
#include "Meteorites.Core/Defs.h"

#include "Meteorites.Core/ISolver.h"

// Uses one-step, two-step or three-step Adams method to solve ordinary differential equations
// May be used to verify more complex solvers with some performance optimizations
class GoldSolver : public ISolver
{
  public:
    enum Method
    {
      ONE_STEP_ADAMS   = 1,
      TWO_STEP_ADAMS   = 2,
      THREE_STEP_ADAMS = 3,
    };

    GoldSolver(Method method);

    // ISolver method
    virtual void Configure(real dt, real timeout) override;

    // ISolver method
    virtual void Solve(const Case &problem,
                       const IFunctional &functional,
                       IResultFormatter &results) override;

    // ISolver method
    virtual void Solve(ICaseGenerator &generator,
                       const IFunctional &functional,
                       IResultFormatter &results) override;

  private:
    size_t adams_steps_;
    real dt_;
    real timeout_;
};
