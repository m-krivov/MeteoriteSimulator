#pragma once
#include "Meteorites.Core/Defs.h"

#include "Meteorites.Core/Solvers/BasicSolver.h"


// Uses one-step, two-step or three-step Adams method to solve ordinary differential equations
// May be used to verify more complex solvers with some performance optimizations
class GoldSolver : public BasicSolver
{
  public:
    GoldSolver() = default;

  protected:
    // BasicSolver method
    virtual size_t BatchSize() const override final { return 1; }

    // BasicSolver method
    virtual void SolveAny(MeteoroidEnumerator &problems,
                          const IFunctional &functional,
                          ISimulationRecorder &results) override final;
};
