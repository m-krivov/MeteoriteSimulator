#pragma once
#include "Meteorites.Core/Defs.h"

#include "Meteorites.Core/Solvers/BasicSolver.h"

// This solver implements the same algorithm as GoldSolver but leverages OpenMP parallelization
class OpenmpSolver : public BasicSolver
{
  public:
    OpenmpSolver(size_t cores = 0);

  protected:
    // BasicSolver method
    virtual size_t BatchSize() const override final
    { assert(cores_ > 0); return cores_; }

    // BasicSolver method
    virtual void SolveAny(MeteoroidEnumerator &problems,
                          const IFunctional &functional,
                          ISimulationRecorder &results) override final;

  private:
    size_t cores_{1};
};
