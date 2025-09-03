#pragma once
#include "Meteorites.Core/Defs.h"
#include "Meteorites.Core/BasicSolver.h"


// Implements some common logic that is used by CPU solvers
class CudaSolver : public BasicSolver
{
  public:
    CudaSolver() = default;

    // ISolver method
    virtual void Solve(const Case& problem, const IFunctional& functional, IResultFormatter& results) override final;

    // ISolver method
    virtual void Solve(const std::vector<Case> &problems, const IFunctional &functional, IResultFormatter &results) override final;

    // ISolver method
    virtual void Solve(ICaseGenerator& generator, const IFunctional& functional, IResultFormatter& results) override final;
};
