#pragma once
#include "Meteorites.Core/Defs.h"

#include "BasicCpuSolver.h"

// Uses one-step, two-step or three-step Adams method to solve ordinary differential equations
// May be used to verify more complex solvers with some performance optimizations
class GoldSolver : public BasicCpuSolver
{
  public:
    GoldSolver() = default;

    // ISolver method
    virtual void Solve(const Case &problem,
                       const IFunctional &functional,
                       IResultFormatter &results) final;

    // ISolver method
    virtual void Solve(ICaseGenerator &generator,
                       const IFunctional &functional,
                       IResultFormatter &results) final;
};
