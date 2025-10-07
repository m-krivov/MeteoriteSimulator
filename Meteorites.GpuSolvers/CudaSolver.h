#pragma once
#include "Meteorites.Core/Defs.h"
#include "Meteorites.Core/BasicSolver.h"


// Version for CUDA kernel debugging
class CudaSolver : public BasicSolver
{
  public:
    CudaSolver() = default;

    // ISolver method
    virtual void Solve(const Case &problem,
                       const IFunctional &functional,
                       IResultFormatter &results) override final;

    // ISolver method
    virtual void Solve(const std::vector<Case> &problems,
                       const IFunctional &functional,
                       IResultFormatter &results) override final;

    // ISolver method
    virtual void Solve(ICaseGenerator &generator,
                       const IFunctional &functional,
                       IResultFormatter &results) override final;
  private:
    // How many meteorites should be extracted from the generator at one time
    constexpr static size_t PORTION_SIZE = 4 * 1024;
};
