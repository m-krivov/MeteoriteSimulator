#pragma once
#include "Defs.h"

#include "Case.h"
#include "Functionals.h"
#include "ResultFormatters.h"
#include "ICaseGenerator.h"

// Basic interface for the solvers of any type
// They uses some numerical algorithm to compute the body's trajectory and parameters
class ISolver
{
  public:
    ISolver(const ISolver &) = delete;
    ISolver &operator =(const ISolver &) = delete;

    // Specializes our solver
    // This method is optional, some values are already used by default
    virtual void Configure(real dt, real timeout) = 0;

    // Solves single problem, computes functional value and sends solution to formatter
    virtual void Solve(const Case &problem,
                       const IFunctional &functional,
                       IResultFormatter &results) = 0;

    // Finds solutions for the set of problems
    // This version may be parallelized
    virtual void Solve(ICaseGenerator &generator,
                       const IFunctional &functional,
                       IResultFormatter &results) = 0;

  protected:
    ISolver() = default;
};
