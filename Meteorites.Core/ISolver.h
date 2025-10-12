#pragma once
#include "Defs.h"

#include "Case.h"
#include "Functionals.h"
#include "ResultFormatters.h"
#include "ICaseGenerator.h"

// The known numerical methods that may be used to solve ODEs
enum class NumericalAlgorithm : uint32_t
{
  ONE_STEP_ADAMS,
  TWO_STEP_ADAMS,
  THREE_STEP_ADAMS
};


// Basic interface for the solvers of any type
// They uses some numerical algorithm to compute the body's trajectory and parameters
class ISolver
{
  public:
    ISolver(const ISolver &) = delete;
    ISolver &operator =(const ISolver &) = delete;
    virtual ~ISolver() = default;

    // Specializes our solver
    // This method is optional, some values are already used by default
    virtual void Configure(NumericalAlgorithm alg, real dt, real timeout) = 0;

    // Solves a single problem, computes functional value and sends solution to formatter
    virtual void Solve(const Case &problem,
                       const IFunctional &functional,
                       IResultFormatter &results) = 0;

    // Finds solutions for the fixed-size set of problems
    // This version may be parallelized
    virtual void Solve(const std::vector<Case> &problems,
                       const IFunctional &functional,
                       IResultFormatter &results) = 0;

    // Finds solutions for the problems represented by a stream
    // This version may be parallelized
    virtual void Solve(ICaseGenerator &generator,
                       const IFunctional &functional,
                       IResultFormatter &results) = 0;

  protected:
    ISolver() = default;
};
