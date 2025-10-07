#pragma once
#include "Meteorites.Core/Defs.h"
#include "Meteorites.Core/ISolver.h"

// Implements some common logic
class BasicSolver : public ISolver
{
  public:
    // ISolver method
    virtual void Configure(NumericalAlgorithm alg, real dt, real timeout) override final;

    // ISolver method
    virtual void Solve(const std::vector<Case> &problems,
                       const IFunctional &functional,
                       IResultFormatter &results) override;

    // ISolver method
    virtual void Solve(ICaseGenerator &generator,
                       const IFunctional &functional,
                       IResultFormatter &results) override;

  protected:
    BasicSolver() = default;

    NumericalAlgorithm Algorithm() const { return algorithm_; }

    size_t Steps() const { return steps_; }

    real Dt() const { return dt_; }

    real Timeout() const { return timeout_; }

  private:
    NumericalAlgorithm algorithm_ = NumericalAlgorithm::ONE_STEP_ADAMS;
    size_t steps_ = 1;
    real dt_ = (real)0.001;
    real timeout_ = (real)1000.0;
};
