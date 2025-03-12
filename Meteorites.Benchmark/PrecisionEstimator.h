#pragma once
#include "Meteorites.Core/Case.h"
#include "Meteorites.CpuSolvers/GoldSolver.h"

// Simulates some virtual meteorite and considers time point 't_end * 0.9' to avoid division by zero
// After that, prints relative numerical errors for 'm', 'V', 'h' and 'Gamma' as CSV table
class PrecisionEstimator
{
  public:
    PrecisionEstimator() = default;
    PrecisionEstimator(const PrecisionEstimator &) = delete;
    PrecisionEstimator &operator =(const PrecisionEstimator &) = delete;

    // Reconstructs trajectory of a virtual meteorite using 1-step, 2-step and 3-step methods
    void CompareMethods(const Case &problem, real dt, std::ostream &str);

    // Estimates the effect of 'dt' on the numerical error
    void CompareSteps(const Case &problem, NumericalAlgorithm method, std::ostream &str);

    // Estimates the effect of peturbation for each initial parameter
    void ComparePerturbations(const Case &problem, NumericalAlgorithm method, real dt, std::ostream &str);
};
