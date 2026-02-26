#pragma once
#include "Meteorites.Core/Defs.h"

#include "BasicFunctional.h"
#include "Meteorites.Core/IMeteorite.h"

// Computes the L1 norm (Manhattan distance) of the difference between observed and predicted trajectories
class L1Functional : public BasicFunctional
{
  public:
    L1Functional(const std::shared_ptr<const IMeteorite> &meteorite,
                 real lambda_v = (real)1.0,
                 real lambda_h = (real)1.0)
      : BasicFunctional(meteorite, lambda_v, lambda_h)
    { }

    // Constructor with custom weights for each measurement point
    L1Functional(const std::shared_ptr<const IMeteorite> &meteorite,
                 real lambda_v,
                 real lambda_h,
                 const std::vector<real> &weights)
      : BasicFunctional(meteorite, lambda_v, lambda_h, weights)
    { }

    // The member of 'IFunctional'
    virtual std::string Name() const override final
    { return HasWeights() ? "L1-Weighted" : "L1"; }

    // The member of 'IFunctional'
    virtual double Compute(size_t num, const real *v, const real *h) const override final;
};
