#pragma once
#include "Meteorites.Core/Defs.h"

#include "BasicFunctional.h"
#include "Meteorites.Core/IMeteorite.h"

// Computes the C functional (maximum error) between observed and predicted trajectories
class CFunctional : public BasicFunctional
{
  public:
    CFunctional(const IMeteorite &meteorite,
                real lambda_v = (real)1.0,
                real lambda_h = (real)1.0)
      : BasicFunctional(meteorite, lambda_v, lambda_h)
    { }

    // The member of 'IFunctional'
    virtual std::string Name() const override final
    { return "C"; }

    // The member of 'IFunctional'
    virtual double Compute(size_t num, const real *v, const real *h) const override final;
};
