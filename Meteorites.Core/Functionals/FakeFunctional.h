#pragma once
#include "Meteorites.Core/Defs.h"

#include "IFunctional.h"

// Does nothing, just a placeholder for testing and debugging
class FakeFunctional : public IFunctional
{
  public:
    FakeFunctional() = default;

    // The member of 'IFunctional'
    virtual std::string Name() const override final
    { return "Fake"; }

    // The member of 'IFunctional'
    virtual void GetTimeStamps(size_t &num, const real *&values) const override
    {
      num = num_stub;
      values = &values_stub;
    }

    // The member of 'IFunctional'
    virtual double Compute(size_t num, const real *v, const real *h) const override
    { return std::numeric_limits<double>::max(); }

  private:
    size_t num_stub   = 1;
    real values_stub  = (real)0.0;
};
