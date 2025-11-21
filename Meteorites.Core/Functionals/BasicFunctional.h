#pragma once
#include "Meteorites.Core/Defs.h"

#include "IFunctional.h"
#include "Meteorites.Core/IMeteorite.h"

// Basic class for the classical functionals like C, L1 or L2
class BasicFunctional : public IFunctional
{
  public:

    // The member of 'IFunctional'
    virtual void GetTimeStamps(size_t &num, const real *&values) const override final;

  protected:
    BasicFunctional(const IMeteorite &meteorite,
                    real lambda_v, real lambda_h);

    real LambdaV() const { return lambda_v_; }

    real LambdaH() const { return lambda_h_; }

    const std::vector<real> &Time() const { return time_; }

    const std::vector<real> &Velocity() const { return v_; }

    const std::vector<real> &Height() const { return h_; }

  private:
    real lambda_v_ = (real)1.0;
    real lambda_h_ = (real)1.0;
    std::vector<real> time_, v_, h_;
};
