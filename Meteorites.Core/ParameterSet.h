#pragma once
#include "Defs.h"

// Defines ranges for possible values of the unknown parameters
// We want to vary them and consider a different virtual meteorites
// So it is important to correctly define such boundaries
class ParameterSet
{
  public:
    ParameterSet()
      : H_(std::make_pair((real)0.0, (real)0.0)),
        ch_(std::make_pair((real)0.0, (real)0.0)),
        rho_(std::make_pair((real)0.0, (real)0.0)),
        cd_(std::make_pair((real)0.0, (real)0.0)),
        cl_(std::make_pair((real)0.0, (real)0.0)),
        m0_(std::make_pair((real)0.0, (real)0.0)),
        gamma0_(std::make_pair((real)0.0, (real)0.0))
    { }

    ParameterSet(const std::pair<real, real> &H,
                 const std::pair<real, real> &ch,
                 const std::pair<real, real> &rho,
                 const std::pair<real, real> &cd,
                 const std::pair<real, real> &cl,
                 const std::pair<real, real> &m0,
                 const std::pair<real, real> &gamma0)
      : H_(H), ch_(ch), rho_(rho), cd_(cd), cl_(cl), m0_(m0), gamma0_(gamma0)
    { }

    ParameterSet(const ParameterSet &) = default;
    ParameterSet &operator =(const ParameterSet &) = default;

    const std::pair<real, real> &H() const { return H_; }
    const std::pair<real, real> &Ch() const { return ch_; }
    const std::pair<real, real> &Rho() const { return rho_; }
    const std::pair<real, real> &Cd() const { return cd_; }
    const std::pair<real, real> &Cl() const { return cl_; }
    const std::pair<real, real> &M0() const { return m0_; }
    const std::pair<real, real> &Gamma0() const { return gamma0_; }

  private:
    std::pair<real, real> H_, ch_, rho_;
    std::pair<real, real> cd_, cl_;
    std::pair<real, real> m0_, gamma0_;
};
