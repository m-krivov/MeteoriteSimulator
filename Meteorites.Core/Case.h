#pragma once

#include "Defs.h"

// Containes the complete information about a single virtual case
// Where:
//    M0, V0, h0 - the initial mass, velocity and height of the body
//    Gamma0     - the initial angle between trajectory and horizon line
//
//    H          - enthalpy of destruction
//    ch         - coefficient of heat exchange with the environment
//    Rho        - average density of the body
//    cl, cd     - the aerodynamic coefficients
class Case
{
  public:
    Case()
      : H_((real)0.0f), ch_((real)0.0f), rho_((real)0.0f),
        cd_((real)0.0f), cl_((real)0.0f),
        m0_((real)0.0f), v0_((real)0.0f),
        h0_((real)0.0f), gamma0_((real)0.0f) { }

    Case(real H, real ch, real rho, real cd, real cl,
         real m0, real v0, real h0, real gamma0)
      : H_(H), ch_(ch), rho_(rho), cd_(cd), cl_(cl),
        m0_(m0), v0_(v0), h0_(h0), gamma0_(gamma0) { }

    Case(const Case &) = default;
    Case &operator =(const Case &) = default;

    real H() const { return H_; }
    real Ch() const { return ch_; }
    real Rho() const { return rho_; }
    real Cd() const { return cd_; }
    real Cl() const { return cl_; }

    real M0() const { return m0_; }
    real V0() const { return v0_; }
    real h0() const { return h0_; }
    real l0() const { return (real)0.0; }
    real Gamma0() const { return gamma0_; }

  private:
    real H_, ch_, rho_, cd_, cl_;
    real m0_, v0_, h0_, gamma0_;
};
