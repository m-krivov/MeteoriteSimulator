#pragma once

#include "Defs.h"

// Containes the complete information about a single virtual case
// Where:
//    M0, V0, h0 - the initial mass, velocity and height of the body
//    l0         - just a stub for 'path length', always equal to zero
//    Gamma0     - the initial angle between trajectory and horizon line
//
//    H          - enthalpy of destruction
//    ch         - coefficient of heat exchange with the environment
//    Rho        - average density of the body
//    cl, cd     - the aerodynamic coefficients
// Implemented as 'public' to simplify CUDA kernels
struct Case
{
  public:
    real H                     = (real)0.0;
    real Ch                    = (real)0.0;
    real Rho                   = (real)0.0;
    real Cd                    = (real)0.0;
    real Cl                    = (real)0.0;

    real M0                    = (real)0.0;
    real V0                    = (real)0.0;
    real h0                    = (real)0.0;
    constexpr static real l0   = (real)0.0;
    real Gamma0                = (real)0.0;

    Case() = default;
    Case(real H_, real ch_, real rho_, real cd_, real cl_,
         real m0_, real v0_, real h0_, real gamma0_)
      : H(H_), Ch(ch_), Rho(rho_), Cd(cd_), Cl(cl_),
        M0(m0_), V0(v0_), h0(h0_), Gamma0(gamma0_) { }
    Case(const Case &) = default;
    Case &operator =(const Case &) = default;
};
