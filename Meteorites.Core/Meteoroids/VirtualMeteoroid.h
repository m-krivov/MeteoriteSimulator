#pragma once
#include "Meteorites.Core/Defs.h"

// Contains the complete information about a single virtual meteoroid
// Implemented as 'public' to simplify CUDA kernels
struct VirtualMeteoroid
{
  public:
    real H                     = (real)0.0; // Enthalpy of destruction
    real Ch                    = (real)0.0; // Coefficient of heat exchange with the environment
    real Rho                   = (real)0.0; // Average density of the body
    real Cd                    = (real)0.0; // Aerodynamic drag coefficient
    real Cl                    = (real)0.0; // Aerodynamic lift coefficient
    real M0                    = (real)0.0; // Initial mass of the body
    real V0                    = (real)0.0; // Initial velocity of the body
    real h0                    = (real)0.0; // Initial height of the body
    constexpr static real l0   = (real)0.0; // Stub for path length, always equal to zero
    real Gamma0                = (real)0.0; // Initial angle between trajectory and horizon line

    DEVICE VirtualMeteoroid() = default;
    DEVICE VirtualMeteoroid(real H_, real ch_, real rho_, real cd_, real cl_,
                            real m0_, real v0_, real h0_, real gamma0_)
      : H(H_), Ch(ch_), Rho(rho_), Cd(cd_), Cl(cl_),
        M0(m0_), V0(v0_), h0(h0_), Gamma0(gamma0_) { }
    DEVICE VirtualMeteoroid(const VirtualMeteoroid &) = default;
    DEVICE VirtualMeteoroid &operator =(const VirtualMeteoroid &) = default;
};
