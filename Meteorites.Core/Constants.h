#pragma once
#include "Defs.h"

#if defined(__CUDA_ARCH__)
  #define CONSTANTS_DEVICE __device__
  #define CONSTANTS_STATIC
#else
  #define CONSTANTS_DEVICE
  #define CONSTANTS_STATIC static
#endif

// Defines pseudo-constant values that are used in the model
// Where:
//    R       - radius of the planet, Earth
//    g       - acceleration of gravity, Earth
//    RhoAtm  - atmosphere density, Earth

namespace Constants
{

static inline constexpr real R()
{ return (real)6371000; }

static inline constexpr real g(real h)
{
  real r_km = (6371 + h / 1000);
  return ((real)(6.67428 * 5.9726) / (r_km * r_km) * (real)1e7);
}

static inline constexpr real g() { return g(0); }

CONSTANTS_DEVICE
static inline real RhoAtm(real h)
{
  constexpr real rho_0 = (real)1.125;
  constexpr real alpha = (real)(0.029 * 9.8 / (8.31 * 273));
  return rho_0 * std::exp(-alpha * h);
}

CONSTANTS_DEVICE
static real Midsection(real M, real rho)
{
  CONSTANTS_STATIC const real coeff = (real)(M_PI * std::pow(3.0 / (4.0 * M_PI), 2.0 / 3.0));
  return coeff * std::pow(M / rho, (real)(2.0 / 3.0));
}

} // namespace Constants

#undef CONSTANTS_DEVICE
#undef CONSTANTS_STATIC
