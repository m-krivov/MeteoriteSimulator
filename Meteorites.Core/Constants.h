#pragma once
#include "Defs.h"

// Defines pseudo-constant values that are used in the model
// Where:
//    R      - radius of the planet, Earth
//    g      - acceleration of gravity, Earth
//    rho_a  - atmosphere density, Earth
class Constants
{
  public:
    Constants() = delete;
    Constants(const Constants &) = delete;
    void operator =(const Constants &) = delete;

    static constexpr real R()
    { return (real)6371000; }

    static constexpr real g(real h)
    {
      real r_km = (6371 + h / 1000);
      return ((real)(6.67428 * 5.9726) / (r_km * r_km) * (real)1e7);
    }

    static constexpr real g() { return g(0); }

    static real rho_a(real h);
    static real Midsection(real M, real rho);
};
