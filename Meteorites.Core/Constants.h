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

    static real R();
    static real g(real h);
    static real g() { return g(0); }
    static real rho_a(real h);
    static real Midsection(real M, real rho);
};
