#include "Constants.h"

#define PI (real)3.1415926535897932384626433832795

real Constants::rho_a(real h)
{
  constexpr real rho_0 = (real)1.125;
  constexpr real alpha = (real)(0.029 * 9.8 / (8.31 * 273));
  return rho_0 * std::exp(-alpha * h);
}

real Constants::Midsection(real M, real rho)
{
  static const real coeff = (real)(PI * std::pow(3.0 / (4.0 * PI), 2.0 / 3.0));
  return coeff * std::pow(M / rho, (real)(2.0 / 3.0));
}
