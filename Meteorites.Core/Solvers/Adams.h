#pragma once
#include "Meteorites.Core/Defs.h"

#include "Meteorites.Core/Meteoroids/VirtualMeteoroid.h"
#include "Meteorites.Core/Constants.h"

namespace Adams
{

// Caches meteoroid's parameters and constants that may be considered as unchangeable values
// Also precomputes some frequently used expressions
struct Unchangeable
{
  // Meteoroid's parameters
  real H          = (real)0.0;
  real Ch         = (real)0.0;
  real Cd         = (real)0.0;
  real Cl         = (real)0.0;
  real Rho        = (real)0.0;

  // Constants
  real R          = (real)0.0;

  // Precomputed values
  real inv_H      = (real)0.0;
  real inv_R      = (real)0.0;
  real midsection = (real)0.0;

  Unchangeable() = default;
  DEVICE Unchangeable(const VirtualMeteoroid &meteoroid)
    : H{ meteoroid.H },
      Ch{ meteoroid.Ch },
      Cd{ meteoroid.Cd },
      Cl{ meteoroid.Cl },
      Rho{ meteoroid.Rho },
      R{ Constants::R() }
  {
    inv_H       = (real)1.0 / H;
    inv_R       = (real)1.0 / R;
    midsection  = Constants::Midsection((real)1.0, Rho);
  }
  Unchangeable(const Unchangeable &) = default;
  Unchangeable &operator =(const Unchangeable &) = default;
};


// Represents four primary ODE's variables (velocity, angle, height and mass)
// and one auxiliary (length) at some point in time
struct State
{
  real V       = (real)0.0;
  real Gamma   = (real)0.0;
  real h       = (real)0.0;
  real l       = (real)0.0;
  real M       = (real)0.0;

  State() = default;
  State(const State &) = default;
  State &operator =(const State &) = default;
};


// Computes the right part of Stulov's ODE system
// These values are equal to derivatives
inline DEVICE
void ComputeDerivatives(State &dx, const State &x, const Unchangeable &params)
{
  assert(params.Rho > 0.0f);
  assert(params.H > 1e-3f);
  assert(x.V > 0.0f);

  // Branchless guard for non-positive mass for cases when dt is too large
  const bool alive = x.M > (real)0.0;
  const real M     = std::max(x.M, (real)1e-30);

  auto sin_gamma  = std::sin(x.Gamma);
  auto cos_gamma  = std::cos(x.Gamma);
  auto g          = Constants::g(x.h);
  auto rho_a      = Constants::RhoAtm(x.h);
  auto midsection = params.midsection * std::pow(M, (real)(2.0 / 3.0));
  auto inv_2M     = (real)0.5 / M;

  const real dx_V = - params.Cd * rho_a * x.V * x.V * midsection * inv_2M
                    + g * sin_gamma;
  const real dx_G = + g * cos_gamma / x.V
                    - x.V * cos_gamma * params.inv_R
                    - params.Cl * rho_a * x.V * midsection * inv_2M;
  const real dx_h = - x.V * sin_gamma;
  const real dx_l = x.V * (params.R / (params.R + x.h)) * cos_gamma;
  const real dx_M = - (params.Ch * rho_a * x.V * x.V * x.V
                       * midsection * (real)0.5) * params.inv_H;

  dx.V     = alive ? dx_V : (real)0.0;
  dx.Gamma = alive ? dx_G : (real)0.0;
  dx.h     = alive ? dx_h : (real)0.0;
  dx.l     = alive ? dx_l : (real)0.0;
  dx.M     = alive ? dx_M : (real)0.0;
}


// Iteration for one-step Adam's method
inline DEVICE
void OneStepIteration(State &x1, State &dx0, const State &x0,
                      const Unchangeable &params, real dt)
{
  ComputeDerivatives(dx0, x0, params);

  x1 = { x0.V     + dx0.V     * dt,
         x0.Gamma + dx0.Gamma * dt,
         x0.h     + dx0.h     * dt,
         x0.l     + dx0.l     * dt,
         x0.M     + dx0.M     * dt };
}

// Iteration for two-step Adam's method
inline DEVICE
void TwoStepIteration(State &x2, State &dx1, const State &x1, const State &dx0,
                      const Unchangeable &params, real dt)
{
  constexpr auto c1 =  (real)1.5;
  constexpr auto c0 = -(real)0.5;

  ComputeDerivatives(dx1, x1, params);

  x2 = { x1.V     + ( c1 * dx1.V     + c0 * dx0.V     ) * dt,
         x1.Gamma + ( c1 * dx1.Gamma + c0 * dx0.Gamma ) * dt,
         x1.h     + ( c1 * dx1.h     + c0 * dx0.h     ) * dt,
         x1.l     + ( c1 * dx1.l     + c0 * dx0.l     ) * dt,
         x1.M     + ( c1 * dx1.M     + c0 * dx0.M     ) * dt };
}

// Iteration for three-step Adam's method
inline DEVICE
void ThreeStepIteration(State &x3, State &dx2, const State &x2, const State &dx1, const State &dx0,
                        const Unchangeable &params, real dt)
{
  constexpr auto c2 =  (real)23 / 12;
  constexpr auto c1 = -(real)16 / 12;
  constexpr auto c0 =  (real)5  / 12;

  ComputeDerivatives(dx2, x2, params);
  
  x3 = { x2.V     + ( c2 * dx2.V     + c1 * dx1.V     + c0 * dx0.V     ) * dt,
         x2.Gamma + ( c2 * dx2.Gamma + c1 * dx1.Gamma + c0 * dx0.Gamma ) * dt,
         x2.h     + ( c2 * dx2.h     + c1 * dx1.h     + c0 * dx0.h     ) * dt,
         x2.l     + ( c2 * dx2.l     + c1 * dx1.l     + c0 * dx0.l     ) * dt,
         x2.M     + ( c2 * dx2.M     + c1 * dx1.M     + c0 * dx0.M     ) * dt };
}

// Performs uni-step iteration of Adams' method using a cycled buffer with derivatives
template <unsigned int STEPS>
inline DEVICE
void Iteration(State &x, State dx[STEPS], const Unchangeable &params, size_t nxt, real dt);

template<>
inline DEVICE
void Iteration<1>(State &x, State dx[1], const Unchangeable &params, size_t nxt, real dt)
{ OneStepIteration(x, dx[nxt], x, params, dt); }

template<>
inline DEVICE
void Iteration<2>(State &x, State dx[2], const Unchangeable &params, size_t nxt, real dt)
{ TwoStepIteration(x, dx[nxt], x, dx[(nxt + 1) % 2], params, dt); }

template<>
inline DEVICE
void Iteration<3>(State &x, State dx[3], const Unchangeable &params, size_t nxt, real dt)
{ ThreeStepIteration(x, dx[nxt], x, dx[(nxt + 2) % 3], dx[(nxt + 1) % 3], params, dt); }

} // namespace Adams
