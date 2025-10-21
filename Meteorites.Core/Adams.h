#pragma once
#include "Meteorites.Core/Defs.h"

#include "Meteorites.Core/Case.h"
#include "Meteorites.Core/Constants.h"

namespace Adams
{

// Caches constants and problem's parameters that may be considered as unchangeable values
// This class is a thin wrapper for 'Case' and 'Constants'
struct Unchangeable
{
  real H   = (real)0.0;
  real Ch  = (real)0.0;
  real Cd  = (real)0.0;
  real Cl  = (real)0.0;
  real Rho = (real)0.0;
  real R   = (real)0.0;

  Unchangeable() = default;
  DEVICE Unchangeable(const Case &problem)
    : H{ problem.H },
      Ch{ problem.Ch },
      Cd{ problem.Cd },
      Cl{ problem.Cl },
      Rho{ problem.Rho },
      R{ Constants::R() }
  { }
  Unchangeable(const Unchangeable &) = default;
  Unchangeable &operator =(const Unchangeable &) = default;
};


// Represents four primary parameters (velocity, angle, height and mass)
// and one auxiliary (length) at some point in time
// Also contains precomputed right part of Stulov's ODE
struct Layer
{
  real V       = (real)0.0;
  real Gamma   = (real)0.0;
  real h       = (real)0.0;
  real l       = (real)0.0;
  real M       = (real)0.0;
    
  real fV      = (real)0.0;
  real fGamma  = (real)0.0;
  real fh      = (real)0.0;
  real fl      = (real)0.0;
  real fM      = (real)0.0;

  Layer() = default;
  Layer(const Layer &) = default;
  Layer &operator =(const Layer &) = default;
};


// Computes the right part of Stulov's ODE and simply assign the values
static inline DEVICE
void SetLayer(Layer &layer, const Unchangeable &params,
              real V, real Gamma, real h, real l, real M)
{
  assert(params.Rho > 0.0f);
  assert(params.H > 1e-3f);
  assert(V > 0.0f);

  layer.V     = V;
  layer.Gamma = Gamma;
  layer.h     = h;
  layer.l     = l;
  layer.M     = M;

  if (M <= (real)0.0)   // probably, 'dt' is too large
  {
    layer.fV = layer.fh = layer.fl = layer.fM = (real)0.0;
  }
  else
  {
    auto sin_gamma = std::sin(Gamma);
    auto cos_gamma = std::cos(Gamma);
    auto g = Constants::g(h);
    auto rho_a = Constants::RhoAtm(h);
    auto midsection = Constants::Midsection(M, params.Rho);

    layer.fV = - params.Cd * rho_a * V * V * midsection / (2 * M)
               + g * sin_gamma;
    layer.fGamma =  + g * cos_gamma / V
                    - V * cos_gamma / params.R
                    - params.Cl * rho_a * V * midsection / (2 * M);
    layer.fh = - V * sin_gamma;
    layer.fl = V * (params.R / (params.R + h)) * cos_gamma;
    layer.fM = - (params.Ch * rho_a * V * V * V * midsection / 2) / params.H;
  }
}


// Iteration for one-step Adam's method
static inline DEVICE
void OneStepIteration(Layer &res, const Layer &f0,
                      const Unchangeable &params, real dt)
{
  SetLayer(res, params,
           f0.V     + f0.fV     * dt,
           f0.Gamma + f0.fGamma * dt,
           f0.h     + f0.fh     * dt,
           f0.l     + f0.fl     * dt,
           f0.M     + f0.fM     * dt);
}

// Iteration for two-step Adam's method
static inline DEVICE
void TwoStepIteration(Layer &res, const Layer &f1, const Layer &f0,
                      const Unchangeable &params, real dt)
{
  const auto c1 =  (real)1.5;
  const auto c0 = -(real)0.5;
  
  SetLayer(res, params,
           f1.V     + ( c1 * f1.fV     + c0 * f0.fV     ) * dt,
           f1.Gamma + ( c1 * f1.fGamma + c0 * f0.fGamma ) * dt,
           f1.h     + ( c1 * f1.fh     + c0 * f0.fh     ) * dt,
           f1.l     + ( c1 * f1.fl     + c0 * f0.fl     ) * dt,
           f1.M     + ( c1 * f1.fM     + c0 * f0.fM     ) * dt);
}

// Iteration for three-step Adam's method
static inline DEVICE
void ThreeStepIteration(Layer &res, const Layer &f2, const Layer &f1, const Layer &f0,
                        const Unchangeable &params, real dt)
{
  constexpr auto c2 =  (real)23 / 12;
  constexpr auto c1 = -(real)16 / 12;
  constexpr auto c0 =  (real)5  / 12;
  
  SetLayer(res, params,
           f2.V     + ( c2 * f2.fV     + c1 * f1.fV     + c0 * f0.fV     ) * dt,
           f2.Gamma + ( c2 * f2.fGamma + c1 * f1.fGamma + c0 * f0.fGamma ) * dt,
           f2.h     + ( c2 * f2.fh     + c1 * f1.fh     + c0 * f0.fh     ) * dt,
           f2.l     + ( c2 * f2.fl     + c1 * f1.fl     + c0 * f0.fl     ) * dt,
           f2.M     + ( c2 * f2.fM     + c1 * f1.fM     + c0 * f0.fM     ) * dt);
}


// Helper to implement partial specialization for a function
// Enjoy the beaity of modern C++. No, COME BACK AND ENJOY!
template <unsigned int STEPS>
struct _IterationImpl
{
  template <typename LAYERS>
  static DEVICE void Perform(LAYERS &f, const Unchangeable &params, size_t nxt, real dt);
};
template <>
struct _IterationImpl<1>
{
  template <typename LAYERS>
  static DEVICE void Perform(LAYERS &f, const Unchangeable &params, size_t nxt, real dt)
  { OneStepIteration(f[nxt], f[(nxt + 1) & 1], params, dt); }
};
template <>
struct _IterationImpl<2>
{
  template <typename LAYERS>
  static DEVICE void Perform(LAYERS &f, const Unchangeable &params, size_t nxt, real dt)
  { TwoStepIteration(f[nxt], f[(nxt + 1) % 3], f[(nxt + 2) % 3], params, dt); }
};
template <>
struct _IterationImpl<3>
{
  template <typename LAYERS>
  static DEVICE void Perform(LAYERS &f, const Unchangeable &params, size_t nxt, real dt)
  { ThreeStepIteration(f[nxt], f[(nxt + 1) & 3], f[(nxt + 2) & 3], f[(nxt + 3) & 3], params, dt); }
};


// Performs uni-step iteration of Adams' method using a cycled buffer with layers
template <unsigned int STEPS, typename LAYERS = std::array<Layer, STEPS + 1>>
static inline DEVICE
void Iteration(LAYERS &f, const Unchangeable &params, size_t nxt, real dt)
{ _IterationImpl<STEPS>::Perform(f, params, nxt, dt); }


} // namespace Adams
