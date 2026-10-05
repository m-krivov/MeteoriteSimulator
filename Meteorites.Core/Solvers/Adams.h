#pragma once
#include "Meteorites.Core/Defs.h"

#include "Meteorites.Core/Meteoroids/VirtualMeteoroid.h"
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

  // Precomputed per-meteoroid constants:
  // inv_H             = 1 / H
  // inv_R             = 1 / R
  // midsection_factor = coeff * Rho^(-2/3), so that
  //     Midsection(M, Rho) == midsection_factor * M^(2/3)
  real inv_H             = (real)0.0;
  real inv_R             = (real)0.0;
  real midsection_factor = (real)0.0;

  Unchangeable() = default;
  DEVICE Unchangeable(const VirtualMeteoroid &problem)
    : H{ problem.H },
      Ch{ problem.Ch },
      Cd{ problem.Cd },
      Cl{ problem.Cl },
      Rho{ problem.Rho },
      R{ Constants::R() }
  {
    inv_H             = (real)1.0 / H;
    inv_R             = (real)1.0 / R;
    midsection_factor = Constants::Midsection((real)1.0, Rho);
  }
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

  Layer() = default;
  Layer(const Layer &) = default;
  Layer &operator =(const Layer &) = default;
};


// Computes the right part of Stulov's ODE (derivatives of Layer fields)
static inline DEVICE
void ComputeDLayer(Layer &dlayer, const Layer &layer, const Unchangeable &params)
{
  assert(params.Rho > 0.0f);
  assert(params.H > 1e-3f);
  assert(layer.V > 0.0f);

  // Branchless guard for non-positive mass ('dt is too large' case).
  // In the steady-state hot loop layer.M is always positive (the caller
  // checks M <= 0.01 after every step), so this path is almost always
  // inactive. A data-dependent select is cheaper than a branch: it removes
  // the branch-resolution stall and warp divergence at this point.
  const bool alive = layer.M > (real)0.0;
  const real M     = fmax(layer.M, (real)1e-30);

  auto sin_gamma = std::sin(layer.Gamma);
  auto cos_gamma = std::cos(layer.Gamma);
  auto g = Constants::g(layer.h);
  auto rho_a = Constants::RhoAtm(layer.h);
  auto midsection = params.midsection_factor * std::pow(M, (real)(2.0 / 3.0));
  auto inv_2M = (real)0.5 / M;

  const real dl_V = - params.Cd * rho_a * layer.V * layer.V * midsection * inv_2M
                    + g * sin_gamma;
  const real dl_G = + g * cos_gamma / layer.V
                    - layer.V * cos_gamma * params.inv_R
                    - params.Cl * rho_a * layer.V * midsection * inv_2M;
  const real dl_h = - layer.V * sin_gamma;
  const real dl_l = layer.V * (params.R / (params.R + layer.h)) * cos_gamma;
  const real dl_M = - (params.Ch * rho_a * layer.V * layer.V * layer.V
                       * midsection * (real)0.5) * params.inv_H;

  dlayer.V     = alive ? dl_V : (real)0.0;
  dlayer.Gamma = alive ? dl_G : (real)0.0;
  dlayer.h     = alive ? dl_h : (real)0.0;
  dlayer.l     = alive ? dl_l : (real)0.0;
  dlayer.M     = alive ? dl_M : (real)0.0;
}


// Iteration for one-step Adam's method
static inline DEVICE
void OneStepIteration(Layer &layer0, Layer &dlayer0,
                      const Unchangeable &params, real dt)
{
  ComputeDLayer(dlayer0, layer0, params);

  layer0 = { layer0.V     + dlayer0.V     * dt,
             layer0.Gamma + dlayer0.Gamma * dt,
             layer0.h     + dlayer0.h     * dt,
             layer0.l     + dlayer0.l     * dt,
             layer0.M     + dlayer0.M     * dt };
}

// Iteration for two-step Adam's method
static inline DEVICE
void TwoStepIteration(Layer &layer1, Layer &dlayer1, const Layer &dlayer0,
                      const Unchangeable &params, real dt)
{
  constexpr auto c1 =  (real)1.5;
  constexpr auto c0 = -(real)0.5;

  ComputeDLayer(dlayer1, layer1, params);

  layer1 = { layer1.V     + ( c1 * dlayer1.V     + c0 * dlayer0.V     ) * dt,
             layer1.Gamma + ( c1 * dlayer1.Gamma + c0 * dlayer0.Gamma ) * dt,
             layer1.h     + ( c1 * dlayer1.h     + c0 * dlayer0.h     ) * dt,
             layer1.l     + ( c1 * dlayer1.l     + c0 * dlayer0.l     ) * dt,
             layer1.M     + ( c1 * dlayer1.M     + c0 * dlayer0.M     ) * dt };
}

// Iteration for three-step Adam's method
static inline DEVICE
void ThreeStepIteration(Layer &layer2, Layer &dlayer2, const Layer &dlayer1, const Layer &dlayer0,
                        const Unchangeable &params, real dt)
{
  constexpr auto c2 =  (real)23 / 12;
  constexpr auto c1 = -(real)16 / 12;
  constexpr auto c0 =  (real)5  / 12;

  ComputeDLayer(dlayer2, layer2, params);
  
  layer2 = { layer2.V     + ( c2 * dlayer2.V     + c1 * dlayer1.V     + c0 * dlayer0.V     ) * dt,
             layer2.Gamma + ( c2 * dlayer2.Gamma + c1 * dlayer1.Gamma + c0 * dlayer0.Gamma ) * dt,
             layer2.h     + ( c2 * dlayer2.h     + c1 * dlayer1.h     + c0 * dlayer0.h     ) * dt,
             layer2.l     + ( c2 * dlayer2.l     + c1 * dlayer1.l     + c0 * dlayer0.l     ) * dt,
             layer2.M     + ( c2 * dlayer2.M     + c1 * dlayer1.M     + c0 * dlayer0.M     ) * dt };
}


// Helper to implement partial specialization for a function
// Enjoy the beaity of modern C++. No, COME BACK AND ENJOY!
template <unsigned int STEPS>
struct _IterationImpl
{
  template <typename LAYER>
  static DEVICE void Perform(LAYER &l, LAYER (&f)[STEPS], const Unchangeable &params, size_t nxt, real dt);
};
template <>
struct _IterationImpl<1>
{
  template <typename LAYER>
  static DEVICE void Perform(LAYER &l, LAYER f[], const Unchangeable &params, size_t nxt, real dt)
  { OneStepIteration(l, f[nxt], params, dt); }
};
template <>
struct _IterationImpl<2>
{
  template <typename LAYER>
  static DEVICE void Perform(LAYER &l, LAYER f[], const Unchangeable &params, size_t nxt, real dt)
  { TwoStepIteration(l, f[nxt], f[(nxt + 1) % 2], params, dt); }
};
template <>
struct _IterationImpl<3>
{
  template <typename LAYER>
  static DEVICE void Perform(LAYER &l, LAYER f[], const Unchangeable &params, size_t nxt, real dt)
  { ThreeStepIteration(l, f[nxt], f[(nxt + 2) % 3], f[(nxt + 1) % 3], params, dt); }
};


// Performs uni-step iteration of Adams' method using a cycled buffer with layers
template <unsigned int STEPS, typename LAYER>
static inline DEVICE
void Iteration(LAYER &l, LAYER f[STEPS], const Unchangeable &params, size_t nxt, real dt)
{ _IterationImpl<STEPS>::Perform(l, f, params, nxt, dt); }


} // namespace Adams
