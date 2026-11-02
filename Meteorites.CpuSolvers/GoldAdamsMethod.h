#pragma once
#include "Meteorites.Core/Defs.h"

#include "Meteorites.Core/Solvers/Adams.h"
#include "Meteorites.Core/Solvers/ISolver.h"

// An unified implementation for one-step, two-step and three-step Adams method
template <unsigned int STEPS>
void GoldAdamsMethod(const VirtualMeteoroid &problem, real dt, real timeout,
                     const IFunctional &functional, ISimulationRecorder &results)
{
  assert(1u <= STEPS && STEPS <= 3u);   // not adapted for other steps
  assert(problem.M0 > (real)0.0);
  assert(problem.V0 > (real)0.0);
  auto t_next = results.Started(problem);

  // Prepare the initial state and coefficients, verify them
  Adams::Unchangeable params(problem);
  Adams::State curr_layer = { problem.V0, problem.Gamma0, problem.h0, problem.l0, problem.M0 };
  Adams::State steps[STEPS];
  real t = (real)0.0;

  // Compute values for initial steps
  t_next = results.Store(t, curr_layer.M, curr_layer.V, curr_layer.h, curr_layer.l, curr_layer.Gamma);

  Adams::OneStepIteration(curr_layer, steps[0], curr_layer, params, dt);
  t += dt;
  t_next = results.Store(t, curr_layer.M, curr_layer.V, curr_layer.h, curr_layer.l, curr_layer.Gamma);
  
  if constexpr (STEPS >= 2) {
    Adams::TwoStepIteration(curr_layer, steps[1], curr_layer, steps[0], params, dt);
    t += dt;
    t_next = results.Store(t, curr_layer.M, curr_layer.V, curr_layer.h, curr_layer.l, curr_layer.Gamma);
  }

  if constexpr (STEPS >= 3) {
    Adams::ThreeStepIteration(curr_layer, steps[2], curr_layer, steps[1], steps[0], params, dt);
    t += dt;
    t_next = results.Store(t, curr_layer.M, curr_layer.V, curr_layer.h, curr_layer.l, curr_layer.Gamma);
  }

  // Prepare buffers for values that are wanted by functional
  const real *timestamps = nullptr;
  size_t n_timestamps = 0;
  functional.GetTimeStamps(n_timestamps, timestamps);
  assert(n_timestamps > 0);
  assert(timestamps != nullptr);
  
  size_t timestamp = 0;
  std::vector<real> V_arg(n_timestamps, (real)0.0f), h_arg(n_timestamps, (real)0.0f);

  // The main loop: perform simulation until meteorite is not burnt, collided or timeouted
  size_t nxt = 0;
  while (t < timeout)
  {
    // If necessery, update the functional's arguments
    if (timestamp < n_timestamps && t >= timestamps[timestamp])
    {
      V_arg[timestamp] = curr_layer.V;
      h_arg[timestamp] = curr_layer.h;
      timestamp += 1;
    }

    // Compute values for the next step, store them
    Adams::Iteration<STEPS>(curr_layer, steps, params, nxt, dt);

    t += dt;
    if (t >= t_next)
    { t_next = results.Store(t, curr_layer.M, curr_layer.V, curr_layer.h, curr_layer.l, curr_layer.Gamma); }
    nxt = (nxt + 1) % STEPS;
   
    // Check, should we stop the simulation?
    if (curr_layer.M <= (real)0.01)
    {
      results.Finished(ISimulationRecorder::Reason::Burnt,
                       functional.Compute(timestamp, &V_arg[0], &h_arg[0]));
      return;
    }
    if (curr_layer.h <= (real)0.0)
    {
      results.Finished(ISimulationRecorder::Reason::Collided,
                       functional.Compute(timestamp, &V_arg[0], &h_arg[0]));
      return;
    }
  }

  // Looks like something goes wrong
  results.Finished(ISimulationRecorder::Reason::Timeouted,
                   functional.Compute(timestamp, &V_arg[0], &h_arg[0]));
}

inline auto CreateGoldAdamsMethod(NumericalAlgorithm alg) -> auto
{
  switch (alg)
  {
    case NumericalAlgorithm::ONE_STEP_ADAMS:
      return GoldAdamsMethod<1>;

    case NumericalAlgorithm::TWO_STEP_ADAMS:
      return GoldAdamsMethod<2>;

    case NumericalAlgorithm::THREE_STEP_ADAMS:
      return GoldAdamsMethod<3>;

    default:
      throw std::invalid_argument("Unknown numerical algorithm");
  }
}
