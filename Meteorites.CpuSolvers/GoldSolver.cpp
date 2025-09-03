#include "GoldSolver.h"

#include "Meteorites.Core/Adams.h"

namespace
{

// An unified implementation for one-step, two-step and three-step Adams method
template <unsigned int STEPS>
void AdamsMethod(const Case &problem, const IFunctional &functional, real dt, real timeout,
                 IResultFormatter &results)
{
  assert(1u <= STEPS && STEPS <= 3u);   // not adapted for other steps
  assert(problem.M0 > (real)0.0);
  assert(problem.V0 > (real)0.0);
  auto t_next = results.Started(problem);

  // Prepare the initial state and coefficients, verify them
  Adams::Unchangeable params(problem);
  std::array<Adams::Layer, STEPS + 1> steps;
  real t = (real)0.0;

  // Compute values for initial steps
  Adams::SetLayer(steps[STEPS], params,
                  problem.V0, problem.Gamma0,
                  problem.h0, problem.l0, problem.M0);
  t_next = results.Store(t, steps[STEPS].M, steps[STEPS].V,
                         steps[STEPS].h, steps[STEPS].l, steps[STEPS].Gamma);

  Adams::OneStepIteration(steps[STEPS - 1], steps[STEPS], params, dt);
  t += dt;
  t_next = results.Store(t, steps[STEPS - 1].M, steps[STEPS - 1].V,
                         steps[STEPS - 1].h, steps[STEPS - 1].l, steps[STEPS - 1].Gamma);
  
  if constexpr (STEPS >= 2) {
    Adams::TwoStepIteration(steps[STEPS - 2], steps[STEPS - 1], steps[STEPS], params, dt);
    t += dt;
    t_next = results.Store(t, steps[STEPS - 2].M, steps[STEPS - 2].V,
                           steps[STEPS - 2].h, steps[STEPS - 2].l, steps[STEPS - 2].Gamma);
  }

  if constexpr (STEPS >= 3) {
    Adams::ThreeStepIteration(steps[STEPS - 3], steps[STEPS - 2], steps[STEPS - 1], steps[STEPS], params, dt);
    t += dt;
    t_next = results.Store(t, steps[STEPS - 3].M, steps[STEPS - 3].V,
                           steps[STEPS - 3].h, steps[STEPS - 3].l, steps[STEPS - 3].Gamma);
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
  size_t nxt = STEPS;
  while (t < timeout)
  {
    // If necessery, update the functional's arguments
    if (timestamp < n_timestamps && t >= timestamps[timestamp])
    {
      const auto &step = steps[(nxt + 1) % (STEPS + 1)];
      V_arg[timestamp] = step.V;
      h_arg[timestamp] = step.h;
      timestamp += 1;
    }

    // Compute values for the next step, store them
    Adams::Iteration<STEPS>(steps, params, nxt, dt);
    auto M = steps[nxt].M;
    auto h = steps[nxt].h;

    t += dt;
    if (t >= t_next)
    { t_next = results.Store(t, M, steps[nxt].V, h, steps[nxt].l, steps[nxt].Gamma); }
    nxt = (nxt + STEPS) % (STEPS + 1);
   
    // Check, should we stop the simulation?
    if (M <= (real)0.01)
    {
      results.Finished(IResultFormatter::Reason::Burnt,
                       functional.Compute(timestamp, &V_arg[0], &h_arg[0]));
      return;
    }
    if (h <= (real)0.0)
    {
      results.Finished(IResultFormatter::Reason::Collided,
                       functional.Compute(timestamp, &V_arg[0], &h_arg[0]));
      return;
    }
  }

  // Looks like something goes wrong
  results.Finished(IResultFormatter::Reason::Timeouted,
                   functional.Compute(timestamp, &V_arg[0], &h_arg[0]));
}

} // unnamed namespace


void GoldSolver::Solve(const Case &problem, const IFunctional &functional, IResultFormatter &results)
{
  switch (Algorithm())
  {
    case NumericalAlgorithm::ONE_STEP_ADAMS:
      AdamsMethod<1>(problem, functional, Dt(), Timeout(), results);
      break;

    case NumericalAlgorithm::TWO_STEP_ADAMS:
      AdamsMethod<2>(problem, functional, Dt(), Timeout(), results);
      break;

    case NumericalAlgorithm::THREE_STEP_ADAMS:
      AdamsMethod<3>(problem, functional, Dt(), Timeout(), results);
      break;

    default:
      throw std::runtime_error("unknown numverical algorithm");
  }
}
