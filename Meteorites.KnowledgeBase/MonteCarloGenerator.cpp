#include "MonteCarloGenerator.h"


MonteCarloGenerator::MonteCarloGenerator(const IMeteorite &meteorite,
                                         const ParameterSet &range,
                                         size_t n_cases, uint32_t seed)
  : gen_(seed), dist_((real)0.0f, (real)1.0f),
    range_(range), cur_(0), n_cases_(n_cases),
    step_(0.0f), threshold_(std::numeric_limits<float>::max())
{
  size_t records = 0;
  const real *time, *v, *h;
  meteorite.Trajectory(records, time, v, h);
  assert(records >= 2);
  assert(v[0] > (real)0.0f);
  assert(time[1] - time[0] > (real)0.0f);

  v0_ = v[0];
  h0_ = h[0];
}

void MonteCarloGenerator::OnProgress(const std::function<void(float)> &callback, float step)
{
  assert(step >= 0.0f);
  callback_ = callback;
  step_ = step;
  threshold_ = step;
}

namespace
{

real SelectRandom(const std::pair<real, real> &min_max, real weight)
{
  assert(weight >= (real)0.0f);
  assert(weight <= (real)1.0f);

  const auto min_value = min_max.first;
  const auto max_value = min_max.second;
  assert(min_value <= max_value);
  return min_value + (max_value - min_value) * weight;
}

} // unnamed namespace

bool MonteCarloGenerator::Next(Case &problem)
{
  if (cur_ < n_cases_)
  {
    problem = Case(SelectRandom(range_.H(),   dist_(gen_)),
                   SelectRandom(range_.Ch(),  dist_(gen_)),
                   SelectRandom(range_.Rho(), dist_(gen_)),
                   SelectRandom(range_.Cd(),  dist_(gen_)),
                   SelectRandom(range_.Cl(),  dist_(gen_)),
                   SelectRandom(range_.M0(),  dist_(gen_)),
                   v0_, h0_,
                   SelectRandom(range_.Gamma0(), dist_(gen_)));
    cur_ += 1;

    auto ratio = (float)cur_ / n_cases_;
    if (ratio >= threshold_ || cur_ == n_cases_)
    {
      assert(callback_);
      callback_(ratio);
      threshold_ += step_;
    }

    return true;
  }
  else
  { return false; }
}
