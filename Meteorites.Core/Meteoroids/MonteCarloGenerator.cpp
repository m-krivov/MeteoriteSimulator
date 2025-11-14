#include "MonteCarloGenerator.h"


MonteCarloGenerator::MonteCarloGenerator(const IMeteorite &meteorite,
                                         const ParameterSet &range,
                                         size_t n_cases, uint64_t seed)
  : gen_(seed), dist_((real)0.0f, (real)1.0f),
    range_(range), n_cases_(n_cases)
{
  assert(n_cases >= 1);
  
  size_t records = 0;
  const real *time, *v, *h;
  meteorite.Trajectory(records, time, v, h);
  assert(records >= 2);
  assert(v[0] > (real)0.0f);
  assert(time[0] == 0.0f);
  assert(time[1] - time[0] > (real)0.0f);

  v0_ = v[0];
  h0_ = h[0];
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

bool MonteCarloGenerator::MoveNext()
{
  if (case_ < n_cases_)
  {
    current_ = VirtualMeteoroid(SelectRandom(range_.H(),   dist_(gen_)),
                                SelectRandom(range_.Ch(),  dist_(gen_)),
                                SelectRandom(range_.Rho(), dist_(gen_)),
                                SelectRandom(range_.Cd(),  dist_(gen_)),
                                SelectRandom(range_.Cl(),  dist_(gen_)),
                                SelectRandom(range_.M0(),  dist_(gen_)),
                                v0_, h0_,
                                SelectRandom(range_.Gamma0(), dist_(gen_)));
    case_ += 1;
    MovedNext(case_, n_cases_);
    return true;
  }
  else
  { return false; }
}

void MonteCarloGenerator::Reset()
{
  case_ = 0;
  gen_.seed(seed_);
  current_ = VirtualMeteoroid{};
}
