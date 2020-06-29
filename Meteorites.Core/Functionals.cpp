#include "Functionals.h"

L2Functional::L2Functional(const IMeteorite &meteorite)
{
  size_t records = 0;
  const real *time = nullptr, *v = nullptr, *h = nullptr;
  meteorite.Trajectory(records, time, v, h);
  assert(records != 0);
  assert(time != nullptr);
  assert(v != nullptr);
  assert(h != nullptr);

  time_.resize(records);
  v_.resize(records);
  h_.resize(records);
  for (size_t i = 0; i < records; i++)
  {
    time_[i] = time[i];
    v_[i] = v[i];
    h_[i] = h[i];
  }
}

void L2Functional::GetTimeStamps(size_t &num, const real *&values) const
{
  num    = time_.size();
  values = time_.data();
}

double L2Functional::Compute(size_t num, const real *v, const real *h) const
{
  if (num < time_.size())
  { return std::numeric_limits<double>::max(); }

  double v_sum = 0.0, h_sum = 0.0;
  size_t n = std::min(num, time_.size());
  for (size_t i = 0; i < n; i++)
  {
    auto dv = (v_[i] - v[i]) / v_[0];
    v_sum += dv * dv;

    auto dh = (h_[i] - h[i]) / h_[0];
    h_sum += dh * dh;
  }

  return (std::sqrt(v_sum) + std::sqrt(h_sum)) / std::sqrt(time_.size());
}
