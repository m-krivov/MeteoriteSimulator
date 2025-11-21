#include "BasicFunctional.h"

BasicFunctional::BasicFunctional(const IMeteorite &meteorite, real lambda_v, real lambda_h)
  : lambda_h_(lambda_h), lambda_v_(lambda_v)
{
  assert(lambda_h >= 0.0);
  assert(lambda_v >= 0.0);

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

  assert(v_[0] > 0.0);
  assert(h_[0] > 0.0);
}

void BasicFunctional::GetTimeStamps(size_t &num, const real *&values) const
{
  num    = time_.size();
  values = time_.data();
}
