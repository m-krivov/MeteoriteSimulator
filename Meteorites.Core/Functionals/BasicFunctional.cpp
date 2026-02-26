#include "BasicFunctional.h"

BasicFunctional::BasicFunctional(const std::shared_ptr<const IMeteorite> &meteorite,
                                 real lambda_v, real lambda_h)
  : lambda_h_(lambda_h), lambda_v_(lambda_v)
{
  assert(lambda_h >= 0.0);
  assert(lambda_v >= 0.0);

  size_t records = 0;
  const real *time = nullptr, *v = nullptr, *h = nullptr;
  meteorite->Trajectory(records, time, v, h);
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

BasicFunctional::BasicFunctional(const std::shared_ptr<const IMeteorite> &meteorite,
                                 real lambda_v, real lambda_h,
                                 const std::vector<real> &weights)
  : BasicFunctional(meteorite, lambda_v, lambda_h)
{
  assert(weights.size() == time_.size());
  weights_ = weights;
  
  // Verify all weights are non-negative
  for (const auto &w : weights_)
  { assert(w >= 0.0); }
}

void BasicFunctional::GetTimeStamps(size_t &num, const real *&values) const
{
  num    = time_.size();
  values = time_.data();
}

std::vector<real> BasicFunctional::GenerateDecayingWeights(size_t count, real decay_factor)
{
  assert(count > 0);
  assert(decay_factor >= 0.0 && decay_factor <= 1.0);
  
  std::vector<real> weights(count);
  
  if (decay_factor == 0.0)
  {
    // Uniform weights
    for (size_t i = 0; i < count; i++)
    { weights[i] = (real)1.0; }
  }
  else
  {
    // Exponential decay: weight = exp(-decay_factor * i / (count - 1))
    for (size_t i = 0; i < count; i++)
    {
      real t = (real)i / (real)(count - 1);
      weights[i] = std::exp(-decay_factor * (real)3.0 * t);  // factor of 3 for reasonable decay
    }
  }
  
  return weights;
}
