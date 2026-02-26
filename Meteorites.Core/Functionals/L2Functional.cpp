#include "L2Functional.h"


double L2Functional::Compute(size_t num, const real *v, const real *h) const
{
  decltype(auto) time = Time();
  decltype(auto) v_ref = Velocity();
  decltype(auto) h_ref = Height();

  if (num < time.size())
  { return std::numeric_limits<double>::max(); }

  double v_sum = 0.0, h_sum = 0.0;
  size_t n = std::min(num, time.size());
  
  if (HasWeights())
  {
    decltype(auto) weights = Weights();
    for (size_t i = 0; i < n; i++)
    {
      auto dv = (v_ref[i] - v[i]) / v_ref[0];
      v_sum += weights[i] * dv * dv;

      auto dh = (h_ref[i] - h[i]) / h_ref[0];
      h_sum += weights[i] * dh * dh;
    }
  }
  else
  {
    for (size_t i = 0; i < n; i++)
    {
      auto dv = (v_ref[i] - v[i]) / v_ref[0];
      v_sum += dv * dv;

      auto dh = (h_ref[i] - h[i]) / h_ref[0];
      h_sum += dh * dh;
    }
  }

  return (std::sqrt(v_sum) * LambdaV() + std::sqrt(h_sum) * LambdaH()) / std::sqrt(time.size());
}
