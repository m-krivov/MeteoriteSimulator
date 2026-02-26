#include "L1Functional.h"


double L1Functional::Compute(size_t num, const real *v, const real *h) const
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
      v_sum += weights[i] * std::abs((v_ref[i] - v[i]) / v_ref[0]);
      h_sum += weights[i] * std::abs((h_ref[i] - h[i]) / h_ref[0]);
    }
  }
  else
  {
    for (size_t i = 0; i < n; i++)
    {
      v_sum += std::abs((v_ref[i] - v[i]) / v_ref[0]);
      h_sum += std::abs((h_ref[i] - h[i]) / h_ref[0]);
    }
  }

  return (v_sum * LambdaV() + h_sum * LambdaH()) / time.size();
}
