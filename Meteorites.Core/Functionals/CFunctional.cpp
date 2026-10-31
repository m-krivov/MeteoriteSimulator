#include "CFunctional.h"


double CFunctional::Compute(size_t num, const real *v, const real *h) const
{
  decltype(auto) time = Time();
  decltype(auto) v_ref = Velocity();
  decltype(auto) h_ref = Height();

  if (num < time.size())
  { return std::numeric_limits<double>::max(); }

  double max_v_error = 0.0, max_h_error = 0.0;
  size_t n = std::min(num, time.size());
  for (size_t i = 0; i < n; i++)
  {
    double dv = std::abs((v_ref[i] - v[i]) / v_ref[0]);
    double dh = std::abs((h_ref[i] - h[i]) / h_ref[0]);
    
    max_v_error = std::max(max_v_error, dv);
    max_h_error = std::max(max_h_error, dh);
  }

  return (max_v_error * LambdaV() + max_h_error * LambdaH());
}
