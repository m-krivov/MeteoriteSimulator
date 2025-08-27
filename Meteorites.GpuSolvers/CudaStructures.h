#pragma once
#include "Meteorites.GpuSolvers/CudaDefs.h"

#include "Meteorites.Core/Case.h"
#include "Meteorites.Core/Constants.h"


class Layer
{
  public:
  real V_, gamma_, h_, l_, M_;
  real fV_, fgamma_, fh_, fl_, fM_;

  __device__ void Set(real new_V, real new_gamma, real new_h, real new_l, real new_M, const Case& params)
  {
    V_ = new_V;
    gamma_ = new_gamma;
    h_ = new_h;
    l_ = new_l;
    M_ = new_M;

    if (new_M <= (real)0.0) // probably, 'dt' is too large
    {
      fV_ = fh_ = fl_ = fM_ = (real)0.0;
    } else {
      auto sin_gamma = sin(gamma_);
      auto cos_gamma = cos(gamma_);
      auto g = Constants::g(h_);
      auto rho_a = Constants::RhoAtm(h_);
      auto midsection = Constants::Midsection(M_, params.Rho);

      fV_ = -params.Cd * rho_a * V_ * V_ * midsection / (2 * M_) + g * sin_gamma;
      fgamma_ =
          +g * cos_gamma / V_ - V_ * cos_gamma / Constants::R() - params.Cl * rho_a * V_ * midsection / (2 * M_);
      fh_ = -V_ * sin_gamma;
      fl_ = V_ * (Constants::R() / (Constants::R() + h_)) * cos_gamma;
      fM_ = -(params.Ch * rho_a * V_ * V_ * V_ * midsection / 2) / params.H;
    }
  }
};

class Record
{
  public:
  real t_, m_, v_, h_, l_, gamma_;

  /*test*/ void print() { printf("rec: %f %f %f %f %f %f\n", t_, m_, v_, h_, l_, gamma_); }
};

template <unsigned int STEPS, unsigned int ITERS> class ThreadContext
{
  public:
  Layer steps[STEPS + 1];
  Case params;
  size_t nxt;
  real t;
  size_t timestamp;
  size_t curr_case_num;
  bool ended;
};
