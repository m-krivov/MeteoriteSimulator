#pragma once

#include "Meteorites.Core/Case.h"
#include "GPUParameters.h"

/*test*/ #include<iostream>

__device__ inline double cudaExp(double x) { return exp(x); }
__device__ inline float cudaExp(float x) { return expf(x); }
__device__ inline double cudaPow(double x, double y) { return pow(x, y); }
__device__ inline float cudaPow(float x, float y) { return powf(x, y); }
__device__ inline double cudaSin(double x) { return sin(x); }
__device__ inline float cudaSin(float x) { return sinf(x); }
__device__ inline double cudaCos(double x) { return cos(x); }
__device__ inline float cudaCos(float x) { return cosf(x); }

class CudaConstants
{
  public:
    __device__ static constexpr real R() { return (real)6371000; }

    __device__ static constexpr real g(real h)
    {
      return ((real)(6.67428 * 5.9726) / ((6371 + h / 1000) * (6371 + h / 1000)) * (real)1e7);
    }

    __device__ static constexpr real g() { return g(0); }

    __device__ static real rho_a(real h)
    {
      constexpr real rho_0 = (real)1.125;
      constexpr real alpha = (real)(0.029 * 9.8 / (8.31 * 273));
      return rho_0 * cudaExp(-alpha * h);
    }

    __device__ static real Midsection(real M, real rho)
    {
      return (real)(M_PI * pow(3.0 / (4.0 * M_PI), 2.0 / 3.0)) * pow(M / rho, (real)(2.0 / 3.0));
    }
};

class CudaCase
{
  public:
    real H_, ch_, rho_, cd_, cl_;
    real m0_, v0_, h0_, gamma0_;
};

class Layer
{
  public:
    real V_, gamma_, h_, l_, M_;
    real fV_, fgamma_, fh_, fl_, fM_;

    __device__ void Set(real new_V, real new_gamma, real new_h, real new_l, real new_M, const CudaCase& params)
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
        auto sin_gamma = cudaSin(gamma_);
        auto cos_gamma = cudaCos(gamma_);
        auto g = CudaConstants::g(h_);
        auto rho_a = CudaConstants::rho_a(h_);
        auto midsection = CudaConstants::Midsection(M_, params.rho_);

        fV_ = -params.cd_ * rho_a * V_ * V_ * midsection / (2 * M_) + g * sin_gamma;
        fgamma_ =
            +g * cos_gamma / V_ - V_ * cos_gamma / CudaConstants::R() - params.cl_ * rho_a * V_ * midsection / (2 * M_);
        fh_ = -V_ * sin_gamma;
        fl_ = V_ * (CudaConstants::R() / (CudaConstants::R() + h_)) * cos_gamma;
        fM_ = -(params.ch_ * rho_a * V_ * V_ * V_ * midsection / 2) / params.H_;
      }
    }
};

class Record
{
  public:
    real t_, m_, v_, h_, l_, gamma_;
};

template <unsigned int STEPS, unsigned int ITERS> class ThreadContext
{
  public:
    Layer steps[STEPS + 1];
    CudaCase params;
    size_t nxt;
    real t;
    size_t timestamp;
    size_t curr_case_num;
    bool ended;
};
