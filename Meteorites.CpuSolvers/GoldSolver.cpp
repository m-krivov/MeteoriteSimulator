#include "GoldSolver.h"

#include "Meteorites.Core/Constants.h"

#include <iostream>
#include <chrono>
#include <thread>

namespace
{

// Caches constants and problem's parameters that may be considered as unchangeable values
// This class is a thin wrapper for 'Case' and 'Constants'
class Unchangeable
{
  public:
    Unchangeable(const Case &problem)
      : H_{ problem.H() },
        Ch_{ problem.Ch() },
        Cd_{ problem.Cd() },
        Cl_{ problem.Cl() },
        rho_{ problem.Rho() },
        R_{ Constants::R() }
    { }
    Unchangeable(const Unchangeable &) = delete;
    Unchangeable &operator =(const Unchangeable &) = delete;

    real H() const { return H_; }
    real Ch() const { return Ch_; }
    real Cd() const { return Cd_; }
    real Cl() const { return Cl_; }
    real Rho() const { return rho_; }
    real R() const { return R_; }

  private:
    const real H_, Ch_, Cd_, Cl_, rho_;
    const real R_;
};

// Represents four primary parameters (velocity, angle, height and mass) at some point in time
// Also precomputes right parts of the ODE
class Layer
{
  public:
    Layer()
      : params_(nullptr),
        V_(0.0f), gamma_(0.0f), h_(0.0f), l_(0.0f), M_(0.0f),
        fV_(0.0f), fgamma_(0.0f), fh_(0.0f), fl_(0.0f), fM_(0.0f)
    { }
    Layer(const Layer &) = delete;
    Layer &operator =(const Layer &) = delete;

    void Init(const Unchangeable *params)
    {
      assert(params->Rho() > 0.0f);
      assert(params->H() > 1e-3f);
      params_ = params;
    }

    void Set(real new_V, real new_gamma, real new_h, real new_l, real new_M)
    {
      assert(params_ != nullptr);
      assert(new_V > 0.0f);

      V_ = new_V;
      gamma_ = new_gamma;
      h_ = new_h;
      l_ = new_l;
      M_ = new_M;

      if (new_M <= (real)0.0)   // probably, 'dt' is too large
      {
        fV_ = fh_ = fl_ = fM_ = (real)0.0;
      }
      else
      {
        auto sin_gamma = std::sin(gamma_);
        auto cos_gamma = std::cos(gamma_);
        auto g = Constants::g(h_);
        auto rho_a = Constants::rho_a(h_);
        auto midsection = Constants::Midsection(M_, params_->Rho());

        fV_ = -params_->Cd() * rho_a * V_ * V_ * midsection / (2 * M_)
              + g * sin_gamma;
        fgamma_ =  + g * cos_gamma / V_
                   - V_ * cos_gamma / params_->R()
                   - params_->Cl() * rho_a * V_ * midsection / (2 * M_);
        fh_ = -V_ * sin_gamma;
        fl_ = V_ * (params_->R() / (params_->R() + h_)) * cos_gamma;
        fM_ = - (params_->Ch() * rho_a * V_ * V_ * V_ * midsection / 2) / params_->H();
      }
    }

    real V() const { return V_; }
    real Gamma() const { return gamma_; }
    real h() const { return h_; }
    real l() const { return l_; }
    real M() const { return M_; }
    
    real fV() const { return fV_; }
    real fGamma() const { return fgamma_; }
    real fh() const { return fh_; }
    real fl() const { return fl_; }
    real fM() const { return fM_; }

  private:
    const Unchangeable *params_;
    real V_, gamma_, h_, l_, M_;
    real fV_, fgamma_, fh_, fl_, fM_;
};

template <unsigned int STEPS>
void UniStepAdams(std::array<Layer, STEPS + 1> &f, size_t nxt, real dt);

void OneStepAdams(Layer &res, const Layer &f0, real dt)
{
  res.Set(f0.V()     + f0.fV()     * dt,
          f0.Gamma() + f0.fGamma() * dt,
          f0.h()     + f0.fh()     * dt,
          f0.l()     + f0.fl()     * dt,
          f0.M()     + f0.fM()     * dt);
}

template <>
void UniStepAdams<1>(std::array<Layer, 2> &f, size_t nxt, real dt) {
  OneStepAdams(f[nxt], f[(nxt + 1) & 1], dt);
}

void TwoStepAdams(Layer &res, const Layer &f1, const Layer &f0, real dt)
{
  const auto c1 =  (real)1.5;
  const auto c0 = -(real)0.5;
  
  res.Set(f1.V()     + ( c1 * f1.fV()     + c0 * f0.fV()     ) * dt,
          f1.Gamma() + ( c1 * f1.fGamma() + c0 * f0.fGamma() ) * dt,
          f1.h()     + ( c1 * f1.fh()     + c0 * f0.fh()     ) * dt,
          f1.l()     + ( c1 * f1.fl()     + c0 * f0.fl()     ) * dt,
          f1.M()     + ( c1 * f1.fM()     + c0 * f0.fM()     ) * dt);
}

template <>
void UniStepAdams<2>(std::array<Layer, 3> &f, size_t nxt, real dt) {
  TwoStepAdams(f[nxt], f[(nxt + 1) % 3], f[(nxt + 2) % 3], dt);
}

void ThreeStepAdams(Layer &res, const Layer &f2, const Layer &f1, const Layer &f0, real dt)
{
  const auto c2 =  (real)23 / 12;
  const auto c1 = -(real)16 / 12;
  const auto c0 =  (real)5  / 12;
  
  res.Set(f2.V()     + ( c2 * f2.fV()     + c1 * f1.fV()     + c0 * f0.fV()     ) * dt,
          f2.Gamma() + ( c2 * f2.fGamma() + c1 * f1.fGamma() + c0 * f0.fGamma() ) * dt,
          f2.h()     + ( c2 * f2.fh()     + c1 * f1.fh()     + c0 * f0.fh()     ) * dt,
          f2.l()     + ( c2 * f2.fl()     + c1 * f1.fl()     + c0 * f0.fl()     ) * dt,
          f2.M()     + ( c2 * f2.fM()     + c1 * f1.fM()     + c0 * f0.fM()     ) * dt);
}

template <>
void UniStepAdams<3>(std::array<Layer, 4> &f, size_t nxt, real dt) {
  ThreeStepAdams(f[nxt], f[(nxt + 1) & 3], f[(nxt + 2) & 3], f[(nxt + 3) & 3], dt);
}

// An unified implementation for one-step, two-step and three-step Adams method
template <unsigned int STEPS>
void AdamsMethod(const Case &problem, const IFunctional &functional, real dt, real timeout,
                 IResultFormatter &results)
{
  assert(1u <= STEPS && STEPS <= 3u);   // not adapted for other steps
  assert(problem.M0() > (real)0.0);
  assert(problem.V0() > (real)0.0);
  auto t_next = results.Started(problem);

  // Prepare the initial state and coefficients, verify them
  Unchangeable params(problem);
  std::array<Layer, STEPS + 1> steps;
  for (size_t i = 0; i < steps.size(); i++) {
    steps[i].Init(&params);
  }
  real t = (real)0.0;

  // Compute values for initial steps
  steps[STEPS].Set(problem.V0(), problem.Gamma0(),
                   problem.h0(), problem.l0(), problem.M0());
  t_next = results.Store(t, steps[STEPS].M(), steps[STEPS].V(),
                         steps[STEPS].h(), steps[STEPS].l(), steps[STEPS].Gamma());
  /*test_output*//*printf("init_1_store: %f %f %f %f %f %f\n", t, steps[STEPS].M(), steps[STEPS].V(),
                         steps[STEPS].h(), steps[STEPS].l(), steps[STEPS].Gamma());
  std::cout << "t_next: " << t_next << "\n";*/

  OneStepAdams(steps[STEPS - 1], steps[STEPS], dt);
  t += dt;
  t_next = results.Store(t, steps[STEPS - 1].M(), steps[STEPS - 1].V(),
                         steps[STEPS - 1].h(), steps[STEPS - 1].l(), steps[STEPS - 1].Gamma());
  /*test_output*//*printf("init_2_store: %f %f %f %f %f %f\n", t, steps[STEPS - 1].M(), steps[STEPS - 1].V(),
                         steps[STEPS - 1].h(), steps[STEPS - 1].l(), steps[STEPS - 1].Gamma());
  std::cout << "t_next: " << t_next << "\n";*/
  
  if (STEPS >= 2) {
    TwoStepAdams(steps[STEPS - 2], steps[STEPS - 1], steps[STEPS], dt);
    t += dt;
    t_next = results.Store(t, steps[STEPS - 2].M(), steps[STEPS - 2].V(),
                           steps[STEPS - 2].h(), steps[STEPS - 2].l(), steps[STEPS - 2].Gamma());
    /*test_output*//*printf("init_3_store: %f %f %f %f %f %f\n", t, steps[STEPS - 2].M(), steps[STEPS - 2].V(),
                           steps[STEPS - 2].h(), steps[STEPS - 2].l(), steps[STEPS - 2].Gamma());
    std::cout << "t_next: " << t_next << "\n";*/
  }

  if (STEPS >= 3) {
    ThreeStepAdams(steps[STEPS - 3], steps[STEPS - 2], steps[STEPS - 1], steps[STEPS], dt);
    t += dt;
    t_next = results.Store(t, steps[STEPS - 3].M(), steps[STEPS - 3].V(),
                           steps[STEPS - 3].h(), steps[STEPS - 3].l(), steps[STEPS - 3].Gamma());
    /*test_output*//*printf("init_4_store: %f %f %f %f %f %f\n", t, steps[STEPS - 3].M(), steps[STEPS - 3].V(),
                           steps[STEPS - 3].h(), steps[STEPS - 3].l(), steps[STEPS - 3].Gamma());
    std::cout << "t_next: " << t_next << "\n";*/
  }

  // Prepare buffers for values that are wanted by functional
  const real *timestamps = nullptr;
  size_t n_timestamps = 0;
  functional.GetTimeStamps(n_timestamps, timestamps);
  assert(n_timestamps > 0);
  assert(timestamps != nullptr);
  
  size_t timestamp = 0;
  std::vector<real> V_arg(n_timestamps, (real)0.0f), h_arg(n_timestamps, (real)0.0f);

  // The main loop: perform simulation until meteorite is not burnt, collided or timeouted
  size_t nxt = STEPS;
  while (t < timeout)
  {
    // If necessery, update the functional's arguments
    if (timestamp < n_timestamps && t >= timestamps[timestamp])
    {
      const auto &step = steps[(nxt + 1) % (STEPS + 1)];
      V_arg[timestamp] = step.V();
      h_arg[timestamp] = step.h();
      timestamp += 1;
    }

    // Compute values for the next step, store them
    UniStepAdams<STEPS>(steps, nxt, dt);
    auto M = steps[nxt].M();
    auto h = steps[nxt].h();

    t += dt;

    //std::this_thread::sleep_for(std::chrono::milliseconds(100));
    //if (t >= 0.03) return;

    if (t >= t_next)
    { t_next = results.Store(t, M, steps[nxt].V(), h, steps[nxt].l(), steps[nxt].Gamma());
      /*test_output*//*printf("store: %f %f %f %f %f %f\n", t, M, steps[nxt].V(), h,
                             steps[nxt].l(), steps[nxt].Gamma());*/}
    nxt = (nxt + STEPS) % (STEPS + 1);
   
    // Check, should we stop the simulation?
    if (M <= (real)0.01)
    {
      results.Finished(IResultFormatter::Reason::Burnt,
                       functional.Compute(timestamp, &V_arg[0], &h_arg[0]));
      return;
    }
    if (h <= (real)0.0)
    {
      results.Finished(IResultFormatter::Reason::Collided,
                       functional.Compute(timestamp, &V_arg[0], &h_arg[0]));
      return;
    }
  }

  // Looks like something goes wrong
  results.Finished(IResultFormatter::Reason::Timeouted,
                   functional.Compute(timestamp, &V_arg[0], &h_arg[0]));
}

} // unnamed namespace


void GoldSolver::Solve(const Case &problem, const IFunctional &functional, IResultFormatter &results)
{
  switch (AdamsSteps())
  {
    case 1:
      AdamsMethod<1>(problem, functional, Dt(), Timeout(), results);
      break;

    case 2:
      AdamsMethod<2>(problem, functional, Dt(), Timeout(), results);
      break;

    case 3:
      AdamsMethod<3>(problem, functional, Dt(), Timeout(), results);
      break;

    default:
      assert(false);
  }
}

void GoldSolver::Solve(ICaseGenerator &generator,
                       const IFunctional &functional,
                       IResultFormatter &results)
{
  Case problem;
  //generator.Next(problem);
  while (generator.Next(problem))
  {
    /*test_output*//*std::cout << "gold case:" << problem.Cd() << " "
                  << problem.Ch() << " "
                  << problem.Cl() << " "
                  << problem.Gamma0() << " "
                  << problem.h0() << " "
                  << problem.H() << " "
                  << problem.l0() << " "
                  << problem.M0() << " "
                  << problem.Rho() << " "
                  << problem.V0() << std::endl;*/
    Solve(problem, functional, results);
  }
}
