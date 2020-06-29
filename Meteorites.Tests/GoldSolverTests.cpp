#include "TestDefs.h"

#include "Meteorites.CpuSolvers/GoldSolver.h"


TEST_CLASS(GoldSolverTests)
{
  public:
	  TEST_METHOD(Vacuum_VerticalSpeed)
	  {
      GoldSolver solver(GoldSolver::ONE_STEP_ADAMS);
      solver.Configure(dt_sim, 3600.0f);
      Case problem(1.0f, 0.0f, 2000.0f, 0.0f, 0.0f,
                   1.0f, 0.5f, 1000.0f, (real)std::_Pi / 2);
      FakeFunctional f;
      BufferingFormatter fmt(0.1f);

      solver.Solve(problem, f, fmt);
      decltype(auto) log = ExtractLogPtrs<1>(fmt)[0];
      decltype(auto) last_record = log->records[log->records.size() - 1];
      
      Assert::IsTrue(std::abs(last_record.M - problem.M0()) < (real)1e-3);
      auto a = -Constants::g() / 2.0;
      auto b = -problem.V0();
      auto c = problem.h0() - last_record.h;
      auto d = b * b - 4 * a * c;
      auto t = (-b - std::sqrt(d)) / (2 * a);
      Assert::IsTrue(std::abs(last_record.t - t) < 1e-2f);
      Assert::IsTrue(std::abs(last_record.V - (problem.V0() + Constants::g() * t)) < 1e-1f);
	  }

    TEST_METHOD(Vacuum_HorizontalSpeed)
    {
      GoldSolver solver(GoldSolver::ONE_STEP_ADAMS);
      solver.Configure(dt_sim, 3600.0f);
      Case problem(1.0f, 0.0f, 2000.0f, 0.0f, 0.0f,
                   1.0f, 2.5f, 500.0f, 0.0f);
      FakeFunctional f;
      BufferingFormatter fmt(0.1f);

      solver.Solve(problem, f, fmt);
      decltype(auto) log = ExtractLogPtrs<1>(fmt)[0];
      decltype(auto) last_record = log->records[log->records.size() - 1];

      Assert::IsTrue(std::abs(last_record.M - problem.M0()) < (real)1e-3);
      auto a = -Constants::g() / 2.0;
      auto b = 0.0f;
      auto c = problem.h0() - last_record.h;
      auto d = b * b - 4 * a * c;
      auto t = (-b - std::sqrt(d)) / (2 * a);
      auto dist = problem.V0() * std::cos(problem.Gamma0()) * t;
      Assert::IsTrue(std::abs(last_record.t - t) < 1e-2f);
      Assert::IsTrue(std::abs(last_record.V - Constants::g() * t) < 1e-1f);
      Assert::IsTrue(std::abs(last_record.l - dist) < 2e-1f);
    }

    TEST_METHOD(Atmosphere_BrakingForce)
    {
      GoldSolver solver(GoldSolver::ONE_STEP_ADAMS);
      solver.Configure(dt_sim, 3600.0f);
      FakeFunctional f;
      BufferingFormatter fmt(0.0f);

      for (auto Cd : { 0.0f, 0.5f, 1.0f, 1.5f, 2.0f })
      {
        Case problem(1.0f, 0.0f, 2000.0f, Cd, 0.0f,
                     1.0f, 50.0f, 500.0f, -(real)std::_Pi / 4);
        solver.Solve(problem, f, fmt);
      }
      
      decltype(auto) logs = ExtractLogPtrs<5>(fmt);
      std::array<std::pair<real, real>, 5> peaks;
      for (size_t i = 0; i < logs.size(); i++)
      {
        decltype(auto) records = logs[i]->records;
        auto peak = std::make_pair(records[0].h, records[0].t);
        bool peak_reached = false;
        for (size_t j = 1; j < records.size(); j++)
        {
          if (!peak_reached)
          {
            if (records[j].h < peak.first)
            { peak_reached = true; }
            else
            { peak = std::make_pair(records[j].h, records[j].t); }
          }
          else
          { Assert::IsTrue(records[j].h < peak.first); }
        } // for j

        peaks[i] = peak;
      } // for i

      for (size_t i = 1; i < peaks.size(); i++)
      {
        Assert::IsTrue(peaks[i - 1].first > peaks[i].first * 1.005f);
        Assert::IsTrue(peaks[i - 1].second > peaks[i].second * 1.005f);
      }
    }

    TEST_METHOD(Atmosphere_LiftingForce)
    {
      GoldSolver solver(GoldSolver::ONE_STEP_ADAMS);
      solver.Configure(dt_sim, 3600.0f);
      FakeFunctional f;
      BufferingFormatter fmt(0.0f);

      for (auto Cl : { 0.0f, 0.05f, 0.1f, 0.15f })
      {
        Case problem(1.0f, 0.0f, 2000.0f, 0.0f, Cl,
                     1.0f, 100.0f, 500.0f, 0.0f);
        solver.Solve(problem, f, fmt);
      }

      decltype(auto) logs = ExtractLogPtrs<4>(fmt);
      for (size_t i = 0; i < logs.size() - 1; i++)
      {
        Assert::IsTrue(logs[i]->reason     == IResultFormatter::Reason::Collided);
        Assert::IsTrue(logs[i + 1]->reason == IResultFormatter::Reason::Collided);


        Assert::IsTrue(logs[i]->records.size() > 1000);
        for (size_t j = 1000; j < logs[i]->records.size(); j++)
        {
          if (j >= logs[i + 1]->records.size()) { break; }
          Assert::IsTrue(logs[i + 1]->records[j].h > logs[i]->records[j].h * 1.001f);
        } // for j
      } // for i
    }

    TEST_METHOD(Atmosphere_HeatExchange)
    {
      GoldSolver solver(GoldSolver::ONE_STEP_ADAMS);
      solver.Configure(dt_sim, 3600.0f);
      FakeFunctional f;
      BufferingFormatter fmt(0.0f);

      for (auto coeffs : { std::make_pair(1e6f, 0.3f),
                           std::make_pair(2e6f, 0.2f),
                           std::make_pair(3e6f, 0.1f) })
      {
        Case problem(coeffs.first, coeffs.second, 2000.0f, 0.0f, 0.0f,
                     10.0f, 100.0f, 2000.0f, 0.0f);
        solver.Solve(problem, f, fmt);
      }

      decltype(auto) logs = ExtractLogPtrs<3>(fmt);
      for (size_t i = 0; i < logs.size(); i++)
      {
        decltype(auto) records = logs[i]->records;
        for (size_t j = 1; j < records.size(); j++)
        { Assert::IsTrue(records[j - 1].M > records[j].M); }
      } // for i

      auto n = std::min(logs[0]->records.size(), logs[1]->records.size());
      n = std::min(n, logs[2]->records.size());
      Assert::IsTrue(n > 1000);
      for (size_t i = 1000; i < n; i++)
      {
        Assert::IsTrue(logs[1]->records[i].M > logs[0]->records[i].M * 1.0001f);
        Assert::IsTrue(logs[2]->records[i].M > logs[1]->records[i].M * 1.0001f);
      }
    }

    TEST_METHOD(Method_TwoThreeSteps)
    {
      GoldSolver one_step(GoldSolver::ONE_STEP_ADAMS),
                 two_step(GoldSolver::TWO_STEP_ADAMS),
                 three_step(GoldSolver::THREE_STEP_ADAMS);
      for (auto solver : { &one_step, &two_step, &three_step })
      { solver->Configure(dt_sim, 3600.0f); }
      FakeFunctional f;

      std::vector<Case> problems;
      problems.emplace_back(Case(1e6f, 0.35f, 3500.0f, 1.1f, 0.1f,
                                 10.0f, 14e3f, 60e3f, 0.0f));
      problems.emplace_back(Case(0.5e6f, 0.5f, 2000.0f, 1.5f, 0.05f,
                                 5.0f, 10e3f, 55e3f, (real)std::_Pi / 9));
      for (auto &problem : problems)
      {
        BufferingFormatter fmt(0.01f);
        for (auto solver : { &one_step, &two_step, &three_step })
        { solver->Solve(problem, f, fmt); }
        auto logs = ExtractLogPtrs<3>(fmt);
        auto n = std::min(logs[0]->records.size(),
                          logs[1]->records.size());
        n = std::min(n, logs[2]->records.size());

        for (size_t i = 0; i < n; i++)
        {
          AreAlmostEqual(logs[0]->records[i].gamma,
                         logs[1]->records[i].gamma,
                         logs[2]->records[i].gamma, 1e-2f);
          AreAlmostEqual(logs[0]->records[i].h,
                         logs[1]->records[i].h,
                         logs[2]->records[i].h, 1e1f);
          AreAlmostEqual(logs[0]->records[i].M,
                         logs[1]->records[i].M,
                         logs[2]->records[i].M, 1e-1f);
          AreAlmostEqual(logs[0]->records[i].V,
                         logs[1]->records[i].V,
                         logs[2]->records[i].V, 1e1f);
        } // for j
      } // for problem
    }

    TEST_METHOD(Method_Precision)
    {
      GoldSolver solver(GoldSolver::ONE_STEP_ADAMS);
      Case problem(0.5e6f, 0.1f, 4000.0f, 1.0f, 0.01f,
                   25.0f, 5e3f, 30e3f, (real)std::_Pi / 4);
      FakeFunctional f;
      BufferingFormatter fmt(0.01f);
      for (auto dt : { 1e-3f, 1e-4f, 1e-5f })
      {
        solver.Configure(dt, 3600.0f);
        solver.Solve(problem, f, fmt);
      }

      auto logs = ExtractLogPtrs<3>(fmt);
      auto n = std::min(logs[0]->records.size(),
                        logs[1]->records.size());
      n = std::min(n, logs[2]->records.size());

      for (size_t i = 0; i < n; i++)
      {
        AreAlmostEqual(logs[0]->records[i].gamma,
                       logs[1]->records[i].gamma,
                       logs[2]->records[i].gamma, 1e-2f);
        AreAlmostEqual(logs[0]->records[i].h,
                       logs[1]->records[i].h,
                       logs[2]->records[i].h, 1e1f);
        AreAlmostEqual(logs[0]->records[i].M,
                       logs[1]->records[i].M,
                       logs[2]->records[i].M, 1e-1f);
        AreAlmostEqual(logs[0]->records[i].V,
                       logs[1]->records[i].V,
                       logs[2]->records[i].V, 1e1f);
      } // for j
    }

  private:
    template <size_t N>
    static decltype(auto) ExtractLogPtrs(const BufferingFormatter &fmt)
    {
      decltype(auto) logs = fmt.Logs();
      Assert::AreEqual(logs.size(), N);

      std::array<const BufferingFormatter::Log *, N> res;
      for (size_t i = 0; i < logs.size(); i++)
      {
        decltype(auto) log = logs[i];
        Assert::IsFalse(log.records.empty());
        Assert::IsTrue(log.reason != IResultFormatter::Reason::NA);
        res[i] = &log;
      }
      return res;
    }

    void AreAlmostEqual(real val1, real val2, real val3, real eps)
    {
      Assert::IsTrue(std::abs(val1 - val2) <= eps);
      Assert::IsTrue(std::abs(val1 - val3) <= eps);
      Assert::IsTrue(std::abs(val2 - val3) <= eps);
    }

    real dt_sim  = (real)1e-3;
};
