#include "PrecisionEstimator.h"

#include "Meteorites.Core/ResultFormatters.h"

namespace
{

constexpr real TIMEOUT      = (real)3600.0;
constexpr real TIME_POINT   = (real)0.9;

std::string ToString(float value)
{
  char buf[1024];
  buf[sprintf_s(buf, "%e", value)] = '\0';
  return std::string(buf);
}

// For each log, determines the number of record that corresponds to normalized 'time'
std::vector<size_t> GetTimePoints(const std::vector<BufferingFormatter::Log> &logs, real time)
{
  auto t_end = std::numeric_limits<real>::max();
  for (const auto &log : logs)
  {
    const auto &records = log.records;
    assert(!records.empty());
    t_end = std::min(t_end, records[records.size() - 1].t);
  }
  assert(t_end > (real)0.0);

  std::vector<size_t> res(logs.size());
  auto t = t_end * time;
  for (size_t i = 0; i < logs.size(); i++)
  {
    const auto &records = logs[i].records;
    size_t lo = 0, hi = records.size();
    while (lo + 1 < hi)
    {
      auto mid = (lo + hi) >> 1;
      if (records[mid].t > t)
      { hi = mid; }
      else
      { lo = mid; }
    }

    assert(lo + 1 == hi);
    if (hi == records.size())
    {
      res[i] = lo;
    }
    else
    {
      res[i] = (t - records[lo].t) > (records[hi].t - t) ? hi : lo;
    }
  }

  return res;
}

std::string ToString(double value)
{
  char buf[1024];
  buf[sprintf_s(buf, "%le", value)] = '\0';
  return std::string(buf);
}

real Eps(real value, real reference)
{
  assert(reference > (real)0.0);
  return std::abs(value - reference) / reference * 100;
}

real ToDegrees(real rads)
{
  return (real)(rads / std::_Pi * 180);
}

} // unnamed namespace

void PrecisionEstimator::CompareMethods(const Case &problem, real dt, std::ostream &str)
{
  assert(dt <= 1e-1f);

  // Prepare data that we want to analyze
  GoldSolver solver;
  BufferingFormatter fmt(0.0f);
  FakeFunctional f;

  solver.Configure(NumericalAlgorithm::ONE_STEP_ADAMS, dt, TIMEOUT);
  solver.Solve(problem, f, fmt);
  solver.Configure(NumericalAlgorithm::TWO_STEP_ADAMS, dt, TIMEOUT);
  solver.Solve(problem, f, fmt);
  solver.Configure(NumericalAlgorithm::THREE_STEP_ADAMS, dt, TIMEOUT);
  solver.Solve(problem, f, fmt);

  decltype(auto) logs = fmt.Logs();
  assert(logs.size() == 3);
  assert(!logs[0].records.empty());
  assert(!logs[1].records.empty());
  assert(!logs[2].records.empty());

  // Print tables
  auto n = GetTimePoints(logs, TIME_POINT);
  std::array<std::string, 3> methods = { "One-step", "Two-step", "Three-step" };
  
  str << "Method,m (kg),V (m/s),h (m), l (m),gamma (deg)," << std::endl;
  for (size_t i = 0; i < methods.size(); i++)
  {
    const auto &record = logs[i].records[n[i]];
    str << methods[i] << ','
        << ToString(record.M) << ','
        << ToString(record.V) << ','
        << ToString(record.h) << ','
        << ToString(record.l) << ','
        << ToString(ToDegrees(record.gamma)) << ','
        << std::endl;
  }
  str << std::endl;
  str << "Method,m (%),V (%),h (%),l (%),gamma (%)," << std::endl;
  for (size_t i = 0; i < methods.size() - 1; i++)
  {
    const auto &record = logs[i].records[n[i]];
    const auto &ref_record = logs[logs.size() - 1].records[n[n.size() - 1]];
    str << methods[i] << ','
        << ToString(Eps(record.M, ref_record.M)) << ','
        << ToString(Eps(record.V, ref_record.V)) << ','
        << ToString(Eps(record.h, ref_record.h)) << ','
        << ToString(Eps(record.l, ref_record.l)) << ','
        << ToString(Eps(record.gamma, ref_record.gamma))
        << ',' << std::endl;
  }
}

void PrecisionEstimator::CompareSteps(const Case &problem,
                                      NumericalAlgorithm method,
                                      std::ostream &str)
{
  // Prepare data that we want to analyze
  FakeFunctional f;
  GoldSolver solver;
  std::array<real, 6> steps{ (real)1e-1, (real)1e-2, (real)1e-3,
                             (real)1e-4, (real)1e-5, (real)1e-6 };

  BufferingFormatter fmt((real)0.0);
  for (real dt : steps)
  {
    solver.Configure(method, dt, TIMEOUT);
    solver.Solve(problem, f, fmt);
  }
  decltype(auto) logs = fmt.Logs();
  assert(logs.size() == steps.size());

  // Print tables
  auto n = GetTimePoints(logs, TIME_POINT);
  str << "dt (s),m (kg),V (m/s),h (m),l (m),gamma (deg)," << std::endl;
  for (size_t i = 0; i < steps.size(); i++)
  {
    const auto &record = logs[i].records[n[i]];
    str << steps[i] << ','
        << ToString(record.M) << ','
        << ToString(record.V) << ','
        << ToString(record.h) << ','
        << ToString(record.l) << ','
        << ToString(ToDegrees(record.gamma)) << ','
        << std::endl;
  }
  str << std::endl;
  str << "dt (s),m (%),V (%),h (%),l (%),gamma (%)," << std::endl;
  for (size_t i = 0; i < steps.size() - 1; i++)
  {
    const auto &record = logs[i].records[n[i]];
    const auto &ref_record = logs[logs.size() - 1].records[n[n.size() - 1]];
    str << steps[i] << ','
        << Eps(record.M, ref_record.M) << ','
        << Eps(record.V, ref_record.V) << ','
        << Eps(record.h, ref_record.h) << ','
        << Eps(record.l, ref_record.l) << ','
        << Eps(record.gamma, ref_record.gamma)
        << ',' << std::endl;
  }
}

void PrecisionEstimator::ComparePerturbations(const Case &p,
                                              NumericalAlgorithm method,
                                              real dt, std::ostream &str)
{
  // Prepare data that we want to analyze
  FakeFunctional f;
  GoldSolver solver;
  solver.Configure(method, dt, TIMEOUT);
  std::array<real, 5> peturbations{ (real)0.01, (real)0.02, (real)0.04, (real)0.08, (real)0.16 };

  std::array<real, 9> params{ p.H(), p.Ch(), p.Rho(), p.Cd(), p.Cl(),
                              p.M0(), p.V0(), p.h0(), p.Gamma0() };
  std::array<std::string, 9> names{ "H", "Ch", "Rho", "Cd", "Cl",
                                    "M0", "V0", "h0", "Gamma0" };

  BufferingFormatter fmt(dt);
  for (size_t i = 0; i < params.size(); i++)
  {
    for (size_t j = 0; j < peturbations.size(); j++)
    {
      auto new_params = params;
      new_params[i] += new_params[i] * peturbations[j];
      Case problem(new_params[0], new_params[1], new_params[2], new_params[3], new_params[4],
                   new_params[5], new_params[6], new_params[7], new_params[8]);
      solver.Solve(problem, f, fmt);
    }
  }
  solver.Solve(p, f, fmt);
  decltype(auto) logs = fmt.Logs();
  assert(logs.size() == params.size() * peturbations.size() + 1);

  // Print tables
  str << "parameter,perturbation (%),t_flight (s),m (kg),V (m/s),h (m),l (m),gamma (deg)," << std::endl;
  for (size_t i = 0; i < params.size(); i++)
  {
    for (size_t j = 0; j < peturbations.size(); j++)
    {
      auto idx = j + i * peturbations.size();
      const auto &records = logs[idx].records;
      const auto &record = records[records.size() - 1];
      str << names[i] << ","
          << peturbations[j] * 100 << ','
          << record.t << ','
          << record.M << ','
          << record.V << ','
          << record.h << ','
          << record.l << ','
          << record.gamma << ','
          << std::endl;
    } // for j
  } // for i
  str << std::endl;

  str << "parameter,perturbation (%),t_flight (%)," << std::endl;
  const auto &ref_log = logs[logs.size() - 1];
  auto ref_time = ref_log.records[ref_log.records.size() - 1].t;
  for (size_t i = 0; i < params.size(); i++)
  {
    for (size_t j = 0; j < peturbations.size(); j++)
    {
      auto idx = j + i * peturbations.size();
      const auto &records = logs[idx].records;
      const auto &record = records[records.size() - 1];
      str << names[i] << ","
          << peturbations[j] * 100 << ','
          << Eps(record.t, ref_time) << ','
          << std::endl;
    } // for j
  } // for i
}
