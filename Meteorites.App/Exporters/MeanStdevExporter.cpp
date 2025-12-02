#include "MeanStdevExporter.h"

namespace
{

constexpr char SEPARATOR = ';';
constexpr real RAD_TO_DEG = (real)(180.0 / M_PI);

std::pair<double, double> CalculateMeanAndStdev(const std::vector<real>::iterator &begin,
                                                const std::vector<real>::iterator &end)
{
  if (begin == end)
  { return std::make_pair(0.0, 0.0); }
  
  size_t size = std::distance(begin, end);
  if (size == 1)
  { return std::make_pair(*begin, 0.0); }

  double sum = 0.0;
  for (auto iter = begin; iter != end; iter++)
  { sum += *iter; }
  double mean = sum / size;

  double acc = 0.0;
  for (auto iter = begin; iter != end; iter++)
  {
    double diff = *iter - mean;
    acc += diff * diff;
  }
  double stdev = std::sqrt(acc / (size - 1));

  return std::make_pair(mean, stdev);
}

void ExportParameter(std::ofstream &mean_file,
                     std::ofstream &mean_stdev_file,
                     const std::string &name,
                     const std::vector<size_t> &group_sizes,
                     const std::vector<MeteoroidTrajectory> &trajectories,
                     const std::function<real(const MeteoroidTrajectory &)> &get_value)
{
  assert(mean_file.good());
  assert(mean_stdev_file.good());
  assert(!group_sizes.empty());

  std::vector<real> values;
  for (const auto &trajectory : trajectories)
  { values.emplace_back(get_value(trajectory)); }

  mean_file << name << SEPARATOR;
  mean_stdev_file << name << SEPARATOR;
  for (auto group_size : group_sizes)
  {
    assert(group_size <= values.size());
    auto mean_stdev = CalculateMeanAndStdev(values.begin(), values.begin() + group_size);
    mean_file << mean_stdev.first << SEPARATOR;
    mean_stdev_file << mean_stdev.first << "+-" << mean_stdev.second << SEPARATOR;
  }
  mean_file << std::endl;
  mean_stdev_file << std::endl;
}

} // unnamed namespace


void MeanStdevExporter::Export(const std::vector<MeteoroidTrajectory> &trajectories)
{
  if (trajectories.empty())
  { return; }
  
  std::ofstream mean_file(Directory() / "mean.csv", std::ios::out);
  std::ofstream mean_stdev_file(Directory() / "mean_stdev.csv", std::ios::out);

  // Split all meteoroids into subsets of different size
  std::vector<size_t> group_sizes(groups_);
  for (size_t i = 0; i < group_sizes.size(); i++)
  {
    group_sizes[i] = (size_t)(trajectories.size() * ((i + 1) / (double)groups_));
    group_sizes[i] = std::max((size_t)1, group_sizes[i]);
  }
  assert(groups_ >= 1);
  group_sizes[groups_ - 1] = trajectories.size();

  // Print the header of CSV table
  mean_file << "Parameter name" << SEPARATOR;
  mean_stdev_file << "Parameter name" << SEPARATOR;
  for (auto group_size : group_sizes)
  {
    mean_file << "Top" << group_size << SEPARATOR;
    mean_stdev_file << "Top" << group_size << SEPARATOR;
  }
  mean_file << std::endl;
  mean_stdev_file << std::endl;

  // Print mean and standard deviation for each important parameter
  std::vector<std::pair<std::string, std::function<real(const MeteoroidTrajectory &)>>> parameters =
  {
    { "Initial mass, kg", [](const MeteoroidTrajectory &m) -> real { return m.Meteoroid().M0; } },
    { "Residual mass, kg", [](const MeteoroidTrajectory &m) -> real { return m.LastRecord().M; } },
    { "Density, kg/m^3", [](const MeteoroidTrajectory &m) -> real { return m.Meteoroid().Rho; } },
    { "Flight distance, m", [](const MeteoroidTrajectory &m) -> real { return m.LastRecord().l; } },
    { "Entry angle, degrees", [](const MeteoroidTrajectory &m) -> real { return m.Meteoroid().Gamma0 * RAD_TO_DEG; } },
    { "Enthalpy of destruction, J/kg", [](const MeteoroidTrajectory &m) -> real { return m.Meteoroid().H; } },
    { "Drag force coefficient", [](const MeteoroidTrajectory &m) -> real { return m.Meteoroid().Cd; } },
    { "Lift force coefficient", [](const MeteoroidTrajectory &m) -> real { return m.Meteoroid().Cl; } },
    { "Heat transfer coefficient", [](const MeteoroidTrajectory &m) -> real { return m.Meteoroid().Ch; } }
  };

  // In addition, classify meteoroids by the reason why simulation was ended
  parameters.push_back
  ({
    "Probability of combustion, %", [](const MeteoroidTrajectory &m) -> real
    { return m.Reason() == ISimulationRecorder::Reason::Burnt ? (real)100.0 : (real)0.0; }
  });
  parameters.push_back
  ({
    "Probability of collision, %", [](const MeteoroidTrajectory &m) -> real
    { return m.Reason() == ISimulationRecorder::Reason::Collided ? (real)100.0 : (real)0.0; }
  });

  for (size_t i = 0; i < parameters.size(); i++)
  {
    ExportParameter(mean_file, mean_stdev_file, parameters[i].first,
                    group_sizes, trajectories, parameters[i].second);
    UpdateProgress(i + 1, parameters.size());
  }

  if (mean_file.fail() || mean_stdev_file.fail())
  {
    std::ostringstream oss;
    oss << "Failed to save mean+-stdev values to the directory '" << Directory() << '\'';
    throw std::runtime_error(oss.str());
  }
}
