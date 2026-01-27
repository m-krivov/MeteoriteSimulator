#include "HistogramBuilder.h"

#if defined(METEORITES_GNUPLOT)

#define NOMINMAX
#include <matplot/matplot.h>

namespace
{

void PlotHistogram(const std::string &caption,
                   const matplot::axes_handle &ax,
                   std::ofstream &csv,
                   const std::vector<MeteoroidTrajectory> &trajectories,
                   const std::function<real(const MeteoroidTrajectory &)> &get_value,
                   size_t decimical_places)
{
  constexpr size_t N_TICKS = 7;
  constexpr size_t N_BINS_PER_TICK  = 3;
  constexpr size_t N_BINS = (N_TICKS - 1) * N_BINS_PER_TICK;

  assert(ax);
  assert(!trajectories.empty());
  assert(N_BINS_PER_TICK >= 1);
  assert(N_TICKS >= 2);
  assert(N_BINS >= 1);

  // Extract and prepare raw values
  double round_mask = 1.0;
  for (size_t i = 0; i < decimical_places; i++)
  { round_mask /= 10; }
  auto round = [round_mask](real x) -> double
               { return round_mask * std::round(x / round_mask); };

  std::vector<double> values;
  values.reserve(trajectories.size());
  for (const auto &tr : trajectories)
  { values.emplace_back(get_value(tr)); }

  // Construct ticks for histogram, because auto-generated ones are ugly
  double min_value = 0.0, max_value = 0.0;
  {
    auto minmax = std::minmax_element(values.begin(), values.end());
    min_value = *minmax.first  - std::numeric_limits<double>::epsilon();
    max_value = *minmax.second + std::numeric_limits<double>::epsilon();
    if (min_value == max_value)
    {
      min_value -= round_mask * (N_TICKS - 1) / 2;
      max_value += round_mask * (N_TICKS - 1) / 2;
    }
  }
  std::vector<double> tick_positions(N_BINS + 1);
  for (size_t i = 0; i < tick_positions.size(); i++)
  {
    double tick = min_value + (max_value - min_value) * i / (tick_positions.size() - 1);
    tick_positions[i] = tick;
  }

  std::vector<std::string> tick_labels(tick_positions.size());
  for (size_t i = 0; i < tick_positions.size(); i += N_BINS_PER_TICK)
  {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(decimical_places) << round(tick_positions[i]);
    tick_labels[i] = oss.str();
  }

  // Draw them as a histogram
  using namespace matplot;
  auto hist = ax->hist(values, tick_positions);
  hist->normalization(histogram::normalization::probability);
  assert(hist);

  ax->xlabel(caption);
  ax->xlim({ min_value, max_value });
  ax->xticks(tick_positions);
  ax->xticklabels(tick_labels);
  ax->ylabel("Probability");
  ax->draw();

  // Additionally, save these data to a CSV file
  decltype(auto) bins = hist->values();
  assert(tick_positions.size() == bins.size() + 1);
  
  csv << caption << ";probability" << std::endl;
  for (size_t i = 0; i < bins.size(); i++)
  {
    csv << tick_positions[i] << ";" << bins[i] << std::endl;
  }
  csv << tick_positions[bins.size()] << ";" << std::endl;
  csv << std::endl;
}

} // unnamed namespace


void HistogramBuilder::Export(const std::vector<MeteoroidTrajectory> &trajectories)
{
  std::string csv_filename = (Directory() / "parameters.csv").string();
  std::string png_filename = (Directory() / "parameters.png").string();

  std::ofstream csv(csv_filename);
  if (!csv.good())
  {
    std::ostringstream oss;
    oss << "Could not create file '" << csv_filename << "'";
    throw std::runtime_error(oss.str());
  }

  using namespace matplot;
  auto f = figure(true);
  assert(f);
  f->size(1500, 1500);
  f->title(meteorite_->Name());
  {
    PlotHistogram("Density, g/cm^3", subplot(f, 3, 3, 0), csv, trajectories,
                  [](const MeteoroidTrajectory &tr) -> double { return tr.Meteoroid().Rho * 1e-3; }, 1);
    PlotHistogram("Initial mass, kg", subplot(f, 3, 3, 1), csv, trajectories,
                  [](const MeteoroidTrajectory &tr) -> double { return tr.Meteoroid().M0; }, 0);
    PlotHistogram("Initial volume, cm^3", subplot(f, 3, 3, 2), csv, trajectories,
                  [](const MeteoroidTrajectory &tr) -> double
                  { return tr.Meteoroid().M0 / tr.Meteoroid().Rho * 1e3; }, 0);

    PlotHistogram("Entry angle, degrees", subplot(f, 3, 3, 3), csv, trajectories,
                  [](const MeteoroidTrajectory &tr) -> double
                  { return std::round(tr.Meteoroid().Gamma0 * 180.0 / M_PI * 10) / 10; }, 1);
    PlotHistogram("Residual mass, kg", subplot(f, 3, 3, 4), csv, trajectories,
                  [](const MeteoroidTrajectory &tr) -> double
                  { return tr.Reason() == ISimulationRecorder::Reason::Burnt ? 0.0 : tr.LastRecord().M; }, 1);
    PlotHistogram("Flight distance, km", subplot(f, 3, 3, 5), csv, trajectories,
                  [](const MeteoroidTrajectory &tr) -> double { return tr.LastRecord().l * 1e-3; }, 1);

    PlotHistogram("Braking coefficient (Cd)", subplot(f, 3, 3, 6), csv, trajectories,
                  [](const MeteoroidTrajectory &tr) -> double { return tr.Meteoroid().Cd; }, 1);
    PlotHistogram("Lift coefficient (Cl)", subplot(f, 3, 3, 7), csv, trajectories,
                  [](const MeteoroidTrajectory &tr) -> double { return tr.Meteoroid().Cl; }, 2);
  }
  f->draw();
  f->save(png_filename, "png");

  if (!csv.good())
  {
    std::ostringstream oss;
    oss << "Failed to store histogram to '" << csv_filename << "'";
    throw std::runtime_error(oss.str());
  }
}

#endif
