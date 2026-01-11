#include "HistogramBuilder.h"

#if defined(METEORITES_GNUPLOT)

#define NOMINMAX
#include <matplot/matplot.h>

namespace
{

void PlotHistogram(const matplot::axes_handle &ax,
                   const std::string &caption,
                   const std::vector<MeteoroidTrajectory> &trajectories,
                   const std::function<real(const MeteoroidTrajectory &)> &get_value)
{
  assert(ax);
  using namespace matplot;
  std::vector<double> values;
  values.reserve(trajectories.size());
  for (const auto &tr : trajectories)
  { values.emplace_back(get_value(tr)); }

  auto hist = ax->hist(values, histogram::normalization::probability);
  assert(hist);

  hist->num_bins(15);
  ax->xlabel(caption);
  ax->ylabel("Probability");
}

} // unnamed namespace


void HistogramBuilder::Export(const std::vector<MeteoroidTrajectory> &trajectories)
{
  using namespace matplot;
  auto f = figure(true);
  assert(f);

  f->size(1500, 1500);
  f->title(meteorite_.Name());
  {
    PlotHistogram(subplot(f, 3, 3, 0), "Density, g/cm^3", trajectories,
                  [](const MeteoroidTrajectory &tr) -> double { return tr.Meteoroid().Rho * 1e-3; });
    PlotHistogram(subplot(f, 3, 3, 1), "Initial mass, kg", trajectories,
                  [](const MeteoroidTrajectory &tr) -> double { return tr.Meteoroid().M0; });
    PlotHistogram(subplot(f, 3, 3, 2), "Initial volume, cm^3", trajectories,
                  [](const MeteoroidTrajectory &tr) -> double
                  { return tr.Meteoroid().M0 / tr.Meteoroid().Rho * 1e3; });

    PlotHistogram(subplot(f, 3, 3, 3), "Entry angle, degrees", trajectories,
                  [](const MeteoroidTrajectory &tr) -> double { return tr.Meteoroid().Gamma0 * 180.0 / M_PI; });
    PlotHistogram(subplot(f, 3, 3, 4), "Residual mass, kg", trajectories,
                  [](const MeteoroidTrajectory &tr) -> double { return tr.LastRecord().M; });
    PlotHistogram(subplot(f, 3, 3, 5), "Flight distance, km", trajectories,
                  [](const MeteoroidTrajectory &tr) -> double { return tr.LastRecord().l * 1e-3; });

     PlotHistogram(subplot(f, 3, 3, 6), "Braking coefficient (Cd)", trajectories,
                  [](const MeteoroidTrajectory &tr) -> double { return tr.Meteoroid().Cd; });
     PlotHistogram(subplot(f, 3, 3, 7), "Lift coefficient (Cl)", trajectories,
                  [](const MeteoroidTrajectory &tr) -> double { return tr.Meteoroid().Cl; });
  }
  f->draw();
  f->save((Directory() / "parameters.png").string(), "png");
}

#endif
