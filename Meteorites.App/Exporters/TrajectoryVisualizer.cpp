#include "TrajectoryVisualizer.h"

#if defined(METEORITES_GNUPLOT)

#define NOMINMAX
#include <matplot/matplot.h>

namespace
{

matplot::vector_2d ToColor(const std::array<double, 3> &color)
{
  return matplot::vector_2d { { color[0], color[1], color[2] } };
};

void Average(const std::vector<double> &arguments,
             const std::vector<std::vector<double>> &values,
             std::vector<double> &average)
{
  average.resize(arguments.size(), 0.0);
  std::vector<size_t> count(arguments.size(), 0);
  for (const auto &v : values)
  {
    for (size_t i = 0; i < v.size(); i++)
    {
      average[i] += v[i];
      count[i]++;
    }
  }
  for (size_t i = 0; i < average.size(); i++)
  {
    if (count[i] == 0)
    {
      average.resize(i);
      break;
    }
    average[i] /= count[i];
  }
}

// Partially vibe-coded by DeepSeek
void Average(const std::vector<std::vector<double>> &arguments,
             const std::vector<std::vector<double>> &values,
             const std::vector<double> &average_arguments,
             std::vector<double> &average_values)
{
  assert(arguments.size() == values.size());
  
  // Number of sample points for the mean trajectory
  const size_t n_points = 200;
  std::vector<double> sampled_arguments(n_points, 0.0);
  std::vector<double> sampled_values(n_points, 0.0);
  std::vector<size_t> counts(n_points, 0);
  
  // For each trajectory, sample it using normalized progress points
  for (size_t j = 0; j < arguments.size(); j++)
  {
    const auto &cur_arguments = arguments[j];
    const auto &cur_values    = values[j];
    if (cur_arguments.size() < 2) continue;
    
    // Sample this trajectory at normalized progress [0, 1]
    for (size_t i = 0; i < n_points; i++)
    {
      double target_distance = (cur_arguments.back() * i) / (n_points - 1);
      
      // Find segment containing target_distance
      auto it = std::lower_bound(cur_arguments.begin(), cur_arguments.end(), target_distance);
      size_t seg_idx = (it == cur_arguments.begin()) ? 0 : std::distance(cur_arguments.begin(), it) - 1;
      
      if (seg_idx + 1 >= cur_arguments.size())
      {
        sampled_arguments[i] += cur_arguments.back();
        sampled_values[i] += cur_values.back();
        counts[i]++;
        break;
      }
      
      // Linear interpolation within segment
      double d1 = cur_arguments[seg_idx];
      double d2 = cur_arguments[seg_idx + 1];
      double h1 = cur_values[seg_idx];
      double h2 = cur_values[seg_idx + 1];
      
      if (std::abs(d2 - d1) > std::numeric_limits<double>::epsilon() * 10)
      {
        double t = (target_distance - d1) / (d2 - d1);
        double interpolated = h1 + t * (h2 - h1);
        
        sampled_arguments[i] += target_distance;
        sampled_values[i] += interpolated;
        counts[i]++;
      }
      else
      {
        sampled_arguments[i] += d1;
        sampled_values[i] += h1;
        counts[i]++;
      }
    }
  }
  
  // Calculate the mean values and create the mean trajectory
  for (size_t i = 0; i < n_points; i++)
  {
    if (counts[i] > 0)
    {
      sampled_arguments[i] /= counts[i];
      sampled_values[i] /= counts[i];
    }
  }

  // Finally, map sampled averages to the requested arguments
  average_values.clear();
  average_values.reserve(average_arguments.size());
  
  for (size_t i = 0; i < average_arguments.size(); i++)
  {
    double target_arg = average_arguments[i];
    
    size_t seg_idx = 0;
    while (seg_idx + 1 < n_points && sampled_arguments[seg_idx + 1] < target_arg)
    { seg_idx++; }
    
    //  Just interrupt average trajectory if the requested arguments do not match it
    if (seg_idx + 1 >= n_points)
    { break; }
    
    if (target_arg <= sampled_arguments[0])
    {
      average_values.push_back(sampled_values[0]);
      continue;
    }
    
    // Linear interpolation
    double d1 = sampled_arguments[seg_idx];
    double d2 = sampled_arguments[seg_idx + 1];
    double h1 = sampled_values[seg_idx];
    double h2 = sampled_values[seg_idx + 1];
    
    if (std::abs(d2 - d1) > std::numeric_limits<double>::epsilon() * 10)
    {
      double t = (target_arg - d1) / (d2 - d1);
      double interpolated = h1 + t * (h2 - h1);
      average_values.push_back(interpolated);
    }
    else
    {
      average_values.push_back(h1);
    }
  }
}

// Partially vibe-coded by DeepSeek
void Resample(const std::vector<std::vector<double>> &arguments,
              const std::vector<std::vector<double>> &values,
              std::vector<double> &resampled_arguments,
              std::vector<std::vector<double>> &resampled_values)
{
  assert(arguments.size() == values.size());
  if (arguments.empty() || values.empty())
  { return; }

  // Step 1: Find the global maximum across all arguments
  double max_arg = std::numeric_limits<double>::lowest();  
  for (const auto &r : arguments)
  {
    if (r.empty())
    { continue; }
    
    assert(r[0] == 0.0);
    max_arg = std::max(max_arg, r.back());
  }

  // Step 2: Create unified vector with arguments
  // Use the maximum number of points from all trajectories
  {
    size_t max_points = 0;
    for (const auto &r : arguments)
    { max_points = std::max(max_points, r.size()); }
    
    resampled_arguments.reserve(max_points);
    double step = max_arg / (max_points - 1);
    for (size_t i = 0; i < max_points; ++i)
    { resampled_arguments.push_back(i * step); }
  }

  // Step 3: Resample arguments and values
  for (size_t i = 0; i < arguments.size(); i++)
  {
    const auto &cur_arguments = arguments[i];
    const auto &cur_values = values[i];
    assert(cur_arguments.size() == cur_values.size());
    
    std::vector<double> cur_resampled_values;
    cur_resampled_values.reserve(values.size());
    
    for (double arg : resampled_arguments)
    {
      auto it = std::lower_bound(cur_arguments.begin(), cur_arguments.end(), arg);
      if (it == cur_arguments.end())
      { break; }
      
      if (it != cur_arguments.begin())
      {
        size_t idx = std::distance(cur_arguments.begin(), it);
        double prev_arg = cur_arguments[idx - 1];
        double next_arg = cur_arguments[idx];
        double prev_val = cur_values[idx - 1];
        double next_val = cur_values[idx];

        if (std::abs(next_arg - prev_arg) <= std::numeric_limits<double>::epsilon() * 10)
        { cur_resampled_values.push_back(prev_val); }
        else
        {
          double t = (arg - prev_arg) / (next_arg - prev_arg);
          double interpolated = prev_val + t * (next_val - prev_val);
          cur_resampled_values.push_back(interpolated);
        }
      }
      else
      { cur_resampled_values.push_back(cur_values.front()); }
    }
    resampled_values.emplace_back(std::move(cur_resampled_values));
  }
}

// Contains data for visualization h(t), V(t), M(t) and h(l)
struct VisualizationData
{
  std::vector<double> arguments;
  std::vector<std::vector<double>> values;
  std::vector<double> average;
  std::vector<std::pair<double, double>> points;

  std::array<double, 3> value_color{ 0.75, 0.75, 0.75 };
  std::array<double, 3> average_color{ 0.19, 0.18, 0.17 };
  std::array<double, 3> point_color{ 0.82, 0.32, 0.13 };

  float average_width{ 1.5f };
  size_t point_size{ 4 };
  
  // Version for h(t), V(t), M(t)
  VisualizationData(const std::vector<double> &arguments_,
                    std::vector<std::vector<double>> &&values_,
                    std::vector<std::pair<double, double>> &&points_)
    : arguments(arguments_), values(std::move(values_)), points(std::move(points_))
  {
    Average(arguments, values, average);
  }

  // Version for h(l)
  VisualizationData(const std::vector<std::vector<double>> &arguments_,
                    const std::vector<std::vector<double>> &values_,
                    std::vector<std::pair<double, double>> &&points_)
    : points(std::move(points_))
  {
    Resample(arguments_, values_, arguments, values);
    Average(arguments_, values_, arguments, average);
  }

  VisualizationData(const VisualizationData &) = delete;
  VisualizationData(VisualizationData &&) = delete;
};

void Visualize(const matplot::axes_handle &ax,
               const VisualizationData &tr)
{
  assert(!tr.arguments.empty());
  assert(!tr.values.empty());
  assert(!tr.values[0].empty());

  using namespace matplot;
  grid(ax, on);

  // The original trajectories
  colororder(ax, ToColor(tr.value_color));
  plot(ax, tr.arguments, tr.values);

  // Average trajectory
  {
    hold(ax, on);
    colororder(ax, ToColor(tr.average_color));
    auto l = plot(ax, tr.arguments, { tr.average });
    l->line_width(tr.average_width);
    l->line_style("-");
    hold(ax, off);
  }
  
  // Points for reference values (if available)
  if (!tr.points.empty())
  {
    hold(ax, on);
    colororder(ax, ToColor(tr.point_color));
    std::vector<double> point_x, point_y;
    for (const auto &p : tr.points)
    {
      point_x.emplace_back(p.first);
      point_y.emplace_back(p.second);
    }
    auto l = scatter(ax, point_x, point_y, tr.point_size);
    l->marker_face(true);
    hold(ax, off);
  }
}

void Store(const std::filesystem::path &filename,
           const VisualizationData &tr)
{
  std::ofstream f(filename);

  // Titles
  f << "Argument;Mean;";
  for (size_t i = 0; i < tr.values.size(); i++)
  { f << "Value #" << i << ';'; }
  f << "Reference value" << std::endl;

  // Values themselves
  size_t next_point = 0;
  for (size_t i = 0; i < tr.arguments.size(); i++)
  {
    f << tr.arguments[i] << ';';
    if (i < tr.average.size())
    { f << tr.average[i]; }
    f << ';';

    for (const auto &value : tr.values)
    {
      if (i < value.size())
      { f << value[i]; }
      f << ';';
    }

    if (next_point < tr.points.size() && tr.arguments[i] >= tr.points[next_point].first)
    {
      f << tr.points[next_point].second;
      next_point += 1;
    }
    f << ';' << std::endl;
  }

  // Final check
  if (!f.good())
  {
    std::ostringstream oss;
    oss << "Failed to create file '" << filename << "'";
    throw std::runtime_error(oss.str());
  }
}

} // unnamed namespace


void TrajectoryVisualizer::Export(const std::vector<MeteoroidTrajectory> &trajectories)
{
  assert(trajectories.size() > 0);

  // Extract the reference values for velocity and height (will be presented as points)
  real t_end = (real)0.0f;
  std::vector<std::pair<double, double>> ref_velocity, ref_height;
  {
    size_t n_records = 0;
    const real *tab_time = nullptr, *tab_velocity = nullptr, *tab_height = nullptr;
    meteorite_->Trajectory(n_records, tab_time, tab_velocity, tab_height);
    assert(n_records > 0);
    t_end = tab_time[n_records - 1] * multiplier_;
    
    for (size_t i = 0; i < n_records; i++)
    {
      ref_velocity.emplace_back(std::make_pair(tab_time[i], tab_velocity[i] * 1e-3));
      ref_height.emplace_back(std::make_pair(tab_time[i], tab_height[i] * 1e-3));
    }
  }

  // Build the Ox axis for time
  std::vector<double> time_axis;
  {
    auto iter = std::max_element(trajectories.begin(), trajectories.end(),
                [](const MeteoroidTrajectory &t1, const MeteoroidTrajectory &t2) -> bool
                { return t1.Records().size() < t2.Records().size(); });
    assert(iter != trajectories.end());
    time_axis.reserve(iter->Records().size());
    for (const auto &record : iter->Records())
    {
      if (record.t > t_end + std::numeric_limits<real>::epsilon())
      { break; }
      time_axis.emplace_back(record.t);
    }
  }

  // Convert trajectories to a format expected by Matplot++
  std::vector<std::vector<double>> tr_partial_velocity(trajectories.size()),
                                   tr_partial_height(trajectories.size()),
                                   tr_partial_mass(trajectories.size());
  for (size_t i = 0; i < trajectories.size(); i++)
  {
    for (const auto &record : trajectories[i].Records())
    {
      if (record.t > t_end + std::numeric_limits<real>::epsilon())
      { break; }

      tr_partial_velocity[i].emplace_back(record.V * 1e-3);
      tr_partial_height[i].emplace_back(record.h * 1e-3);
      tr_partial_mass[i].emplace_back(record.M);
    }
  }
  std::vector<std::vector<double>> tr_full_distance(trajectories.size()),
                                   tr_full_height(trajectories.size());
  for (size_t i = 0; i < trajectories.size(); i++)
  {
    tr_full_distance[i].reserve(trajectories[i].Records().size());
    tr_full_height[i].reserve(trajectories[i].Records().size());
    for (const auto &record : trajectories[i].Records())
    {
      tr_full_distance[i].emplace_back(record.l * 1e-3);
      tr_full_height[i].emplace_back(record.h * 1e-3);
    }
  }

  // Draw the figures, save them to disk
  using namespace matplot;
  auto f = figure();
  assert(f);

  f->size(1200, 1200);
  f->title(meteorite_->Name());

  {
    auto ax = subplot(f, 3, 2, 0);
    assert(ax);
    xlabel(ax, "Time, seconds");
    xlim(ax, { 0.0, t_end });
    ylabel(ax, "Velocity, km/s");

    VisualizationData tr(time_axis,
                         std::move(tr_partial_velocity),
                         std::move(ref_velocity));
    Visualize(ax, tr);
    Store(Directory() / "trajectories_velocity.csv", tr);
  }
  UpdateProgress(1, 4);

  {
    auto ax = subplot(f, 3, 2, 1);
    assert(ax);
    xlabel(ax, "Time, seconds");
    xlim(ax, { 0.0, t_end });
    ylabel(ax, "Height, km");

    VisualizationData tr(time_axis,
                         std::move(tr_partial_height),
                         std::move(ref_height));
    Visualize(ax, tr);
    Store(Directory() / "trajectories_height.csv", tr);
  }
  UpdateProgress(2, 4);

  {
    auto ax = subplot(f, 3, 2, 2);
    assert(ax);
    xlabel(ax, "Time, seconds");
    xlim(ax, { 0.0, t_end });
    ylabel(ax, "Mass, kg");

    VisualizationData tr(time_axis, std::move(tr_partial_mass), {});
    Visualize(ax, tr);
    Store(Directory() / "trajectories_mass.csv", tr);
  }
  UpdateProgress(3, 4);
  
  {
    auto ax = f->add_subplot(3, 2, { 4, 5 });
    assert(ax);
    xlabel(ax, "Distance, km");
    ylabel(ax, "Height, km");
    VisualizationData tr(tr_full_distance, tr_full_height, {});
    Visualize(ax, tr);
    Store(Directory() / "trajectories_distance.csv", tr);
  }
  UpdateProgress(4, 4);

  f->draw();
  f->save((Directory() / "trajectories.png").string(), "png");
  
  // Save functional structure information if available
  if (functional_ != nullptr)
  {
    std::ofstream f_info(Directory() / "functional_info.txt");
    if (f_info.good())
    {
      f_info << "Functional Structure Information" << std::endl;
      f_info << "================================" << std::endl << std::endl;
      f_info << functional_->GetStructureDescription() << std::endl;
    }
  }
}

#endif
