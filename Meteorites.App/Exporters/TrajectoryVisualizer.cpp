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

// Contains data for visualization h(t), V(t) and M(t)
struct TimeDependentTrajectories
{
  std::vector<double> times;
  std::vector<std::vector<double>> values;
  std::vector<std::pair<double, double>> points;

  std::array<double, 3> line_color{ 0.75, 0.75, 0.75 };
  std::array<double, 3> average_color{ 0.19, 0.18, 0.17 };
  std::array<double, 3> point_color{ 0.82, 0.32, 0.13 };

  float average_width{ 1.5f };
  size_t point_size{ 4 };
  
  TimeDependentTrajectories(const std::vector<double> &times_,
                            std::vector<std::vector<double>> &&values_,
                            std::vector<std::pair<double, double>> &&points_)
    : times(times_), values(std::move(values_)), points(std::move(points_))
  { }
  TimeDependentTrajectories(const TimeDependentTrajectories &) = delete;
  TimeDependentTrajectories(TimeDependentTrajectories &&) = delete;
};

void Visualize(const matplot::axes_handle &ax, const TimeDependentTrajectories &tr)
{
  assert(!tr.times.empty());
  assert(!tr.values.empty());
  assert(!tr.values[0].empty());

  using namespace matplot;
  grid(ax, on);

  // The original trajectories
  colororder(ax, ToColor(tr.line_color));
  plot(ax, tr.times, tr.values);

  // Average trajectory
  {
    std::vector<double> average(tr.times.size(), 0.0);
    std::vector<size_t> count(tr.times.size(), 0);
    for (const auto &v : tr.values)
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
    hold(ax, on);
    colororder(ax, ToColor(tr.average_color));
    auto l = plot(ax, tr.times, { std::move(average) });
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

// More traditional trajectories that represent height-distance dependencies
struct DistanceHeightTrajectories
{
  std::vector<std::vector<double>> distances;
  std::vector<std::vector<double>> heights;

  std::array<double, 3> line_color{ 0.75, 0.75, 0.75 };
  std::array<double, 3> average_color{ 0.19, 0.18, 0.17 };
  float average_width{ 1.5f };

  DistanceHeightTrajectories(std::vector<std::vector<double>> &&distances_,
                             std::vector<std::vector<double>> &&heights_)
    : distances(std::move(distances_)), heights(std::move(heights_))
  { }
  DistanceHeightTrajectories(const DistanceHeightTrajectories &) = delete;
  DistanceHeightTrajectories(DistanceHeightTrajectories &&) = delete;
};

void Visualize(const matplot::axes_handle &ax, const DistanceHeightTrajectories &tr)
{
  assert(tr.distances.size() > 1);   // for interpolation
  assert(tr.distances.size() == tr.heights.size());

  using namespace matplot;
  grid(ax, on);
  hold(ax, on);

  // Visualize individual trajectories
  colororder(ax, ToColor(tr.line_color));

  for (size_t i = 0; i < tr.distances.size(); i++)
  {
    const auto &d = tr.distances[i];
    const auto &h = tr.heights[i];
    assert(!d.empty());
    assert(d.size() == h.size());
    plot(ax, d, h);
  }

  // Vibe-coded by DeepSeek
  {
    // First, find the maximum total distance across all trajectories
    double max_total_distance = 0.0;
    for (const auto &d_vec : tr.distances)
    {
      if (!d_vec.empty())
      {
        max_total_distance = std::max(max_total_distance, d_vec.back());
      }
    }
    
    // Number of sample points for the mean trajectory
    const size_t num_mean_points = 200;
    std::vector<double> mean_distances(num_mean_points, 0.0);
    std::vector<double> mean_heights(num_mean_points, 0.0);
    std::vector<size_t> mean_counts(num_mean_points, 0);
    
    // For each trajectory, sample at normalized progress points
    for (size_t traj_idx = 0; traj_idx < tr.distances.size(); traj_idx++)
    {
      const auto &dist_vec = tr.distances[traj_idx];
      const auto &height_vec = tr.heights[traj_idx];
      
      if (dist_vec.size() < 2) continue;
      
      double traj_total_distance = dist_vec.back();
      
      // Sample this trajectory at normalized progress [0, 1]
      for (size_t i = 0; i < num_mean_points; i++)
      {
        double progress = static_cast<double>(i) / (num_mean_points - 1);
        
        // Scale progress by this trajectory's total distance
        double target_distance = progress * traj_total_distance;
        
        // Find segment containing target_distance
        size_t seg_idx = 0;
        while (seg_idx + 1 < dist_vec.size() && dist_vec[seg_idx + 1] < target_distance)
        {
          seg_idx++;
        }
        
        if (seg_idx + 1 >= dist_vec.size())
        {
          // Use last point if we're beyond the trajectory
          mean_distances[i] += dist_vec.back();
          mean_heights[i] += height_vec.back();
          mean_counts[i]++;
          break;
        }
        
        // Linear interpolation within segment
        double d1 = dist_vec[seg_idx];
        double d2 = dist_vec[seg_idx + 1];
        double h1 = height_vec[seg_idx];
        double h2 = height_vec[seg_idx + 1];
        
        if (std::abs(d2 - d1) > std::numeric_limits<double>::epsilon())
        {
          double t = (target_distance - d1) / (d2 - d1);
          double interpolated_height = h1 + t * (h2 - h1);
          
          mean_distances[i] += target_distance;  // Use the target distance
          mean_heights[i] += interpolated_height;
          mean_counts[i]++;
        }
        else
        {
          // Zero-length segment, use the point
          mean_distances[i] += d1;
          mean_heights[i] += h1;
          mean_counts[i]++;
        }
      }
    }
    
    // Calculate the mean values and create the mean trajectory
    std::vector<double> final_mean_distances;
    std::vector<double> final_mean_heights;
    
    for (size_t i = 0; i < num_mean_points; i++)
    {
      if (mean_counts[i] > 0)
      {
        final_mean_distances.push_back(mean_distances[i] / mean_counts[i]);
        final_mean_heights.push_back(mean_heights[i] / mean_counts[i]);
      }
    }
    
    // Plot the mean trajectory
    if (!final_mean_distances.empty() && final_mean_distances.size() > 1)
    {
      colororder(ax, ToColor(tr.average_color));
      auto mean_line = plot(ax, final_mean_distances, final_mean_heights);
      mean_line->line_width(tr.average_width);
      mean_line->line_style("-");
    }
  }

  hold(ax, off);
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

    TimeDependentTrajectories tr(time_axis,
                                 std::move(tr_partial_velocity),
                                 std::move(ref_velocity));
    Visualize(ax, tr);
  }
  UpdateProgress(1, 4);

  {
    auto ax = subplot(f, 3, 2, 1);
    assert(ax);
    xlabel(ax, "Time, seconds");
    xlim(ax, { 0.0, t_end });
    ylabel(ax, "Height, km");

    TimeDependentTrajectories tr(time_axis,
                                 std::move(tr_partial_height),
                                 std::move(ref_height));
    Visualize(ax, tr);
  }
  UpdateProgress(2, 4);

  {
    auto ax = subplot(f, 3, 2, 2);
    assert(ax);
    xlabel(ax, "Time, seconds");
    xlim(ax, { 0.0, t_end });
    ylabel(ax, "Mass, kg");

    TimeDependentTrajectories tr(time_axis, std::move(tr_partial_mass), {});
    Visualize(ax, tr);
  }
  UpdateProgress(3, 4);
  
  {
    auto ax = f->add_subplot(3, 2, { 4, 5 });
    assert(ax);
    xlabel(ax, "Distance, km");
    ylabel(ax, "Height, km");
    DistanceHeightTrajectories tr(std::move(tr_full_distance),
                                  std::move(tr_full_height));
    Visualize(ax, tr);
  }
  UpdateProgress(4, 4);

  f->draw();
  f->save((Directory() / "trajectories.png").string(), "png");
}

#endif
