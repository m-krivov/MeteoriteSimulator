#include "TrajectoryExporter.h"

#include "Meteorites.Core/Recorders/CsvRecorder.h"

void TrajectoryExporter::Export(const std::vector<MeteoroidTrajectory> &trajectories)
{
  CsvRecorder recorder(Directory(), Meteorite().Name(), dt_);
  for (size_t i = 0; i < trajectories.size(); i++)
  {
    const auto &trajectory = trajectories[i];

    recorder.Started(trajectory.Meteoroid());
    for (const auto &point : trajectory.Records())
    { recorder.Store(point.t, point.M, point.V, point.h, point.l, point.gamma); }
    recorder.Finished(trajectory.Reason(), trajectory.Accuracy());

    UpdateProgress(i + 1, trajectories.size());
  }
}
