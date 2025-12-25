#include "TrajectoryExporter.h"

#include "Meteorites.Core/Recorders/CsvRecorder.h"


void TrajectoryExporter::Export(const std::vector<MeteoroidTrajectory> &trajectories)
{
  CsvRecorder recorder(Directory(), Meteorite().Name(), dt_);
  for (size_t i = 0; i < trajectories.size(); i++)
  {
    trajectories[i].ExportTo(recorder);
    UpdateProgress(i + 1, trajectories.size());
  }
}
