#include "TrajectoryExporter.h"

#include "Meteorites.Core/Recorders/CsvRecorder.h"


void TrajectoryExporter::Export(const std::vector<MeteoroidTrajectory> &trajectories)
{
  std::filesystem::path dir = Directory() / "trajectories";
  if (!std::filesystem::create_directories(dir))
  {
    std::ostringstream oss;
    oss << "Failed to create directory '" << dir << '\'';
    throw std::runtime_error(oss.str());
  }

  CsvRecorder recorder(dir, Meteorite().Name(), dt_);
  for (size_t i = 0; i < trajectories.size(); i++)
  {
    trajectories[i].ExportTo(recorder);
    UpdateProgress(i + 1, trajectories.size());
  }
}
