#pragma once
#include "Meteorites.Core/Defs.h"

#include "BasicExporter.h"
#include "Meteorites.Core/IMeteorite.h"

class TrajectoryVisualizer : public BasicExporter
{
  public:
    TrajectoryVisualizer() = delete;
    TrajectoryVisualizer(const std::shared_ptr<const IMeteorite> &meteorite, real multiplier = 2)
      : meteorite_(meteorite), multiplier_(multiplier)
    { }

    // The member of 'IExporter'
    virtual void Export(const std::vector<MeteoroidTrajectory> &trajectories) override final;

  private:
    const std::shared_ptr<const IMeteorite> &meteorite_;
    real multiplier_ = (real)1.0;
};
