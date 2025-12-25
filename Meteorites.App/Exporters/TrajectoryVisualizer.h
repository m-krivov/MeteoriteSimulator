#pragma once
#include "Meteorites.Core/Defs.h"

#include "BasicExporter.h"

class TrajectoryVisualizer : public BasicExporter
{
  public:
    TrajectoryVisualizer() = delete;
    TrajectoryVisualizer(real dt) : dt_(dt) { }

    // The member of 'IExporter'
    virtual void Export(const std::vector<MeteoroidTrajectory> &trajectories) override final;

  private:
    real dt_ = (real)0.0;
};
