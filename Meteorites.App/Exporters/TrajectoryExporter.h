#pragma once
#include "Meteorites.Core/Defs.h"

#include "BasicExporter.h"

// Saves all provided trajectories to disk, creates one *.csv file per trajectory
class TrajectoryExporter : public BasicExporter
{
  public:
    TrajectoryExporter() = delete;
    TrajectoryExporter(real dt) : dt_(dt) { }

    // The member of 'IExporter'
    virtual void Export(const std::vector<MeteoroidTrajectory> &trajectories) override final;

  private:
    real dt_ = (real)0.0;
};
