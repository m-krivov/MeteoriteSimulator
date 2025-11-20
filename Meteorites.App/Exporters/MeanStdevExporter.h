#pragma once
#include "Meteorites.Core/Defs.h"

#include "BasicExporter.h"

// Averages parameters of the simulated meteoroids and stores them as mean and standard deviation
// For each parameter, creates a separate *.csv file
class MeanStdevExporter : public BasicExporter
{
  public:
    MeanStdevExporter() = delete;
    MeanStdevExporter(size_t groups) : groups_(groups) { assert(groups >= 1); }

    // The member of 'IExporter'
    virtual void Export(const std::vector<MeteoroidTrajectory> &trajectories) override final;

  private:
    size_t groups_{1};
};
