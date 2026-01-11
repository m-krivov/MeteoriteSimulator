#pragma once
#include "Meteorites.Core/Defs.h"

#include "BasicExporter.h"
#include "Meteorites.Core/IMeteorite.h"

class HistogramBuilder : public BasicExporter
{
  public:
    HistogramBuilder() = delete;
    HistogramBuilder(const IMeteorite &meteorite) : meteorite_(meteorite) { }

    // The member of 'IExporter'
    virtual void Export(const std::vector<MeteoroidTrajectory> &trajectories) override final;
    
  private:
    const IMeteorite &meteorite_;
};
