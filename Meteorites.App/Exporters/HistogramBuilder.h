#pragma once
#include "Meteorites.Core/Defs.h"

#include "BasicExporter.h"
#include "Meteorites.Core/IMeteorite.h"

class HistogramBuilder : public BasicExporter
{
  public:
    HistogramBuilder() = delete;
    HistogramBuilder(const std::shared_ptr<const IMeteorite> &meteorite) : meteorite_(meteorite) { }

    // The member of 'IExporter'
    virtual void Export(const std::vector<MeteoroidTrajectory> &trajectories) override final;
    
  private:
    const std::shared_ptr<const IMeteorite> &meteorite_;
};
