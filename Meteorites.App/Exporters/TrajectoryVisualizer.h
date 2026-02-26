#pragma once
#include "Meteorites.Core/Defs.h"

#include "BasicExporter.h"
#include "Meteorites.Core/IMeteorite.h"
#include "Meteorites.Core/Functionals/IFunctional.h"

class TrajectoryVisualizer : public BasicExporter
{
  public:
    TrajectoryVisualizer() = delete;
    TrajectoryVisualizer(const std::shared_ptr<const IMeteorite> &meteorite, real multiplier = 2)
      : meteorite_(meteorite), multiplier_(multiplier)
    { }

    // Optional: Set the functional to include its structure in visualization
    void SetFunctional(const IFunctional *functional)
    { functional_ = functional; }

    // The member of 'IExporter'
    virtual void Export(const std::vector<MeteoroidTrajectory> &trajectories) override final;

  private:
    const std::shared_ptr<const IMeteorite> &meteorite_;
    const IFunctional *functional_ = nullptr;
    real multiplier_ = (real)1.0;
};
