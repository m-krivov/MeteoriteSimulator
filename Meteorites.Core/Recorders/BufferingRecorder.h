#pragma once
#include "Meteorites.Core/Defs.h"

#include "ISimulationRecorder.h"
#include "MeteoroidTrajectory.h"


// Stores all data that were passed to this recorder
class BufferingRecorder : public ISimulationRecorder
{
  public:
    // Value 'dt' defines how often the incoming records must be stored, in seconds
    BufferingRecorder(real dt) : dt_(dt) { assert(dt >= (real)0.0); }

    // The member of 'ISimulationRecorder'
    virtual real Started(const VirtualMeteoroid &meteoroid) override;

    // The member of 'ISimulationRecorder'
    virtual real Store(real t, real m, real v, real h, real l, real gamma) override;

    // The member of 'ISimulationRecorder'
    virtual bool NeedTrajectory() const override { return true;}
    
    // The member of 'ISimulationRecorder'
    virtual void Finished(Reason reason, double accuracy) override;

    // Information about all trajectories that were passed to this recorder
    const std::vector<MeteoroidTrajectory> &Trajectories() const { return trajectories_; }

    // Deletes information about all recorded trajectories
    void Reset() { trajectories_.clear(); }

    // Moves all recorded trajectories to the specified container
    // After that, resets the internal state
    void MoveTo(std::vector<MeteoroidTrajectory> &trajectories);

  private:
    const real dt_{};
    real t_next_{};
    std::optional<std::tuple<VirtualMeteoroid,
                             ISimulationRecorder::Reason,
                             double,
                             std::vector<MeteoroidTrajectory::Record>>> current_;
    std::vector<MeteoroidTrajectory> trajectories_;
};
