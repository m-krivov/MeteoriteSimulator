#pragma once
#include "Meteorites.Core/Defs.h"

#include "Meteorites.Core/Meteoroids/VirtualMeteoroid.h"


// Interface that saves the numerically computed values
class ISimulationRecorder
{
  public:
    // Due to which reason the simulation was stopped?
    enum class Reason : uint32_t
    {
      // Field is not initialized
      NA          = 0,

      // No mass left, meteorite is burnt
      Burnt       = 1,

      // Height is equal or lesser than zero, B-O-O-O-M!
      Collided    = 2,

      // Looks like something is wrong, meteorite flies for too long time
      // Did we process some UFO trajectory instead of meteorite?
      Timeouted   = 3
    };

    // A helper that determines, why do we stop simulation of a meteorite
    // TODO: implement as a second version of 'Finished()'
    static Reason Classify(real t, real M, real h);

    ISimulationRecorder(const ISimulationRecorder &) = delete;
    ISimulationRecorder &operator =(const ISimulationRecorder &) = delete;
    virtual ~ISimulationRecorder() { }

    // Notifies that simulation is started, data will be coming soon
    // Provides problem that are going to be simulated
    // Returns the expected time of the first record (take a look at 'Store()')
    virtual real Started(const VirtualMeteoroid &problem) = 0;

    // Makes single record about meteorite state
    // Returns the proposed time of the next record, e.g.
    //    '0.0f'     - accept any,
    //    't + 0.1f' - skip '0.1' seconds,
    //    'MAX_FLT'  - don't need new information
    // It's a recommendation: you may ignore this value and continue spamming
    virtual real Store(real t, real m, real v, real h, real l, real gamma) = 0;

    // Returns true if the recorder needs information about the trajectory of a meteorite
    // If not, the 'Store()' methods always returns 'false'. So, you can avoid calling it
    virtual bool NeedTrajectory() const = 0;

    // Notifies that computations are ended due to some reason
    // Value 'accuracy' defines how accurately this virtual meteorite describes the real one
    virtual void Finished(Reason reason, double accuracy) = 0;

  protected:
    ISimulationRecorder() = default;
};
