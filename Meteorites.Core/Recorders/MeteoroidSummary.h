#pragma once
#include "Meteorites.Core/Defs.h"

#include "ISimulationRecorder.h"
#include "Meteorites.Core/Meteoroids/VirtualMeteoroid.h"


// Contains only metadata about a completed simulation
// Its trajectory is skipped, only the final state is presented
class MeteoroidSummary
{
  public:
    MeteoroidSummary() = default;
    MeteoroidSummary(MeteoroidSummary &&) noexcept = default;
    MeteoroidSummary(const MeteoroidSummary &) noexcept = default;
    MeteoroidSummary &operator =(MeteoroidSummary &&) noexcept = default;
    MeteoroidSummary &operator =(const MeteoroidSummary &) noexcept = default;

    // Returns the original virtual meteoroid used for the simulation
    const VirtualMeteoroid &Meteoroid() const { return meteoroid_; }

    // Returns the reason why the simulation was stopped
    ISimulationRecorder::Reason Reason() const { return reason_; }

    // Returns the value of the used functional
    // In the other words, accuracy of the simulation compared to the real meteorite
    double Accuracy() const { return accuracy_; }

    // Exports meta information to some instance of 'ISimulationRecorder'
    // Warning: 'recorder->NeedTrajectory()' must be equal to false!
    void ExportTo(ISimulationRecorder &recorder) const
    {
      assert(!recorder.NeedTrajectory());
      recorder.Started(meteoroid_);
      recorder.Finished(reason_, accuracy_);
    }

  protected:
    MeteoroidSummary(const VirtualMeteoroid &meteoroid,
                     ISimulationRecorder::Reason reason,
                     double accuracy)
      : meteoroid_(meteoroid), reason_(reason), accuracy_(accuracy)
    { assert(accuracy >= 0.0); }

  private:
    VirtualMeteoroid meteoroid_{};
    ISimulationRecorder::Reason reason_{ ISimulationRecorder::Reason::NA };
    double accuracy_{ std::numeric_limits<double>::max() };
  
  friend class MetaRecorder;
};
