#pragma once
#include "Meteorites.Core/Defs.h"

#include "ISimulationRecorder.h"


// The complete information about the trajectory of a single virtual meteoroid
class MeteoroidTrajectory
{
  public:
    // Just a named analogue of 'std::tuple<real, real, real, real, real>'
    struct Record
    {
      real t{}, M{}, V{}, h{}, l{}, gamma{};

      Record(real t_, real m_, real v_, real h_, real l_, real gamma_)
        : t(t_), M(m_), V(v_), h(h_), l(l_), gamma(gamma_) { }
      Record(const Record &) = default;
      Record &operator =(const Record &) = default;
    };

    MeteoroidTrajectory() = delete;
    MeteoroidTrajectory(const MeteoroidTrajectory &) noexcept = default;
    MeteoroidTrajectory(MeteoroidTrajectory &&) noexcept = default;
    MeteoroidTrajectory &operator =(const MeteoroidTrajectory &) noexcept = default;
    MeteoroidTrajectory &operator =(MeteoroidTrajectory &&) noexcept = default;

    // Returns the original virtual meteoroid used for the simulation
    const VirtualMeteoroid &Meteoroid() const { return meteoroid_; }

    // Returns all recorded trajectory points sorted by time
    const std::vector<Record> &Records() const { return records_; }

    // Returns the last record of the trajectory
    const Record &LastRecord() const { assert(!records_.empty()); return records_[records_.size() - 1]; }

    // Returns the reason why the simulation was stopped
    ISimulationRecorder::Reason Reason() const { return reason_; }

    // Returns the value of the used functional
    // In the other words, accuracy of the simulation compared to the real meteorite
    double Accuracy() const { return accuracy_; }

    // Export all trajectory data to another instance of 'ISimulationRecorder'
    void ExportTo(ISimulationRecorder &recorder) const;
  
  private:
    MeteoroidTrajectory(const VirtualMeteoroid &meteoroid) : meteoroid_(meteoroid) { }

    VirtualMeteoroid meteoroid_;
    std::vector<Record> records_;
    ISimulationRecorder::Reason reason_{ISimulationRecorder::Reason::NA};
    double accuracy_{std::numeric_limits<double>::max()};

  friend class BufferingRecorder;
};


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
    // Equivalent to the following code:
    //   trajectories = recorder.Trajectories();
    //   recorder.Reset();
    void MoveTo(std::vector<MeteoroidTrajectory> &trajectories);

  private:
    real dt_{}, t_next_{};
    std::optional<MeteoroidTrajectory> current_;
    std::vector<MeteoroidTrajectory> trajectories_;
};
