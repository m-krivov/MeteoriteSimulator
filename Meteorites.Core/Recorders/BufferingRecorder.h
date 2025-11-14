#pragma once
#include "Meteorites.Core/Defs.h"

#include "ISimulationRecorder.h"

// Stores all data that were passed to this formatter
class BufferingRecorder : public ISimulationRecorder
{
  public:
    // The named analogue of 'std::tuple<real, real, real, real, real>'
    struct Record
    {
      real t{}, M{}, V{}, h{}, l{}, gamma{};

      Record(real t_, real m_, real v_, real h_, real l_, real gamma_)
        : t(t_), M(m_), V(v_), h(h_), l(l_), gamma(gamma_) { }
      Record(const Record &) = default;
      Record &operator =(const Record &) = default;
    };

    // The complete information about single meteorite 
    struct Log
    {
      Log() : reason(Reason::NA), accuracy(std::numeric_limits<double>::max()) { }
      Log(const VirtualMeteoroid &problem_)
        : problem(problem_), reason(Reason::NA), accuracy(std::numeric_limits<double>::max()) { }
      Log(const Log &) = default;
      Log &operator =(const Log &) = default;

      VirtualMeteoroid problem;
      std::vector<Record> records;
      Reason reason;
      double accuracy;
    };

    // Value 'dt' defines how often the incoming records must be stored, in seconds
    BufferingRecorder(real dt);

    // The member of 'ISimulationRecorder'
    virtual real Started(const VirtualMeteoroid &problem) override;

    // The member of 'ISimulationRecorder'
    virtual real Store(real t, real m, real v, real h, real l, real gamma) override;

    // The member of 'ISimulationRecorder'
    virtual bool NeedTrajectory() const override { return true;}
    
    // The member of 'ISimulationRecorder'
    virtual void Finished(Reason reason, double accuracy) override;

    // Information about all trajectories that were passed to this formatter
    const std::vector<Log> &Logs() const { return logs_; }

    // Skips information about all recorded informations
    void Reset() { logs_.clear(); }

  private:
    real dt_, t_next_;
    std::vector<Log> logs_;
};
