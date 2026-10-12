#pragma once
#include "Meteorites.Core/Defs.h"

#include "MeteoroidSummary.h"


// The complete information about the trajectory of a single virtual meteoroid
class MeteoroidTrajectory : public MeteoroidSummary
{
  public:
    // Just a named analogue of 'std::tuple<real, real, real, real, real>'
    struct Record
    {
      real t{}, M{}, V{}, h{}, l{}, gamma{};

      Record(real t_, real m_, real v_, real h_, real l_, real gamma_) noexcept
        : t(t_), M(m_), V(v_), h(h_), l(l_), gamma(gamma_) { }
      Record(const Record &) noexcept = default;
      Record &operator =(const Record &) noexcept = default;
    };

    MeteoroidTrajectory() noexcept = default;
    MeteoroidTrajectory(MeteoroidTrajectory &&) noexcept = default;
    MeteoroidTrajectory(const MeteoroidTrajectory &) noexcept = default;

    MeteoroidTrajectory &operator =(MeteoroidTrajectory &&m) noexcept
    {
      if (&m != this)
      {
        MeteoroidSummary::operator=(m);
        records_ = std::move(m.records_);
      }
      return *this;
    }

    MeteoroidTrajectory &operator =(const MeteoroidTrajectory &m) noexcept
    {
      if (&m != this)
      {
        MeteoroidSummary::operator=(m);
        records_ = m.records_;
      }
      return *this;
    }

    // Returns all recorded trajectory points sorted by time
    const std::vector<Record> &Records() const { return records_; }

    // Returns the last record of the trajectory
    const Record &LastRecord() const { assert(!records_.empty()); return records_.back(); }

    // Exports all trajectory data to another instance of 'ISimulationRecorder'
    void ExportTo(ISimulationRecorder &recorder) const
    {
      recorder.Started(Meteoroid());
      if (recorder.NeedTrajectory())
      {
        for (const auto &record : records_)
        { recorder.Store(record.t, record.M, record.V, record.h, record.l, record.gamma); }
      }
      recorder.Finished(Reason(), Accuracy());
    }
  
  private:
    explicit MeteoroidTrajectory(std::tuple<VirtualMeteoroid,
                                            ISimulationRecorder::Reason,
                                            double,
                                            std::vector<Record>> &data)
      : MeteoroidSummary(std::get<0>(data), std::get<1>(data), std::get<2>(data)),
        records_(std::move(std::get<3>(data)))
    { }

    std::vector<Record> records_;

  friend class BufferingRecorder;
};
