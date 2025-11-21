#include "BufferingRecorder.h"

//---------------------------
//--- MeteoroidTrajectory ---
//---------------------------

void MeteoroidTrajectory::ExportTo(ISimulationRecorder &recorder) const
{
  recorder.Started(meteoroid_);
  if (recorder.NeedTrajectory())
  {
    for (const auto &record : records_)
    { recorder.Store(record.t, record.M, record.V, record.h, record.l, record.gamma); }
  }
  recorder.Finished(reason_, accuracy_);
}

//-------------------------
//--- BufferingRecorder ---
//-------------------------

real BufferingRecorder::Started(const VirtualMeteoroid &meteoroid)
{
  assert((bool)current_ == false);
  t_next_   = (real)0.0;
  current_  = MeteoroidTrajectory(meteoroid);
  return t_next_;
}

real BufferingRecorder::Store(real t, real m, real v, real h, real l, real gamma)
{
  assert((bool)current_ == true);
  if (t >= t_next_)
  {
    current_->records_.emplace_back(MeteoroidTrajectory::Record(t, m, v, h, l, gamma));
    t_next_ =  std::max(t_next_ + dt_, t);
  }
  return t_next_;
}

void BufferingRecorder::Finished(Reason reason, double accuracy)
{
  assert((bool)current_ == true);
  current_->reason_   = reason;
  current_->accuracy_ = accuracy;
  trajectories_.emplace_back(std::move(current_.value()));
  current_.reset();
}

void BufferingRecorder::MoveTo(std::vector<MeteoroidTrajectory> &trajectories)
{
  trajectories = std::move(trajectories_);
}
