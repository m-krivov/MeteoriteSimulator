#include "BufferingRecorder.h"


//-------------------------
//--- BufferingRecorder ---
//-------------------------

real BufferingRecorder::Started(const VirtualMeteoroid &meteoroid)
{
  assert((bool)current_ == false);
  t_next_   = (real)0.0;
  current_.emplace(std::make_tuple(meteoroid, ISimulationRecorder::Reason::NA,
                                   0.0, std::vector<MeteoroidTrajectory::Record>()));
  return t_next_;
}

real BufferingRecorder::Store(real t, real m, real v, real h, real l, real gamma)
{
  assert((bool)current_ == true);
  if (t >= t_next_)
  {
    std::get<3>(current_.value()).emplace_back(MeteoroidTrajectory::Record(t, m, v, h, l, gamma));
    t_next_ =  std::max(t_next_ + dt_, t);
  }
  return t_next_;
}

void BufferingRecorder::Finished(Reason reason, double accuracy)
{
  assert((bool)current_ == true);
  std::get<1>(current_.value())  = reason;
  std::get<2>(current_.value())  = accuracy;
  trajectories_.emplace_back(MeteoroidTrajectory(current_.value()));
  current_.reset();
}

void BufferingRecorder::MoveTo(std::vector<MeteoroidTrajectory> &trajectories)
{
  assert((bool)current_ == false);
  trajectories = std::move(trajectories_);
}
