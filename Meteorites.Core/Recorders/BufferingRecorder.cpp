#include "BufferingRecorder.h"

BufferingRecorder::BufferingRecorder(real dt)
  : dt_(dt), t_next_((real)0.0)
{
  if (dt < (real)0.0)
  { throw std::runtime_error("'dt' must be declared as positive number or zero"); }
}

real BufferingRecorder::Started(const VirtualMeteoroid &problem)
{
  t_next_   = (real)0.0;
  logs_.emplace_back(Log(problem));
  return t_next_;
}

real BufferingRecorder::Store(real t, real m, real v, real h, real l, real gamma)
{
  if (t >= t_next_)
  {
    assert(!logs_.empty());
    auto &log = logs_[logs_.size() - 1];
    log.records.emplace_back(Record(t, m, v, h, l, gamma));
    t_next_ += dt_;
  }

  return t_next_;
}

void BufferingRecorder::Finished(Reason reason, double accuracy)
{
  assert(!logs_.empty());
  auto &log = logs_[logs_.size() - 1];
  log.reason = reason;
  log.accuracy = accuracy;
}
