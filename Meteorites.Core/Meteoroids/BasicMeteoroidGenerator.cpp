#include "BasicMeteoroidGenerator.h"

void BasicMeteoroidGenerator::OnProgress(const std::function<void(float)> &callback, float step)
{
  assert(step <= 1.0f);
  callback_  = callback;
  step_      = step;
  threshold_ = step;
}

void BasicMeteoroidGenerator::MovedNext(size_t current, size_t total)
{
  assert(current <= total);
  if (!callback_)
  { return; }

  // Generate the required number of events
  auto ratio = (float)current / total;
  while (ratio >= threshold_)
  {
    callback_(threshold_);
    threshold_ += step_;
  }

  // Process the last portion without taking into account the 'step_' increment
  if (ratio < threshold_ && current == total)
  { callback_(ratio); }
}
