#include "ISimulationRecorder.h"

ISimulationRecorder::Reason ISimulationRecorder::Classify(real t, real M, real h)
{
  if (M <= (real)0.01)
  {
    return Reason::Burnt;
  }
  else if (h <= (real)0.0)
  {
    return Reason::Collided;
  }
  else
  {
    return Reason::Timeouted;
  }
}
