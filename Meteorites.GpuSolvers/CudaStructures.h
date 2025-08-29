#pragma once
#include "Meteorites.GpuSolvers/CudaDefs.h"

#include "Meteorites.Core/Adams.h"
#include "Meteorites.Core/Case.h"
#include "Meteorites.Core/Constants.h"

class Record
{
  public:
  real t_, m_, v_, h_, l_, gamma_;

  /*test*/ void print() { printf("rec: %f %f %f %f %f %f\n", t_, m_, v_, h_, l_, gamma_); }
};

template <unsigned int STEPS> class ThreadContext
{
  public:
  Adams::Layer steps[STEPS + 1];
  Case params;
  size_t nxt;
  real t;
  size_t timestamp;
  size_t curr_case_num;
  bool ended;
};
