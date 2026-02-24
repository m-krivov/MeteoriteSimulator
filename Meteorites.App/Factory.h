#pragma once
#include "Meteorites.Core/Defs.h"

#include "Meteorites.Core/Solvers/ISolver.h"

class Factory
{
  public:
#if defined(METEORITES_CUDA)
    Factory(bool use_gpu = false) : use_gpu_(use_gpu) {}
#else
    Factory() : use_gpu_(false) {}
#endif
    Factory(const Factory&) = delete;
    Factory& operator=(const Factory&) = delete;

    std::unique_ptr<ISolver> StageOneSolver() const;

    std::unique_ptr<ISolver> StageTwoSolver() const;

  private:
    bool use_gpu_{};
};
