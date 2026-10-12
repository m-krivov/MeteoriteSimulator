#include "Factory.h"

#include "Meteorites.CpuSolvers/GoldSolver.h"
#include "Meteorites.CpuSolvers/OpenmpSolver.h"

#if defined(METEORITES_CUDA)
  #include "Meteorites.CudaSolvers/PedanticCudaSolver.h"
#endif

std::unique_ptr<ISolver> Factory::StageOneSolver() const
{
#if defined(METEORITES_CUDA)
  if (use_gpu_)
  {
    // Default config should be good
    return std::make_unique<PedanticCudaSolver>();
  }
#endif
  return std::make_unique<OpenmpSolver>();
}

std::unique_ptr<ISolver> Factory::StageTwoSolver() const
{
  return std::make_unique<OpenmpSolver>();
}
