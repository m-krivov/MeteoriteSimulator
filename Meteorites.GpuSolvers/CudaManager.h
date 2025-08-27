#pragma once
#include "Meteorites.GpuSolvers/CudaDefs.h"

#include "Meteorites.Core/Case.h"
#include "Meteorites.Core/Functionals.h"
#include "Meteorites.Core/ResultFormatters.h"

template <unsigned int STEPS, unsigned int ITERS>
void CudaManager(const std::vector<Case>& problems_vector, const IFunctional& functional, real dt, real timeout,
                 IResultFormatter& results);
