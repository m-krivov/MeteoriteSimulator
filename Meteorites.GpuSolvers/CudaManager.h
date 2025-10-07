#pragma once
#include "Meteorites.GpuSolvers/CudaDefs.h"

#include "Meteorites.Core/Case.h"
#include "Meteorites.Core/Functionals.h"
#include "Meteorites.Core/ResultFormatters.h"

template <uint32_t STEPS>
void CudaManager(const std::vector<Case> &problems, const IFunctional &functional,
                 real dt, real timeout, IResultFormatter &results);
