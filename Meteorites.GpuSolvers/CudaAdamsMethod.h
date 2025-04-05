#pragma once

#include "Meteorites.Core/Case.h"
#include "Meteorites.Core/Functionals.h"
#include "Meteorites.Core/ResultFormatters.h"

template <unsigned int STEPS, unsigned int ITERS>
void CudaAdamsMethod(const std::vector<Case> &problem_vector, const IFunctional &functional, real dt, real timeout,
                     IResultFormatter &results);
