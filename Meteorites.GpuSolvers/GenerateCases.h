#pragma once

#include "CudaStructures.h"
#include "Meteorites.Core/Case.h"
#include <curand_kernel.h>

void CudaRestoreCasesLauncher(Case *cases, const curandState *states, int N, real v0, real h0);