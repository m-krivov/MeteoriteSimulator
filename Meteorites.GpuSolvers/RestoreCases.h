#pragma once

#include "Meteorites.Core/Case.h"
#include "Meteorites.Core/Defs.h"

void CudaRestoreCasesLauncher(Case *cases, const curandState *states, int N, real v0, real h0);