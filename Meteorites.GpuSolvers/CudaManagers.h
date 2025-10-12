#pragma once

#include "Meteorites.Core/IMeteorite.h"
#include "Meteorites.Core/Case.h"
#include "Meteorites.Core/Defs.h"

template <uint32_t STEPS>
void CudaAdamsMethodManager(const IMeteorite &meteorite, uint64_t *seeds, int N,
                 real dt, real timeout, real border_dev,
                 uint32_t BLOCKS_NUM, uint32_t THREADS_PER_BLOCK,
                 uint32_t BEST_CASES_BUFFER_SIZE, uint32_t CASES_PER_THREAD,
                 curandState *result_curand_states, real *result_deviations);

std::vector<Case> CudaRestoreCasesManager(const curandState *states, int N, real v0, real h0);