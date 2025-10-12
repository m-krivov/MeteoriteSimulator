#pragma once

#include "Meteorites.Core/Defs.h"

constexpr int BEST_CASES_BUFFER_SIZE = 4;

template <unsigned STEPS>
void CudaAdamsMethodLauncher(const uint64_t *seeds,
                             const real *timestamps, size_t n_timestamps,
                             const real *ref_v, const real *ref_h,
                             real dt, real timeout,
                             real *best_cases_devs, curandState *best_cases_states, real border_dev,
                             real *v_args_arr, real *h_args_arr,
                             uint32_t CASES_PER_THREAD,
                             uint32_t BLOCKS_NUM, uint32_t THREADS_PER_BLOCK);
