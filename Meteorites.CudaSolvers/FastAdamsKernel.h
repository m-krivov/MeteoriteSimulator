#pragma once

#include "Meteorites.CudaSolvers/CudaDefs.h"

#include "Meteorites.Core/Meteoroids/VirtualMeteoroid.h"

void RestoreMeteoroidsKernel(VirtualMeteoroid *meteoroids, const curandState *states,
                             const size_t size, const real v0, const real h0);

template <unsigned int STEPS, size_t BEST_METEOROIDS_BUFFER_SIZE>
void FastAdamsKernel(const uint64_t *seeds, const size_t n_problems,
                     const real dt, const real timeout,
                     const real *timestamps, size_t n_timestamps,
                     real *functional_args, const real *functional_args_references,
                     real *best_deviations, curandState *best_curand_states, const real border_deviation,
                     const size_t meteoroids_per_thread,
                     const size_t blocks_num, const size_t threads_per_block);