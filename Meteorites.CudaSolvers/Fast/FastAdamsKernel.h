#pragma once

#include "Meteorites.CudaSolvers/CudaDefs.h"

#include "Meteorites.Core/Meteoroids/VirtualMeteoroid.h"

template <uint32_t STEPS, uint32_t BEST_METEOROIDS_PER_THREAD>
void FastAdamsKernel(const uint64_t *seeds,
                     const real *timestamps, const uint32_t n_timestamps,
                     const TrajectoryPoint *reference_points, TrajectoryPoint *context_points,
                     const real dt, const real timeout,
                     MeteoroidDeviation *global_best_meteoroids, const real border_dev,
                     const uint32_t meteoroids_per_thread,
                     const size_t blocks_num, const size_t threads_per_block);