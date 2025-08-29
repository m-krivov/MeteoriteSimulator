#pragma once
#include "Meteorites.GpuSolvers/CudaDefs.h"

#include "CudaStructures.h"

template <unsigned int STEPS, unsigned int ITERS>
void CudaLauncher(ThreadContext<STEPS>* dev_thread_sandbox_arr, Case* dev_problems,
                  size_t promblems_vector_size, real dt, real timeout, real* dev_timestamps, size_t n_timestamps,
                  real* dev_functional_args, Record* dev_records, unsigned int* dev_active_threads,
                  unsigned int launching_blocks, unsigned int launching_threads_per_block);
