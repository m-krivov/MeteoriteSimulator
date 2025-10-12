#include "RestoreCases.h"

__global__ void RestoreCasesKernel(Case* cases, const curandState *states, int N, real v0, real h0)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid < N)
    {
      curandState state = states[tid];
        
      cases[tid] = Case(1e5 + curand_uniform(&state) * (5e6 - 1e5),
                        0.1 + curand_uniform(&state) * (0.9 - 0.1),
                        2000.0 + curand_uniform(&state) * (5000.0 - 2000.0),
                        0.5 + curand_uniform(&state) * (2.5 - 0.5),
                        0.0 + curand_uniform(&state) * (0.25 - 0.0),
                        10.0 + curand_uniform(&state) * (500.0 - 10.0),
                        v0, h0, // maybe add deviations
                        curand_uniform(&state) * (M_PI / 2));
    }
}

void CudaRestoreCasesLauncher(Case *cases, const curandState *states, int N, real v0, real h0)
{
  RestoreCasesKernel<<<((N-1) / 256 + 1), 256>>>(cases, states, N, v0, h0);
}