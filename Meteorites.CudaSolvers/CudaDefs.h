#pragma once
#include "Meteorites.Core/Defs.h"

#include <cuda_runtime.h>
#include <curand_kernel.h>
//#include <thrust/device_ptr.h>
//#include <thrust/sort.h>

// Retarget helper functions to device
#if defined(__CUDA_ARCH__)
  #if defined(DEVICE)
    #undef DEVICE
  #endif
  #define DEVICE __device__
#endif


static inline void HandleError(cudaError_t err, const char* file, int line)
{
  if (err != cudaSuccess)
  {
    std::ostringstream oss;
    oss << file << ":" << line << ": " << cudaGetErrorString(err);
    throw std::runtime_error(oss.str());
  }
}
#define HANDLE_ERROR(err) (HandleError(err, __FILE__, __LINE__))


struct CudaDeleter
{
  void operator()(void* ptr) const
  {
    if (ptr != nullptr)
    { cudaFree(ptr); }
  }
};

template <typename T>
using CudaPtr = std::unique_ptr<T, CudaDeleter>;

template <typename T>
cudaError CudaAlloc(CudaPtr<T> &ptr, size_t count)
{
  T *tmp{};
  auto ret = cudaMalloc(&tmp, count * sizeof(T));
  ptr.reset(tmp);
  return ret;
}

// For target parameters of meteoroid. Easy to add new parameters.
// Maybe we should use this throughout the code
class TrajectoryPoint
{
  public:
    real v;
    real h;
};

// Needed to avoid storage and sorting curandState in kernel.
// Seeds for generator can be defined after kernel execution
// but we have to add them to this struct before global sorting
class MeteoroidDeviation
{
  public:
    real     dev;
    // Raw offset for curand_init. DOESN'T displays the number
    // of generated meteoroid
    uint64_t curand_offset; 
};
