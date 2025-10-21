#pragma once
#include "Meteorites.Core/Defs.h"

#include <cuda_runtime.h>

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
