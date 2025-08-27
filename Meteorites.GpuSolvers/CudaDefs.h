#pragma once
#include "Meteorites.Core/Defs.h"

#include <cuda_runtime.h>


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
    if (ptr)
    { cudaFree(ptr); }
  }
};

template <typename T>
using cuda_unique_ptr = std::unique_ptr<T, CudaDeleter>;

template <typename T>
cuda_unique_ptr<T> cuda_make_unique(size_t count)
{
  T *ptr;
  HANDLE_ERROR(cudaMalloc(&ptr, count * sizeof(T)));
  return cuda_unique_ptr<T>(ptr);
}
