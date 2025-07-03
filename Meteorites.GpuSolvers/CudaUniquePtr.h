#pragma once

#include "CudaHandleError.h"

#include <memory>
#include <cuda_runtime.h>

struct CudaDeleter {
  void operator()(void* ptr) const {
    if (ptr) {
      cudaFree(ptr);
    }
  }
};

template <typename T>
using cuda_unique_ptr = std::unique_ptr<T, CudaDeleter>;

template <typename T>
cuda_unique_ptr<T> cuda_make_unique(size_t count) {
    T* ptr;
    HANDLE_ERROR(cudaMalloc(&ptr, count * sizeof(T)));
    return cuda_unique_ptr<T>(ptr);
}