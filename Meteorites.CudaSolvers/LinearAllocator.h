#pragma once
#include "Meteorites.Core/Defs.h"

// Linear allocator with expandable preallocated buffer. 
// Provides UniquePtr for T objects.
// Allocates memory blocks aligned with ALIGNMENT.
template <typename T = uint8_t, size_t ALIGNMENT = 1>
class LinearAllocator
{
  private:
    class Deleter
    {
      public:
        Deleter(LinearAllocator *allocator) : allocator_(allocator) {}
        
        void operator()(T *ptr) const
        {
          if (allocator_) allocator_->Free();
        }
        
      private:
        LinearAllocator *allocator_;
    };

  public:
    using UniquePtr = std::unique_ptr<T, Deleter>;

    LinearAllocator() : buffers_({nullptr}), ptrs_count_(0), allocated_bytes_(0), preallocated_bytes_(0) {}
    ~LinearAllocator() { if (buffers_[0]) free(buffers_[0]);}

    UniquePtr Alloc(size_t count)
    {
      size_t aligned_size = (((count * sizeof(T)) - 1) / ALIGNMENT + 1) * ALIGNMENT;
      if (preallocated_bytes_ >= aligned_size)
      {
        ptrs_count_++;
        preallocated_bytes_ -= aligned_size;
        allocated_bytes_ += aligned_size;
        return UniquePtr(reinterpret_cast<T*>(buffers_[0] + preallocated_bytes_), Deleter(this));
      }
      ptrs_count_++;
      allocated_bytes_ += aligned_size;
      buffers_.emplace_back(static_cast<uint8_t*>(malloc(aligned_size)));
      if (!buffers_.back()) throw std::bad_alloc();
      return UniquePtr(reinterpret_cast<T*>(buffers_.back()), Deleter(this));
    }

    void Free()
    {
      ptrs_count_--;
      if (!ptrs_count_)
      {
        for (size_t i = 1; i < buffers_.size(); i++)
        {
          free(buffers_[i]);
        }
        buffers_.resize(1);
        buffers_[0] = static_cast<uint8_t*>(realloc(buffers_[0], allocated_bytes_));
        if (!buffers_[0]) throw std::bad_alloc();
        preallocated_bytes_ = allocated_bytes_;
        allocated_bytes_ = 0;
      }
    }

  private:
    std::vector<uint8_t*> buffers_;
    size_t ptrs_count_;
    size_t allocated_bytes_;
    size_t preallocated_bytes_;
};