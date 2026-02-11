#pragma once
#include "Meteorites.Core/Defs.h"

// Some variation of a linear allocator with expandable preallocated buffer:
//   allocate all -> compute -> free all -> repeat
// You may assume that all allocated memory blocks are aligned
// But please note that this allocator is not thread-safe
class IterationAllocator
{
  private:
    struct State;

  public:
    class Deleter
    {
      public:
        Deleter(State *state, size_t size)
          : state_(state), size_(size)
        { assert(state_ != nullptr); }
        
        void operator()(void *) const
        {
          state_->refs          -= 1;
          state_->cur_allocated -= size_;
          state_->mode           = Mode::Releasing;

          if (state_->refs == 0)
          { delete state_; }
        }
        
      private:
        State *state_{};
        size_t size_{};
    };

    IterationAllocator(size_t initial_size = 1024 * 1024, size_t alignment = 8)
      : alignment_(std::max((size_t)1, alignment)), state_(new State(alignment_))
    {
      state_->refs          = 1;
      state_->max_pool_size = std::max((size_t)1, initial_size);
      state_->pool.reset(AlignedNew<uint8_t>(state_->max_pool_size));
    }
    IterationAllocator(const IterationAllocator &) = delete;
    IterationAllocator &operator =(const IterationAllocator &) = delete;
    ~IterationAllocator()
    {
      state_->refs -= 1;
      if (state_->refs == 0)
      { delete state_; }
    }

    template <typename T = uint8_t>
    std::unique_ptr<T, Deleter> Alloc(size_t count)
    {
      assert(state_ != nullptr);
      count = std::max((size_t)1, count);
      size_t aligned_size = (((count * sizeof(T)) - 1) / alignment_ + 1) * alignment_;

      // All previous blocks must be released before allocating a new one
      // If the pool was too small, enlarge it
      if (state_->mode == Mode::Releasing)
      {
        if (state_->cur_allocated != 0)
        { throw std::bad_alloc(); }
        else
        {
          assert(state_->refs == 1);
          if (state_->max_allocated > state_->max_pool_size)
          {
            state_->pool.reset(AlignedNew<uint8_t>(state_->max_allocated));
            state_->max_pool_size = state_->max_allocated;
          }
          state_->surpluses.clear();
          state_->cur_pool_size = 0;
        }
        state_->mode = Mode::Allocating;
      }

      // If we have enough space, use memory from the pool
      if (state_->cur_pool_size + aligned_size <= state_->max_pool_size)
      {
        uint8_t *ptr            = state_->pool.get() + state_->cur_pool_size;
        state_->cur_pool_size  += aligned_size;
        state_->cur_allocated  += aligned_size;
        state_->max_allocated   = std::max(state_->cur_allocated, state_->max_allocated);
        state_->refs           += 1;
        return std::unique_ptr<T, Deleter>((T *)ptr, Deleter(state_, aligned_size));
      }
      // Otherwise, allocate a temporary block using new[]
      else
      {
        uint8_t *ptr            = nullptr;
        state_->cur_allocated  += aligned_size;
        state_->max_allocated   = std::max(state_->cur_allocated, state_->max_allocated);
        state_->refs           += 1;

        auto block = AlignedUniquePtr<uint8_t[]>(ptr = AlignedNew<uint8_t>(aligned_size),
                                                 AlignedDeleter(alignment_));
        state_->surpluses.emplace_back(std::move(block));
        return std::unique_ptr<T, Deleter>((T *)ptr, Deleter(state_, aligned_size));
      }
    }

  private:
    class AlignedDeleter
    {
      public:
        AlignedDeleter(size_t alignment) : alignment_{ alignment } {}

        void operator()(void *ptr) const
        { ::operator delete[](ptr, alignment_ ); }

      private:
        std::align_val_t alignment_;
    };

    template <typename T>
    T *AlignedNew(size_t n)
    { return (T *)::operator new[](n * sizeof(T), std::align_val_t{ alignment_ }); }

    template <typename T>
    using AlignedUniquePtr = std::unique_ptr<T, AlignedDeleter>;

    enum class Mode
    {
      Allocating,
      Releasing
    };

    struct State
    {
      State(size_t alignment) : pool(nullptr, AlignedDeleter(alignment)) { }

      Mode mode{ Mode::Allocating };
      size_t refs{0};           // one for allocator itself, one for each allocated block
      size_t cur_pool_size{0};  // current offset in the memory pool
      size_t max_pool_size{0};  // size of a block used for the memory pool
      size_t cur_allocated{0};  // total size of all memory blocks allocated by user (pool + surpluses)
      size_t max_allocated{0};  // the peak usage

      AlignedUniquePtr<uint8_t[]> pool;
      std::vector<AlignedUniquePtr<uint8_t[]>> surpluses;
    };

    const size_t alignment_{1};
    State *state_{}; // if allocator is deleted, all allocated blocks will be valid
};
