#pragma once

#include "Meteorites.Core/Meteoroids/BasicMeteoroidGenerator.h"


// Represents a C++ container as a generator of meteoroids
// The 'CONTAINER' type must provide 'begin()' + 'end()' iterators and the 'size()' method
template <typename CONTAINER = std::vector<VirtualMeteoroid>>
class CollectionMeteoroidGenerator : public BasicMeteoroidGenerator
{
  public:
    CollectionMeteoroidGenerator() = delete;
    CollectionMeteoroidGenerator(const CollectionMeteoroidGenerator &) = delete;
    CollectionMeteoroidGenerator &operator =(const CollectionMeteoroidGenerator &) = delete;

    CollectionMeteoroidGenerator(CONTAINER &&meteoroids)
      : meteoroids_(std::move(meteoroids)), current_(meteoroids_.begin())
    {}

    // Version suitable for 'MetaFormatter'
    CollectionMeteoroidGenerator(const std::vector<std::pair<VirtualMeteoroid, double>> &meteoroids)
    {
      meteoroids_.reserve(meteoroids.size());
      for (const auto &[meteoroid, loss] : meteoroids)
      { meteoroids_.emplace_back(meteoroid); }
      current_ = meteoroids_.begin();
    }
        
    // The member of 'IMeteoroidGenerator'
    virtual bool MoveNext() override final
    {
      if (current_ != meteoroids_.end()) {
        ++current_;
        MovedNext(++counter_, meteoroids_.size());
      }
      return current_ != meteoroids_.end();
    }

    // The member of 'IMeteoroidGenerator'
    virtual const VirtualMeteoroid &Current() const override final
    {
      assert(current_ != meteoroids_.end());
      return *current_;
    }

    // The member of 'IMeteoroidGenerator'
    virtual void Reset() override final
    {
      current_ = meteoroids_.begin();
      counter_ = 0;
    }

  private:
    CONTAINER meteoroids_{};
    typename CONTAINER::iterator current_{};
    size_t counter_{};
};
