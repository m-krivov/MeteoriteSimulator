#pragma once

#include "Meteorites.Core/Meteoroids/BasicMeteoroidGenerator.h"


// Represents a C++ container as a generator of meteoroids
// The 'CONTAINER' type must provide 'begin()' + 'end()' iterators and the 'size()' method
template <typename CONTAINER = std::vector<VirtualMeteoroid>>
class CollectionGenerator : public BasicMeteoroidGenerator
{
  public:
    CollectionGenerator() = delete;
    CollectionGenerator(const CollectionGenerator &) = delete;
    CollectionGenerator &operator =(const CollectionGenerator &) = delete;

    CollectionGenerator(CONTAINER &&meteoroids)
      : meteoroids_(std::move(meteoroids)), current_(meteoroids_.begin())
    {}

    // Version suitable for 'MetaFormatter'
    CollectionGenerator(const std::vector<std::pair<VirtualMeteoroid, double>> &meteoroids)
    {
      meteoroids_.reserve(meteoroids.size());
      for (const auto &[meteoroid, loss] : meteoroids)
      { meteoroids_.emplace_back(meteoroid); }
      current_ = meteoroids_.begin();
    }
        
    // The member of 'IMeteoroidGenerator'
    virtual bool MoveNext() override final
    {
      if (current_ != meteoroids_.end())
      {
        if (!started_)
        { started_ = true; }
        else
        { ++current_; }
      }
      
      if (current_ != meteoroids_.end())
      {
        MovedNext(++counter_, meteoroids_.size());
        return true;
      }
      else
      { return false; }
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
      started_ = false;
      current_ = meteoroids_.begin();
      counter_ = 0;
    }

  private:
    CONTAINER meteoroids_{};
    bool started_{false};
    typename CONTAINER::iterator current_{};
    size_t counter_{};
};
