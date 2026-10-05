#pragma once
#include "Meteorites.Core/Defs.h"
#include "Meteorites.Core/IMeteorite.h"


// Basic class for 'IMeteorites' that converts height and velocity into the desired format
class BasicMeteorite : public IMeteorite
{
  public:
    virtual void Trajectory(size_t &records, const real *&time,
                            const real *&v, const real *&h) const override final
    {
      assert(!time_.empty());

      records = time_.size();
      time = time_.data();
      v = v_.data();
      h = h_.data();
    }

  protected:
    BasicMeteorite() = default;

    void SetTrajectory(const std::vector<real> &time,
                       const std::vector<real> &v,
                       const std::vector<real> &h)
    {
      assert(!time.empty());
      assert(time.size() == v.size());
      assert(time.size() == h.size());

      time_.reserve(time.size());
      for (auto r : time)
      { time_.emplace_back((real)(r - time[0])); }

      v_.reserve(v.size());
      for (auto r : v)
      { v_.emplace_back((real)r); }

      h_.reserve(h.size());
      for (auto r : h)
      { h_.emplace_back((real)r); }
    }

  private:
    std::vector<real> time_, v_, h_;
};
