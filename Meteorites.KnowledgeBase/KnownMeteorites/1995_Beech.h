#pragma once
#include "Meteorites.Core/Defs.h"

#include "BasicMeteorite.h"

namespace Beech1995
{

class Peekskill  : public BasicMeteorite
{
  public:
    virtual std::string DOI() const override final
    { return "https://doi.org/10.1007/BF00671508"; }

    virtual std::string Name() const override final
    { return "Peekskill"; }

    virtual std::string Date() const override final
    { return "9th of October, 1992"; }

    virtual std::string FallLocation() const override final
    { return "Peekskill, New York, USA"; }

    Peekskill()
    {
      std::vector<real> time = { 25.0,
                                 26.0, 27.0, 28.0, 29.0, 30.0, 
                                 31.0, 32.0, 33.0, 34.0, 35.0,
                                 36.0, 37.0, 38.0, 39.0, 40.0,
                                 41.0, 42.0, 43.0, 44.0, 45.0,
                                 46.0, 47.0, 48.0, 49.0, 50.0,
                                 51.0, 52.0, 53.0, 54.0, 55.0,
                                 56.0, 57.0, 58.0, 59.0, 60.0 };
      std::vector<real> h    = { 48.9,
                                 47.4, 46.3, 45.2, 44.2, 43.0,
                                 42.3, 41.4, 40.6, 39.8, 39.0,
                                 38.1, 37.6, 36.9, 36.2, 35.9,
                                 35.3, 34.9, 34.3, 33.9, 33.6,
                                 33.3, 33.0, 32.8, 32.4, 32.4,
                                 32.1, 31.9, 31.7, 31.5, 31.3,
                                 31.3, 31.1, 30.9, 30.9, 30.9 };
      std::vector<real> v    = { 14.9,
                                 14.5, 14.6, 14.6, 14.6, 14.3,
                                 14.3, 13.9, 14.0, 13.7, 13.7,
                                 13.1, 12.7, 12.4, 11.8, 11.5,
                                 10.7, 10.4,  9.8,  9.6,  8.9,
                                  8.5,  8.3,  7.8,  7.2,  7.1,
                                  6.7,  6.3,  6.4,  6.1,  5.8,
                                  5.5,  5.1,  4.5,  3.9,  3.4 };

      for (auto &r : h) { r *= 1000; }
      for (auto &r : v) { r *= 1000; }
      SetTrajectory(time, v, h);
    }
};

void Populate(std::vector<std::shared_ptr<const IMeteorite>> &records)
{
  records.emplace_back(std::shared_ptr<const IMeteorite>(new Peekskill()));
}

} // namespace Beech1995
