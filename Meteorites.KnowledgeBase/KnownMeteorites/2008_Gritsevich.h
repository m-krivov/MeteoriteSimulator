#pragma once
#include "Meteorites.Core/Defs.h"

#include "BasicMeteorite.h"

namespace Gritsevich2008
{

class Pribram : public BasicMeteorite
{
  public:
    virtual std::string DOI() const override final
    { return "https://doi.org/10.1134/S003809460805002X"; }

    virtual std::string Name() const override final
    { return "Pribram"; }

    virtual std::string Date() const override final
    {return "7th of April, 1959"; }

    virtual std::string FallLocation() const override final
    { return "Pribram, Czechoslovakia (Czech Republic)"; }

    Pribram()
    {
      std::vector<real> time = { 0,      0.856,  1.732,
                                 2.494,  2.692,  3.0676 };

      std::vector<real> h    = { 88.594, 76.318, 63.837,
                                 52.97,  50.164, 44.858 };

      std::vector<real> v    = { 20.887, 20.86,  20.838,
                                 20.773, 20.717, 20.459 };

      for (auto &r : h) { r *= 1000; }
      for (auto &r : v) { r *= 1000; }
      SetTrajectory(time, v, h);
    }
};


class LostCity : public BasicMeteorite
{
  public:
    virtual std::string DOI() const override final
    { return "https://doi.org/10.1134/S003809460805002X"; }

    virtual std::string Name() const override final
    { return "Lost City"; }

    virtual std::string Date() const override final
    { return "3rd of January, 1970"; }

    virtual std::string FallLocation() const override final
    { return "Lost City, Oklahoma, USA"; }

    LostCity()
    {
      std::vector<real> time = { 0.05, 1.05, 2.05, 3.05, 4.05,
                                 5.05, 6.05, 7.05, 8.00, 8.95 };

      std::vector<real> h    = { 85.9, 77.1, 68.5, 59.9, 51.3,
                                 42.8, 34.6, 27.5, 22.6, 19.9 };

      std::vector<real> v    = { 14.2, 14.2, 14.2, 14.1, 14.0,
                                 13.8, 12.9, 10.3, 6.1,  3.4 };

      for (auto &r : h) { r *= 1000; }
      for (auto &r : v) { r *= 1000; }
      SetTrajectory(time, v, h);
    }
};

void Populate(std::vector<std::shared_ptr<const IMeteorite>> &records)
{
  records.emplace_back(std::shared_ptr<const IMeteorite>(new Pribram()));
  records.emplace_back(std::shared_ptr<const IMeteorite>(new LostCity()));
}

} // namespace Gritsevich2008
