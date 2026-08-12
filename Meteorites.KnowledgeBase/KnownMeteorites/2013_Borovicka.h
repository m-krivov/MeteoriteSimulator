#pragma once
#include "Meteorites.Core/Defs.h"

#include "BasicMeteorite.h"

namespace Borovicka2013
{

class Chelyabinsk : public BasicMeteorite
{
  public:
    virtual std::string DOI() const override final
    { return "https://doi.org/10.1038/nature12671"; }

    virtual std::string Name() const override final
    { return "Chelyabinsk"; }

    virtual std::string Date() const override final
    { return "15th of February, 2013"; }

    virtual std::string FallLocation() const override final
    { return "Chebarkul, Chelyabinskaya Oblast, Russia"; }

    Chelyabinsk()
    {
      std::vector<real> time = { 1.07,  6.97,  10.46, 12.24, 13.18, 14.18, 15.17 };
      std::vector<real> v    = { 19.03, 19.05, 19.03, 18.9,  18.0,  14.2,  6 };
      std::vector<real> h    = { 95.0,  60.0,  40.0,  30.0,  25.0,  20.0,  17.2 };

      for (auto &r : v) { r *= 1000; }
      for (auto &r : h) { r *= 1000;}
      SetTrajectory(time, v, h);
    }
};

void Populate(std::vector<std::shared_ptr<const IMeteorite>> &records)
{
  records.emplace_back(std::shared_ptr<const IMeteorite>(new Chelyabinsk()));
}

} // namespace Borovicka2013
