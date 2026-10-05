#pragma once
#include "Meteorites.Core/Defs.h"

#include "BasicMeteorite.h"

namespace Halliday1981
{

class Innisfree : public BasicMeteorite
{
  public:
    virtual std::string DOI() const override final
    { return "https://doi.org/10.1111/j.1945-5100.1981.tb00540.x"; }

    virtual std::string Name() const override final
    { return "Innisfree"; }

    virtual std::string Date() const override final
    { return "5th of February, 1977"; }

    virtual std::string FallLocation() const override final
    { return "Innisfree, Alberta, Canada"; }

    Innisfree()
    {
      std::vector<real> time = { 0.0, 0.2, 0.4, 0.6, 0.8, 1.0,
                                 1.2, 1.4, 1.6, 1.8, 2.0, 2.2,
                                 2.4, 2.6, 2.8, 3.0, 3.2 };

      std::vector<real> h    = { 58800, 56100, 53500, 50800, 48200, 45500,
                                 42800, 40200, 37500, 35000, 32500, 30200,
                                 27900, 25900, 24200, 22600, 21500 };

      std::vector<real> v    = { 14510, 14490, 14470, 14440, 14340, 14230,
                                 14050, 13790, 13420, 12960, 12350, 11540,
                                 10430, 8890,  7240,  5540,  4700 };
      SetTrajectory(time, v, h);
    }
};

void Populate(std::vector<std::shared_ptr<const IMeteorite>> &records)
{
  records.emplace_back(std::shared_ptr<const IMeteorite>(new Innisfree()));
}

} // namespace Halliday1981
