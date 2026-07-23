#pragma once
#include "Meteorites.Core/Defs.h"

#include "BasicMeteorite.h"

namespace Halliday1996
{

class Innisfree : public BasicMeteorite
{
  public:
    virtual std::string DOI() const override final
    { return "https://doi.org/10.1111/j.1945-5100.1996.tb02014.x"; }

    virtual std::string Name() const override final
    { return "Innisfree"; }

    virtual std::string Date() const override final
    { return "6th of February, 1977"; }

    virtual std::string FallLocation() const override final
    { return "Innisfree, Alberta, Canada"; }

    Innisfree()
    {
      // Object No 285
      std::vector<real> time = { 0.00, 0.60, 1.20,
                                 1.80, 2.20, 2.60,
                                 3.00, 3.32, 3.82 };

      std::vector<real> h    = { 58.8, 50.8, 42.8,
                                 35.0, 30.2, 25.9,
                                 22.6, 21.7, 19.8 };

      std::vector<real> v    = { 14.5, 14.4, 14.2,
                                 13.4, 12.4, 10.4,
                                 7.2,  5.3,  2.7 };
      for (auto &r : h) { r *= 1000; }
      for (auto &r : v) { r *= 1000; }
      SetTrajectory(time, v, h);
    }
};

void Populate(std::vector<std::shared_ptr<const IMeteorite>> &records)
{
  records.emplace_back(std::shared_ptr<const IMeteorite>(new Innisfree()));
}

} // namespace Halliday1996
