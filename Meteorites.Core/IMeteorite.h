#pragma once
#include "Defs.h"

// Interface that describes a single meteorite: its official name, measured parameters, trajectory, etc
// Does not include any assumptions and simulated results
class IMeteorite
{
  public:
    IMeteorite(const IMeteorite &) = delete;
    void operator =(const IMeteorite &) = delete;
    virtual ~IMeteorite() { }

    // Official name of the meteorite
    virtual std::string Name() const = 0;

    // The date of occurrence, local time
    virtual std::string Date() const = 0;

    // Some information about the place of occurrence
    virtual std::string FallLocation() const = 0;

    // Provides the basic trajectory information that was collected by observers
    virtual void Trajectory(size_t &records, const real *&time,
                            const real *&v, const real *&h) const = 0;

  protected:
    IMeteorite() = default;
};
