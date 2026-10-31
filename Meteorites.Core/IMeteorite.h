#pragma once
#include "Defs.h"

// Interface that describes a single meteorite: its official name, measured parameters, trajectory, etc
// Does not include any assumptions and simulated results
class IMeteorite
{
  public:
    IMeteorite(const IMeteorite &) = delete;
    IMeteorite &operator =(const IMeteorite &) = delete;
    virtual ~IMeteorite() { }

    // Research paper with details about the meteorite
    virtual std::string DOI() const = 0;

    // Official name of the meteorite
    // An empty value means that the meteorite has been detected but not studied
    virtual std::string Name() const = 0;

    // The date of occurrence, local time
    virtual std::string Date() const = 0;

    // Some information about the place of occurrence
    // An empty value means that the meteorite has been detected but not studied
    virtual std::string FallLocation() const = 0;

    // Provides the basic trajectory information that was collected by observers
    virtual void Trajectory(size_t &records, const real *&time,
                            const real *&v, const real *&h) const = 0;

  protected:
    IMeteorite() = default;
};
