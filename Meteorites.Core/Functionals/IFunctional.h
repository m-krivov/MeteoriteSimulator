#pragma once
#include "Meteorites.Core/Defs.h"

// Interface that defines a functional for measuring distance between virtual and real meteorites
class IFunctional
{
  public:
    IFunctional(const IFunctional &) = delete;
    IFunctional &operator =(const IFunctional &) = delete;
    virtual ~IFunctional() = default;

    // Returns the name of the functional (not unique!)
    // E.g. 'L2' or 'C'
    virtual std::string Name() const = 0;

    // Retrieves the timestamps used as input arguments
    // At least one timestamp is guaranteed to be provided
    virtual void GetTimeStamps(size_t &num, const real *&values) const = 0;

    // Computes the functional's value for the given arguments
    // The actual number of arguments may differ from expected: such values just increases error
    virtual double Compute(size_t num, const real *v, const real *h) const = 0;

  protected:
    IFunctional() = default;
};
