#pragma once
#include "Meteorites.Core/Defs.h"
#include "Meteorites.Core/Meteoroids/VirtualMeteoroid.h"

// A factory that generates virtual meteoroids using some rule
// Each of them must be simulated to predict their trajectory and fall point
class IMeteoroidGenerator
{
  public:
    IMeteoroidGenerator(const IMeteoroidGenerator &) = delete;
    IMeteoroidGenerator &operator =(const IMeteoroidGenerator &) = delete;
    virtual ~IMeteoroidGenerator() { }

    // Sets a callback function that will be invoked after processing each portion of meteoroids
    // For instance, 'step = 0.05f' means the callback should be called after every 5% progress
    virtual void OnProgress(const std::function<void(float)> &callback, float step) = 0;

    // Move to the next meteoroid in the sequence and return it if available
    virtual bool MoveNext() = 0;

    // Get the current meteoroid after a call to MoveNext()
    virtual const VirtualMeteoroid &Current() const = 0;

    // Reset the generator to its initial state, so that the next call to MoveNext() starts from the beginning
    virtual void Reset() = 0;

  protected:
    IMeteoroidGenerator() = default;
};
