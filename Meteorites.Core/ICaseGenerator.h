#pragma once
#include "Defs.h"
#include "Case.h"

// Generates cases for virtual meteorites using some rule
class ICaseGenerator
{
  public:
    ICaseGenerator(const ICaseGenerator &) = delete;
    ICaseGenerator &operator =(const ICaseGenerator &) = delete;
    virtual ~ICaseGenerator() { }

    // Sets callback that will be called after each portion of the processed cases
    // For example, 'step = 0.05f' means that the callback must be called after each 5%
    virtual void OnProgress(const std::function<void(float)> &callback, float step) = 0;

    // Spawns new problem that uniquely describes some virtual meteorite
    // Returns true if some other cases must be considered as well
    virtual bool Next(Case &problem) = 0;

  protected:
    ICaseGenerator() = default;
};
