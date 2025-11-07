#pragma once
#include "Meteorites.Core/Defs.h"

#include "IMeteoroidGenerator.h"

// Extends 'IMeteoroidGenerator' by implementing the 'OnProgress' callback
class BasicMeteoroidGenerator : public IMeteoroidGenerator
{
  public:
    // The member of 'IMeteoroidGenerator'
    virtual void OnProgress(const std::function<void(float)> &callback, float step) final override;

  protected:
    BasicMeteoroidGenerator() = default;

    // Notifies that a new meteoroid has been generated
    // So, the 'OnProgress()' callback should be invoked if necessary
    void MovedNext(size_t current, size_t total);

  private:
    std::function<void(float)> callback_;
    size_t n_meteoroids_{};
    float step_{}, threshold_{};
};
