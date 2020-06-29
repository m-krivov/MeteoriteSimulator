#pragma once
#include "Meteorites.Core/Defs.h"

#include "Meteorites.Core/ParameterSet.h"

class PossibleParameters
{
  public:
    enum SetID
    {
      // Provides parameters from ranges that are wide enough
      // May be used to describe any unknown meteorite
      UNIVERSAL = 0
    };

    PossibleParameters() = delete;
    PossibleParameters(const PossibleParameters &) = delete;
    PossibleParameters &operator =(const PossibleParameters &) = delete;

    static const ParameterSet &Get(SetID id);

    static const std::vector<ParameterSet> &All();

  private:
    static void LazyInit();

    static std::vector<ParameterSet> collection_;
};
