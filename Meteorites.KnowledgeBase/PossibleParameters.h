#pragma once
#include "Meteorites.Core/Defs.h"

#include "Meteorites.Core/ParameterSet.h"

// Defines a distribution of unknown parameters adapted for some special cases
enum class Distribution : uint32_t
{
  // Provides uniform parameters from ranges that are wide enough
  // May be used to describe any unknown meteorite
  UNIFORM_ANY = 0
};

class PossibleParameters
{
  public:
    PossibleParameters() = delete;
    PossibleParameters(const PossibleParameters &) = delete;
    PossibleParameters &operator =(const PossibleParameters &) = delete;

    static const ParameterSet &Get(Distribution id);

    static const std::vector<ParameterSet> &All();

  private:
    static void LazyInit();

    static std::vector<ParameterSet> collection_;
};
