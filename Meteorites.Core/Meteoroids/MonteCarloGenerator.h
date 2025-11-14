#pragma once
#include "Meteorites.Core/Defs.h"

#include "Meteorites.Core/ParameterSet.h"
#include "Meteorites.Core/IMeteorite.h"
#include "BasicMeteoroidGenerator.h"

// Generates the required number of random meteoroids using the Monte-Carlo method
class MonteCarloGenerator : public BasicMeteoroidGenerator
{
  public:
    MonteCarloGenerator(const IMeteorite &meteorite,
                        const ParameterSet &range,
                        size_t n_cases, uint64_t seed);

    // The member of 'IMeteoroidGenerator'
    virtual bool MoveNext() override;

    // The member of 'IMeteoroidGenerator'
    virtual const VirtualMeteoroid &Current() const override { return current_; }

    // The member of 'IMeteoroidGenerator'
    virtual void Reset() override;

  private:
    uint64_t seed_{};
    mutable std::mt19937 gen_{};
    std::uniform_real_distribution<real> dist_{};

    VirtualMeteoroid current_{};

    ParameterSet range_;
    real v0_{}, h0_{};
    size_t case_{}, n_cases_{};
};
