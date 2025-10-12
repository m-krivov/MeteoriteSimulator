#pragma once
#include "Meteorites.Core/Defs.h"
#include "Meteorites.Core/ParameterSet.h"
#include "Meteorites.Core/IMeteorite.h"
#include "Meteorites.Core/ICaseGenerator.h"

// Generates the required number of random cases
class MonteCarloGenerator : public ICaseGenerator
{
  public:
    MonteCarloGenerator(const IMeteorite &meteorite,
                        const ParameterSet &range,
                        size_t n_cases, uint32_t seed);

    virtual void OnProgress(const std::function<void(float)> &callback, float step) override;

    virtual bool Next(Case &problem) override;

    //virtual size_t N_cases() {return n_cases_;}

  private:
    mutable std::mt19937 gen_;
    std::uniform_real_distribution<real> dist_;

    ParameterSet range_;
    real v0_, h0_;
    size_t cur_, n_cases_;

    std::function<void(float)> callback_;
    float step_, threshold_;
};
