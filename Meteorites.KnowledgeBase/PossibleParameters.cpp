#include "PossibleParameters.h"

//--------------------------
//--- PossibleParameters ---
//--------------------------

std::vector<ParameterSet> PossibleParameters::collection_;

const ParameterSet &PossibleParameters::Get(Distribution id)
{
  LazyInit();
  if ((uint32_t)id >= collection_.size())
  { throw std::runtime_error("wrong type of parameter set"); }

  return collection_[(uint32_t)id];
}

const std::vector<ParameterSet> &PossibleParameters::All()
{
  LazyInit();
  return collection_;
}

void PossibleParameters::LazyInit()
{
  if (collection_.empty())
  {
    ParameterSet universal(
      std::make_pair((real)1e5, (real)5e6),           // H
      std::make_pair((real)0.1, (real)0.9),           // Ch
      std::make_pair((real)2000.0, (real)5000.0),     // Rho
      std::make_pair((real)0.5, (real)2.5),           // Cd
      std::make_pair((real)0.0, (real)0.25),          // Cl
      std::make_pair((real)10.0, (real)500.0),        // M0
      std::make_pair((real)0.0, (real)std::_Pi / 2)   // Gamma0
    );
    collection_.emplace_back(universal);
    assert(collection_.size() == (uint32_t)Distribution::UNIFORM_ANY + 1);
  }
}
