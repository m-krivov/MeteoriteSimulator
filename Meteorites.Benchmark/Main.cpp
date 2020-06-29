#include "PrecisionEstimator.h"

#include "Meteorites.KnowledgeBase/PossibleParameters.h"
#include "Meteorites.KnowledgeBase/KnownMeteorites.h"

#include <iostream>

void main()
{
  size_t points = 0;
  const real *t = nullptr, *v = nullptr, *h = nullptr;
  KnownMeteorites::Get(KnownMeteorites::INNISFREE)->Trajectory(points, t, v, h);
  assert(points > 0);

  auto params = PossibleParameters::Get(PossibleParameters::UNIVERSAL);
  Case problem((params.H().first + params.H().second) / 2,
               (params.Ch().first + params.Ch().second) / 2,
               (params.Rho().first + params.Rho().second) / 2,
               (params.Cd().first + params.Cd().second) / 2,
               (params.Cl().first + params.Cl().second) / 2,
               (params.M0().first + params.M0().second) / 2,
               v[0], h[0],
               (params.Gamma0().first + params.Gamma0().second) / 2);

  PrecisionEstimator precision;
  std::cout << "Mode: FP" << (sizeof(real) == 4 ? "32" : "64") << std::endl << std::endl;
  
  std::cout << "-------------------------" << std::endl;
  std::cout << "--- Different methods ---" << std::endl;
  std::cout << "-------------------------" << std::endl << std::endl;
  std::cout << "dt:             " << (real)1e-3 << std::endl << std::endl;
  precision.CompareMethods(problem, (real)1e-3, std::cout);
  std::cout << std::endl;

  std::cout << "----------------------------" << std::endl;
  std::cout << "--- Different time steps ---" << std::endl;
  std::cout << "----------------------------" << std::endl << std::endl;
  std::cout << "Method's steps: " << GoldSolver::THREE_STEP_ADAMS << std::endl << std::endl;
  precision.CompareSteps(problem, GoldSolver::THREE_STEP_ADAMS, std::cout);

  std::cout << "--------------------" << std::endl;
  std::cout << "--- Peturbations ---" << std::endl;
  std::cout << "--------------------" << std::endl << std::endl;
  std::cout << "dt:             " << (real)1e-3 << std::endl;
  std::cout << "Method's steps: " << GoldSolver::THREE_STEP_ADAMS << std::endl << std::endl;
  precision.ComparePerturbations(problem, GoldSolver::THREE_STEP_ADAMS, (real)1e-3, std::cout);
}
