#include "TestDefs.h"

#include "Meteorites.Core/Functionals/L2Functional.h"
#include "Meteorites.Core/Functionals/BasicFunctional.h"
#include "Meteorites.KnowledgeBase/KnownMeteorites.h"

class WeightedFunctionalTests : public testing::Test
{
};

TEST_F(WeightedFunctionalTests, GenerateDecayingWeights)
{
  // Test uniform weights (decay_factor = 0.0)
  auto uniform_weights = BasicFunctional::GenerateDecayingWeights(5, 0.0);
  ASSERT_EQ(uniform_weights.size(), 5);
  for (const auto &w : uniform_weights)
  {
    ASSERT_FLOAT_EQ(w, 1.0f);
  }

  // Test decaying weights (decay_factor = 0.5)
  auto decaying_weights = BasicFunctional::GenerateDecayingWeights(5, 0.5);
  ASSERT_EQ(decaying_weights.size(), 5);
  
  // First weight should be larger than last
  ASSERT_GT(decaying_weights[0], decaying_weights[4]);
  
  // Weights should be monotonically decreasing
  for (size_t i = 1; i < decaying_weights.size(); i++)
  {
    ASSERT_GE(decaying_weights[i - 1], decaying_weights[i]);
  }
  
  // All weights should be positive
  for (const auto &w : decaying_weights)
  {
    ASSERT_GT(w, 0.0f);
  }
}

TEST_F(WeightedFunctionalTests, WeightedL2Functional)
{
  // Get a meteorite for testing
  decltype(auto) meteorite = KnownMeteorites::Ref().Get(KnownMeteorites::ID::INNISFREE);
  
  // Get the number of records
  size_t n_records = 0;
  const real *time = nullptr, *v = nullptr, *h = nullptr;
  meteorite->Trajectory(n_records, time, v, h);
  ASSERT_GT(n_records, 0);
  
  // Create weighted functional
  auto weights = BasicFunctional::GenerateDecayingWeights(n_records, 0.5);
  L2Functional weighted_func(meteorite, 1.0, 1.0, weights);
  
  // Verify name reflects weighted status
  ASSERT_EQ(weighted_func.Name(), "L2-Weighted");
  
  // Create non-weighted functional for comparison
  L2Functional normal_func(meteorite, 1.0, 1.0);
  ASSERT_EQ(normal_func.Name(), "L2");
  
  // Test that both functionals can compute values
  std::vector<real> test_v(n_records);
  std::vector<real> test_h(n_records);
  for (size_t i = 0; i < n_records; i++)
  {
    test_v[i] = v[i];
    test_h[i] = h[i];
  }
  
  // Perfect match should give low error
  double weighted_error = weighted_func.Compute(n_records, test_v.data(), test_h.data());
  double normal_error = normal_func.Compute(n_records, test_v.data(), test_h.data());
  
  ASSERT_LT(weighted_error, 0.001);
  ASSERT_LT(normal_error, 0.001);
  
  // Test with perturbed values - weighted should handle errors differently
  test_h[n_records - 1] *= 0.9;  // Perturb last height value
  
  double weighted_error_pert = weighted_func.Compute(n_records, test_v.data(), test_h.data());
  double normal_error_pert = normal_func.Compute(n_records, test_v.data(), test_h.data());
  
  ASSERT_GT(weighted_error_pert, weighted_error);
  ASSERT_GT(normal_error_pert, normal_error);
  
  // The weighted functional should penalize the late error less (when base errors are non-zero)
  if (weighted_error > 1e-6 && normal_error > 1e-6)
  {
    ASSERT_LT(weighted_error_pert / weighted_error, normal_error_pert / normal_error);
  }
}

TEST_F(WeightedFunctionalTests, GetStructureDescription)
{
  decltype(auto) meteorite = KnownMeteorites::Ref().Get(KnownMeteorites::ID::INNISFREE);
  
  size_t n_records = 0;
  const real *time = nullptr, *v = nullptr, *h = nullptr;
  meteorite->Trajectory(n_records, time, v, h);
  
  auto weights = BasicFunctional::GenerateDecayingWeights(n_records, 0.5);
  L2Functional weighted_func(meteorite, 1.0, 1.0, weights);
  
  std::string desc = weighted_func.GetStructureDescription();
  
  // Check that description contains expected information
  ASSERT_NE(desc.find("L2-Weighted"), std::string::npos);
  ASSERT_NE(desc.find("Lambda_v"), std::string::npos);
  ASSERT_NE(desc.find("Lambda_h"), std::string::npos);
  ASSERT_NE(desc.find("Weighted: Yes"), std::string::npos);
}
