# Accuracy and Stability Improvements

This document describes the improvements made to enhance the accuracy and stability of the meteorite simulator, particularly for challenging "bad" meteorites where trajectory data quality varies.

## Overview

Four major improvements have been implemented:

1. **Weighted Functionals** - Handle varying measurement accuracy
2. **Functional Structure Visualization** - Document functional configuration
3. **Burn Time Penalty** - Match observation duration
4. **Stage 0 Entry Angle Optimization** - Reduce computational requirements

## 1. Weighted Functionals

### Problem
Real meteorite observations often have decreasing accuracy over time. Later trajectory measurements are typically less reliable due to:
- Increasing distance from observation cameras
- Decreasing brightness as the meteorite burns
- Atmospheric distortion effects

The Halliday1996 Innisfree meteorite data is a good example where this pattern is observed.

### Solution
Custom weights can now be assigned to each measurement point in the trajectory. This allows later, less accurate measurements to have reduced influence on the optimization.

### Usage

```cpp
// Get trajectory data
size_t n_records = 0;
const real *time = nullptr, *v = nullptr, *h = nullptr;
meteorite->Trajectory(n_records, time, v, h);

// Generate exponentially decaying weights
// decay_factor: 0.0 = uniform weights, 1.0 = strong decay
auto weights = BasicFunctional::GenerateDecayingWeights(n_records, 0.5);

// Create weighted functional
L2Functional weighted_functional(meteorite, 1.0, 1.0, weights);
```

The `GenerateDecayingWeights` function creates weights that decrease exponentially:
- With `decay_factor = 0.5`, the last weight is ~22% of the first
- With `decay_factor = 1.0`, the last weight is ~5% of the first

Both `L1Functional` and `L2Functional` support weighted operation.

## 2. Functional Structure Visualization

### Problem
Understanding which functional configuration was used for a particular simulation is important for reproducibility and analysis.

### Solution
The `IFunctional` interface now includes `GetStructureDescription()` which returns a human-readable description of the functional's configuration.

### Output
When using `TrajectoryVisualizer`, a `functional_info.txt` file is created containing:
```
Functional Structure Information
================================

Functional: L2-Weighted
  Lambda_v: 1
  Lambda_h: 1
  Measurements: 9
  Weighted: Yes (decaying weights for later measurements)
  Weight range: [0.223, 1]
```

### Usage

```cpp
// Set the functional for visualization
visualizer.SetFunctional(&functional);
visualizer.Export(trajectories);
// This will create functional_info.txt automatically
```

## 3. Burn Time Penalty Functional

### Problem
If camera observations stop detecting the meteorite at time T, but a virtual meteorite continues burning beyond T (or stops before T), the simulation is unrealistic.

### Solution
`BurnTimePenaltyFunctional` wraps any functional and adds a penalty when the burn duration doesn't match expectations.

### Usage

```cpp
// Create base functional
auto base_func = std::make_shared<L2Functional>(meteorite);

// Wrap with burn time penalty
// expected_burn_time: when observations end (seconds)
// penalty_weight: how much to penalize mismatches (default 0.1)
BurnTimePenaltyFunctional penalty_func(
    base_func, 
    expected_burn_time,
    0.1  // penalty weight
);

// Use wrapped functional in solver
solver.Solve(generator, penalty_func, recorder);
```

The penalty is proportional to the time difference between expected and actual burn time.

## 4. Stage 0 Entry Angle Optimization

### Problem
Exploring the full parameter space for entry angle (gamma) can be computationally expensive, requiring simulation of millions of virtual meteorites.

### Solution
A preliminary Stage 0 performs a quick scan with:
- Fewer virtual meteorites (10,000 vs 1,000,000)
- Coarser timestep (0.01s vs 0.001s)
- Simpler numerical method (1-step vs 2-step Adams)

The optimal entry angle range is determined from the best 100 candidates, then refined with 10% margin for subsequent stages.

### Configuration

In `Main.cpp`:
```cpp
constexpr size_t   STAGE0_N_TOTAL    = 10000;    // Virtual meteoroids to test
constexpr size_t   STAGE0_N_TOP      = 100;      // Best candidates to keep
constexpr auto     STAGE0_METHOD     = NumericalAlgorithm::ONE_STEP_ADAMS;
constexpr real     STAGE0_DT         = (real)1e-2;  // Timestep in seconds
constexpr real     STAGE0_GAMMA_MARGIN = (real)0.1; // 10% margin around range
```

### Benefits
- Reduces Stage 1 parameter space by focusing on promising angles
- Minimal computational overhead (~1% of total time)
- Maintains accuracy by using sufficient margin

## Testing

Comprehensive tests verify the weighted functional behavior:
- `WeightedFunctionalTests.GenerateDecayingWeights` - Weight generation
- `WeightedFunctionalTests.WeightedL2Functional` - Basic weighted operation
- `WeightedFunctionalTests.WeightedVsNonWeightedWithSignificantError` - Comparative behavior
- `WeightedFunctionalTests.GetStructureDescription` - Visualization output

All tests pass successfully.

## Recommendations

1. **For Halliday1996-type data** with known decreasing accuracy:
   - Use `decay_factor = 0.5` as a good starting point
   - Experiment with values 0.3-0.7 to tune

2. **For high-quality data** with uniform accuracy:
   - Use `decay_factor = 0.0` (uniform weights)
   - Or use standard functionals without weights

3. **For preliminary exploration**:
   - Enable Stage 0 to optimize entry angle range
   - Adjust `STAGE0_N_TOTAL` based on available compute resources

4. **For burn time validation**:
   - Use `BurnTimePenaltyFunctional` when observation end time is known
   - Set `penalty_weight` based on relative importance (0.1-0.3 recommended)

## References

- Original issue: "Improving accuracy and stability"
- Test data: Halliday1996 Innisfree meteorite
- DOI: https://doi.org/10.1111/j.1945-5100.1996.tb02014.x
