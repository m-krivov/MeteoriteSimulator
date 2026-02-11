#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <future>
#include <memory>
#include <numeric>
#include <optional>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <sstream>
#include <stdint.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>

typedef double real;

#if !defined(M_PI)
  constexpr double M_PI = 3.1415926535897932384626433832795;
#endif

// Stub for CUDA solvers
#if !defined(DEVICE)
  #define DEVICE
#endif
