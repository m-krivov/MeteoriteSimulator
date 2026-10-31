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

#if defined(METEORITES_FP64)
  typedef double real;
  constexpr std::string_view METEORITES_PRECISION = "FP64";
#elif defined(METEORITES_FP32)
  typedef float real;
  constexpr std::string_view METEORITES_PRECISION = "FP32";
#else
  #error The project was not configured properly, METEORITES_PRECISION must be set to either 'fp32' or 'fp64'
#endif

#if !defined(M_PI)
  constexpr double M_PI = 3.1415926535897932384626433832795;
#endif

// Stub for CUDA solvers
#if !defined(DEVICE)
  #define DEVICE
#endif
