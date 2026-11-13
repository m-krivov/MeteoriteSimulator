#pragma once
#include "Meteorites.Core/Defs.h"

// The roughly predicted number of operations issued for some action, can be used to estimate algorithmic FLOPS
// Not a hardware instruction counter!
struct OperationCounter
{
  // add, sub, mul
  uint64_t arithmetic{};

  // div, rcp
  uint64_t division{};

  // sin, cos, exp, log, pow, sqrt
  uint64_t special{};
};

inline constexpr OperationCounter operator +(const OperationCounter &lhs, const OperationCounter &rhs)
{
  return OperationCounter{
      .arithmetic = lhs.arithmetic + rhs.arithmetic,
      .division   = lhs.division   + rhs.division,
      .special    = lhs.special    + rhs.special
  };
}

inline constexpr OperationCounter operator -(const OperationCounter &lhs, const OperationCounter &rhs)
{
  return OperationCounter{
      .arithmetic = lhs.arithmetic - rhs.arithmetic,
      .division   = lhs.division   - rhs.division,
      .special    = lhs.special    - rhs.special
  };
}

inline constexpr OperationCounter operator *(const OperationCounter &counter, uint64_t factor)
{
  return OperationCounter{
      .arithmetic = counter.arithmetic * factor,
      .division   = counter.division   * factor,
      .special    = counter.special    * factor
  };
}

inline constexpr OperationCounter operator *(uint64_t factor, const OperationCounter &counter)
{
  return counter * factor;
}

inline constexpr OperationCounter &operator +=(OperationCounter &lhs, const OperationCounter &rhs)
{
  lhs.arithmetic += rhs.arithmetic;
  lhs.division   += rhs.division;
  lhs.special    += rhs.special;
  return lhs;
}

inline constexpr OperationCounter &operator -=(OperationCounter &lhs, const OperationCounter &rhs)
{
  lhs.arithmetic -= rhs.arithmetic;
  lhs.division   -= rhs.division;
  lhs.special    -= rhs.special;
  return lhs;
}

inline constexpr OperationCounter &operator *=(OperationCounter &lhs, uint64_t factor)
{
  lhs.arithmetic *= factor;
  lhs.division   *= factor;
  lhs.special    *= factor;
  return lhs;
}
