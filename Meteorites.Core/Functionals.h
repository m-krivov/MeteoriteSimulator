#pragma once
#include "Defs.h"
#include "IMeteorite.h"

// Interface that measures some kind of 'distance' between virtual meteorite and a real one
class IFunctional
{
  public:
    IFunctional(const IFunctional &) = delete;
    IFunctional &operator =(const IFunctional &) = delete;
    virtual ~IFunctional() { }

    // Returns timestamps that are used as arguments
    // You may expect that at least one timestamp is provided
    virtual void GetTimeStamps(size_t &num, const real *&values) const = 0;

    // Computes functional's value for the specified arguments
    // The actual number of arguments may differ from the expected: such values just increases error
    virtual double Compute(size_t num, const real *v, const real *h) const = 0;

  protected:
    IFunctional() = default;
};

// Does nothing, just a placeholder for tests and debugging
class FakeFunctional : public IFunctional
{
  public:
    FakeFunctional()
      : num_stub(1), values_stub((real)0.0)
    { }

    // IFunctional member
    virtual void GetTimeStamps(size_t &num, const real *&values) const override
    {
      num = num_stub;
      values = &values_stub;
    }

    // IFunctional member
    virtual double Compute(size_t num, const real *v, const real *h) const override
    { return std::numeric_limits<double>::max(); }

  private:
    size_t num_stub;
    real values_stub;
};

// Computes L2 norm
class L2Functional : public IFunctional
{
  public:
    L2Functional(const IMeteorite &meteorite);

    // IFunctional member
    virtual void GetTimeStamps(size_t &num, const real *&values) const override;

    // IFunctional member
    virtual double Compute(size_t num, const real *v, const real *h) const override;

  private:
    std::vector<real> time_, v_, h_;
};
