#pragma once
#include "Meteorites.Core/Defs.h"

#include "IFunctional.h"
#include "Meteorites.Core/IMeteorite.h"

// Basic class for the classical functionals like C, L1 or L2
class BasicFunctional : public IFunctional
{
  public:

    // The member of 'IFunctional'
    virtual void GetTimeStamps(size_t &num, const real *&values) const override final;

    // Returns a description of the functional structure for visualization
    virtual std::string GetStructureDescription() const override
    {
      std::ostringstream oss;
      oss << "Functional: " << Name() << "\n";
      oss << "  Lambda_v: " << lambda_v_ << "\n";
      oss << "  Lambda_h: " << lambda_h_ << "\n";
      oss << "  Measurements: " << time_.size() << "\n";
      if (HasWeights())
      {
        oss << "  Weighted: Yes (decaying weights for later measurements)\n";
        oss << "  Weight range: [" << *std::min_element(weights_.begin(), weights_.end()) 
            << ", " << *std::max_element(weights_.begin(), weights_.end()) << "]\n";
      }
      else
      {
        oss << "  Weighted: No (uniform weights)\n";
      }
      return oss.str();
    }

    // Generate weights that decrease for later measurements
    // This is useful because later trajectory values are typically less accurate
    // decay_factor: how much weight decreases per measurement (0.0 = uniform, 1.0 = exponential decay)
    static std::vector<real> GenerateDecayingWeights(size_t count, real decay_factor = (real)0.5);

  protected:
    BasicFunctional(const std::shared_ptr<const IMeteorite> &meteorite,
                    real lambda_v, real lambda_h);

    // Constructor with custom weights for each measurement point
    BasicFunctional(const std::shared_ptr<const IMeteorite> &meteorite,
                    real lambda_v, real lambda_h,
                    const std::vector<real> &weights);

    real LambdaV() const { return lambda_v_; }

    real LambdaH() const { return lambda_h_; }

    const std::vector<real> &Time() const { return time_; }

    const std::vector<real> &Velocity() const { return v_; }

    const std::vector<real> &Height() const { return h_; }

    const std::vector<real> &Weights() const { return weights_; }

    bool HasWeights() const { return !weights_.empty(); }

  private:
    real lambda_v_ = (real)1.0;
    real lambda_h_ = (real)1.0;
    std::vector<real> time_, v_, h_;
    std::vector<real> weights_;  // Optional weights for each measurement point
};
