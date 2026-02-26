#pragma once
#include "Meteorites.Core/Defs.h"

#include "IFunctional.h"
#include "Meteorites.Core/IMeteorite.h"

// Wraps another functional and adds a penalty based on burn time mismatch
// The penalty is added when the virtual meteorite stops burning at a different time
// than when the camera stops detecting the real meteorite
class BurnTimePenaltyFunctional : public IFunctional
{
  public:
    BurnTimePenaltyFunctional(std::shared_ptr<IFunctional> base_functional,
                              real expected_burn_time,
                              real penalty_weight = (real)0.1)
      : base_functional_(std::move(base_functional)),
        expected_burn_time_(expected_burn_time),
        penalty_weight_(penalty_weight)
    {
      assert(base_functional_ != nullptr);
      assert(expected_burn_time_ > 0.0);
      assert(penalty_weight_ >= 0.0);
    }

    // The member of 'IFunctional'
    virtual std::string Name() const override final
    {
      return base_functional_->Name() + "+BurnTimePenalty";
    }

    // The member of 'IFunctional'
    virtual void GetTimeStamps(size_t &num, const real *&values) const override final
    {
      base_functional_->GetTimeStamps(num, values);
    }

    // The member of 'IFunctional'
    virtual double Compute(size_t num, const real *v, const real *h) const override final
    {
      double base_error = base_functional_->Compute(num, v, h);
      
      // If the virtual meteorite didn't reach all expected timestamps,
      // it burned too early - add penalty
      size_t expected_num = 0;
      const real *expected_times = nullptr;
      GetTimeStamps(expected_num, expected_times);
      
      if (num < expected_num)
      {
        // Meteorite burned earlier than expected
        real actual_burn_time = (num > 0) ? expected_times[num - 1] : (real)0.0;
        real time_diff = expected_burn_time_ - actual_burn_time;
        double penalty = penalty_weight_ * std::abs(time_diff) / expected_burn_time_;
        return base_error + penalty;
      }
      else
      {
        // Meteorite lasted at least as long as expected (normal case)
        return base_error;
      }
    }

    // Returns a description of the functional structure for visualization
    virtual std::string GetStructureDescription() const override
    {
      std::ostringstream oss;
      oss << "Burn Time Penalty Functional\n";
      oss << "  Base Functional:\n";
      std::string base_desc = base_functional_->GetStructureDescription();
      // Indent each line of base description
      size_t pos = 0;
      while (pos < base_desc.length())
      {
        size_t next_pos = base_desc.find('\n', pos);
        if (next_pos == std::string::npos)
        { next_pos = base_desc.length(); }
        oss << "    " << base_desc.substr(pos, next_pos - pos) << "\n";
        pos = next_pos + 1;
      }
      oss << "  Expected burn time: " << expected_burn_time_ << " seconds\n";
      oss << "  Penalty weight: " << penalty_weight_ << "\n";
      return oss.str();
    }

  private:
    std::shared_ptr<IFunctional> base_functional_;
    real expected_burn_time_;
    real penalty_weight_;
};
