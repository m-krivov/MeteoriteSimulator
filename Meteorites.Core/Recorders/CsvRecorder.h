#pragma once
#include "Meteorites.Core/Defs.h"

#include "ISimulationRecorder.h"

// Stores results into a separate CSV file
class CsvRecorder : public ISimulationRecorder
{
  public:
    CsvRecorder(const std::string &id, real dt);
    ~CsvRecorder();

    // The member of 'ISimulationRecorder'
    virtual real Started(const VirtualMeteoroid &problem) override;

    // The member of 'ISimulationRecorder'
    virtual real Store(real t, real m, real v, real h, real l, real gamma) override;

    // The member of 'ISimulationRecorder'
    virtual bool NeedTrajectory() const override { return true; }

    // The member of 'ISimulationRecorder'
    virtual void Finished(Reason reason, double accuracy) override;

  private:
    real dt_, t_next_;
    std::string id_;
    size_t cur_;
    std::ofstream file_;
};
