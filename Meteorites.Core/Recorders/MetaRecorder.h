#pragma once
#include "Meteorites.Core/Defs.h"

#include "ISimulationRecorder.h"

// Accumulates results from multiple experiments and stores meta data about the best ones
// Skips all information about trajectory: keeps only case and functional's value
// Do you need values for velocity/mass/etc? Just recompute the case with different formatter
class MetaRecorder : public ISimulationRecorder
{
  public:
    MetaRecorder(size_t n_best, size_t buffer_size);

    // The member of 'ISimulationRecorder'
    virtual real Started(const VirtualMeteoroid &problem) override;

    // The member of 'ISimulationRecorder'
    virtual real Store(real t, real m, real v, real h, real l, real gamma) override;

    // The member of 'ISimulationRecorder'
    virtual bool NeedTrajectory() const override { return false;}

    // The member of 'ISimulationRecorder'
    virtual void Finished(Reason reason, double accuracy) override;

    // Extracts N best cases that were reported to this formatter
    // Resets the internal ratings
    void ExportAndReset(std::vector<std::pair<VirtualMeteoroid, double> > &results);

  private:
    size_t n_best_, buffer_size_;
    double accuracy_threshold_;
    std::vector<std::pair<VirtualMeteoroid, double> > problems_;
};
