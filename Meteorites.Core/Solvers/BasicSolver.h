#pragma once
#include "Meteorites.Core/Defs.h"

#include "ISolver.h"


// Implements some common logic
class BasicSolver : public ISolver
{
  public:
    // Iterates through virtual meteoroids provided either as a vector or generated on-the-fly,
    // extracting and returning them in batches for processing
    // Note: MeteoroidEnumerator holds references to the underlying container or generator
    // Do not modify the container or generator while iterating over it!
    class MeteoroidEnumerator
    {
      public:
        // Tries to extract next batch of meteoroids
        // Returns true if this batch is not empty
        bool MoveNext(const VirtualMeteoroid *&batch, size_t &size);

      private:
        MeteoroidEnumerator() = delete;
        MeteoroidEnumerator(const MeteoroidEnumerator &) = delete;
        MeteoroidEnumerator &operator =(const MeteoroidEnumerator &) = delete;

        MeteoroidEnumerator(const VirtualMeteoroid &meteoroid);
        MeteoroidEnumerator(const std::vector<VirtualMeteoroid> &meteoroids, size_t batch_size);
        MeteoroidEnumerator(IMeteoroidGenerator &gen, size_t batch_size);

        std::function<bool(const VirtualMeteoroid *&, size_t &)> impl_;

      friend class BasicSolver;
    };

    // ISolver method
    virtual void Configure(NumericalAlgorithm alg, real dt, real timeout) override final;

    // ISolver method
    virtual void Solve(const VirtualMeteoroid &problem,
                       const IFunctional &functional,
                       ISimulationRecorder &results) override final;

    // ISolver method
    virtual void Solve(const std::vector<VirtualMeteoroid> &problems,
                       const IFunctional &functional,
                       ISimulationRecorder &results) override final;

    // ISolver method
    virtual void Solve(IMeteoroidGenerator &generator,
                       const IFunctional &functional,
                       ISimulationRecorder &results) override final;

  protected:
    BasicSolver() = default;

    virtual size_t BatchSize() const = 0;

    virtual void SolveAny(MeteoroidEnumerator &problems,
                          const IFunctional &functional,
                          ISimulationRecorder &results) = 0;

    NumericalAlgorithm Algorithm() const { return algorithm_; }

    size_t Steps() const { return steps_; }

    real Dt() const { return dt_; }

    real Timeout() const { return timeout_; }

  private:
    NumericalAlgorithm algorithm_ = NumericalAlgorithm::ONE_STEP_ADAMS;
    size_t steps_ = 1;
    real dt_ = (real)0.001;
    real timeout_ = (real)1000.0;
};
