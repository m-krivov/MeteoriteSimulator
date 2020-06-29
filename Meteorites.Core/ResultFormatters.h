#pragma once
#include "Defs.h"
#include "Case.h"

// Interface that saves the numerically computed values
class IResultFormatter
{
  public:
    // Due to which reason the simulation was stopped?
    enum class Reason : uint32_t
    {
      // Field is not initialized
      NA          = 0,

      // No mass left, meteorite is burnt
      Burnt       = 1,

      // Height is equal or lesser than zero, B-O-O-O-M!
      Collided    = 2,

      // Looks like something is wrong, meteorite flies for too long time
      // Did we process some UFO trajectory instead of meteorite?
      Timeouted   = 3
    };

    IResultFormatter(const IResultFormatter &) = delete;
    IResultFormatter &operator =(const IResultFormatter &) = delete;
    virtual ~IResultFormatter() { }

    // Notifies that simulation is started, data will be coming soon
    // Provides problem that are going to be simulated
    // Returns the expected time of the first record (take a look at 'Store()')
    virtual real Started(const Case &problem) = 0;

    // Makes single record about meteorite state
    // Returns the proposed time of the next record, e.g.
    //    '0.0f'     - accept any,
    //    't + 0.1f' - skip '0.1' seconds,
    //    'MAX_FLT'  - don't need new information
    // It's a recommendation: you may ignore this value and continue spamming
    virtual real Store(real t, real m, real v, real h, real l, real gamma) = 0;

    // Notifies that computations are ended due to some reason
    // Value 'accuracy' defines how accurately this virtual meteorite describes the real one
    virtual void Finished(Reason reason, double accuracy) = 0;

  protected:
    IResultFormatter() = default;
};

// Stores results into a separate CSV file
class CsvFromatter : public IResultFormatter
{
  public:
    CsvFromatter(const std::string &id, real dt);
    ~CsvFromatter();

    // IResultFormatter member
    virtual real Started(const Case &problem) override;

    // IResultFormatter member
    virtual real Store(real t, real m, real v, real h, real l, real gamma) override;

    // IResultFormatter member
    virtual void Finished(Reason reason, double accuracy) override;

  private:
    real dt_, t_next_;
    std::string id_;
    size_t cur_;
    std::ofstream file_;
};

// Stores all data that were passed to this formatter
class BufferingFormatter : public IResultFormatter
{
  public:
    // The named analogue of 'std::tuple<real, real, real, real, real>'
    struct Record
    {
      Record() : t(0.0f), M(0.0f), V(0.0f), h(0.0f), l(0.0f), gamma(0.0f) { }
      Record(real t_, real m_, real v_, real h_, real l_, real gamma_)
        : t(t_), M(m_), V(v_), h(h_), l(l_), gamma(gamma_) { }
      Record(const Record &) = default;
      Record &operator =(const Record &) = default;

      real t, M, V, h, l, gamma;
    };

    // The complete information about single meteorite 
    struct Log
    {
      Log() : reason(Reason::NA), accuracy(std::numeric_limits<double>::max()) { }
      Log(const Case &problem_)
        : problem(problem_), reason(Reason::NA), accuracy(std::numeric_limits<double>::max()) { }
      Log(const Log &) = default;
      Log &operator =(const Log &) = default;

      Case problem;
      std::vector<Record> records;
      Reason reason;
      double accuracy;
    };

    // Value 'dt' defines how often the incoming records must be stored, in seconds
    BufferingFormatter(real dt);

    // IResultFormatter member
    virtual real Started(const Case &problem) override;

    // IResultFormatter member
    virtual real Store(real t, real m, real v, real h, real l, real gamma) override;

    // IResultFormatter member
    virtual void Finished(Reason reason, double accuracy) override;

    // Information about all trajectories that were passed to this formatter
    const std::vector<Log> &Logs() const { return logs_; }

    // Skips information about all recorded informations
    void Reset() { logs_.clear(); }

  private:
    real dt_, t_next_;
    std::vector<Log> logs_;
};

// Accumulates results from multiple experiments and stores meta data about the best ones
// Skips all information about trajectory: keeps only case and functional's value
// Do you need values for velocity/mass/etc? Just recompute the case with different formatter
class MetaFormatter : public IResultFormatter
{
  public:
    MetaFormatter(size_t n_best, size_t buffer_size);

    // IResultFormatter member
    virtual real Started(const Case &problem) override;

    // IResultFormatter member
    virtual real Store(real t, real m, real v, real h, real l, real gamma) override;

    // IResultFormatter member
    virtual void Finished(Reason reason, double accuracy) override;

    // Extracts N best cases that were reported to this formatter
    // Resets the internal ratings
    void ExportAndReset(std::vector<std::pair<Case, double> > &results);

  private:
    size_t n_best_, buffer_size_;
    double accuracy_threshold_;
    std::vector<std::pair<Case, double> > problems_;
};
