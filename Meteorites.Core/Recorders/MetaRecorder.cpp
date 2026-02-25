#include "MetaRecorder.h"


MetaRecorder::MetaRecorder(size_t n_best, size_t buffer_size)
  : n_best_(n_best), buffer_size_(buffer_size),
    accuracy_threshold_(std::numeric_limits<double>::max())
{
  if (n_best == 0)
  { throw std::runtime_error("'n_best' must be a positive number"); }
  if (buffer_size < n_best)
  { throw std::runtime_error("'buffer_size' must be greater than 'n_best'"); }
}

real MetaRecorder::Started(const VirtualMeteoroid &problem)
{
  assert(!current_.has_value());
  current_.emplace(problem);
  return std::numeric_limits<real>::max();
}

real MetaRecorder::Store(real, real, real, real, real, real)
{
  return std::numeric_limits<real>::max();
}

namespace
{

void SelectBest(std::vector<MeteoroidSummary> &problems, size_t n_best)
{
  std::sort(problems.begin(), problems.end(),
            [](const MeteoroidSummary &el1,
               const MeteoroidSummary &el2) -> bool
            {
              return el1.Accuracy() < el2.Accuracy();
            });

  if (problems.size() > n_best)
  { problems.resize(n_best); }
}

} // unnamed namespace

void MetaRecorder::Finished(Reason reason, double accuracy)
{
  assert(current_.has_value());

  // If solution is not good enough, simply reject it
  if (problems_.size() > n_best_ && accuracy >= accuracy_threshold_)
  {
    current_.reset();
    return;
  }

  // Otherwise, store solution and update list of the best cases
  problems_.emplace_back(MeteoroidSummary(current_.value(), reason, accuracy));
  if (problems_.size() > buffer_size_)
  {
    SelectBest(problems_, n_best_);
    assert(problems_.size() == n_best_);
    accuracy_threshold_ = problems_.back().Accuracy();
  }
  current_.reset();
}

void MetaRecorder::MoveTo(std::vector<MeteoroidSummary> &results)
{
  if (problems_.size() > n_best_)
  { SelectBest(problems_, n_best_); }

  results = std::move(problems_);
  accuracy_threshold_ = std::numeric_limits<double>::max();
}
