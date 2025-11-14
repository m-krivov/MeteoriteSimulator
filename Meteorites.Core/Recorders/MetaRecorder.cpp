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
  problems_.emplace_back(std::make_pair(problem, std::numeric_limits<double>::max()));
  return std::numeric_limits<real>::max();
}

real MetaRecorder::Store(real, real, real, real, real, real)
{
  return std::numeric_limits<real>::max();
}

namespace
{

void SelectBest(std::vector<std::pair<VirtualMeteoroid, double> > &problems, size_t n_best)
{
  std::sort(problems.begin(), problems.end(),
            [](const std::pair<VirtualMeteoroid, double> &el1,
               const std::pair<VirtualMeteoroid, double> &el2) -> bool
            {
              return el1.second < el2.second;
            });

  if (problems.size() > n_best)
  { problems.resize(n_best); }
}

} // unnamed namespace

void MetaRecorder::Finished(Reason reason, double accuracy)
{
  assert(!problems_.empty());

  // If solution is not good enough, simply reject it
  if (problems_.size() > n_best_ && accuracy >= accuracy_threshold_)
  {
    problems_.pop_back();
    return;
  }

  // Otherwise, store solution and update list of the best cases
  problems_[problems_.size() - 1].second = accuracy;
  if (problems_.size() >= buffer_size_)
  {
    SelectBest(problems_, n_best_);
    assert(problems_.size() == n_best_);
    accuracy_threshold_ = problems_[problems_.size() - 1].second;
  }
}

void MetaRecorder::ExportAndReset(std::vector<std::pair<VirtualMeteoroid, double> > &results)
{
  if (problems_.size() >= n_best_)
  { SelectBest(problems_, n_best_); }

  results = std::move(problems_);
  accuracy_threshold_ = std::numeric_limits<double>::max();
}
