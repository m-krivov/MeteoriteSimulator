#include "ResultFormatters.h"

//------------------------
//--- IResultFormatter ---
//------------------------

IResultFormatter::Reason IResultFormatter::Classify(real t, real M, real h)
{
  if (M <= (real)0.01)
  {
    return IResultFormatter::Reason::Burnt;
  }
  else if (h <= (real)0.0)
  {
    return IResultFormatter::Reason::Collided;
  }
  else
  {
    return IResultFormatter::Reason::Timeouted;
  }
}

//--------------------
//--- CsvFromatter ---
//--------------------

namespace
{

std::string FormatCsvName(const std::string &id, size_t cur)
{
  std::stringstream ss;
  ss << id << "_" << cur << ".csv";
  return ss.str();
}

} // unnamed namespace

CsvFromatter::CsvFromatter(const std::string &id, real dt)
  : dt_(dt), t_next_((real)0.0), id_(id), cur_(0)
{
  assert(dt >= (real)0.0);
}

CsvFromatter::~CsvFromatter()
{
  if (file_.good())
  { file_.close(); }
}

real CsvFromatter::Started(const Case &problem)
{
  assert(!file_.is_open());
  auto name = FormatCsvName(id_, cur_);
  file_.open(name);
  if (!file_.good())
  {
    std::stringstream ss;
    ss << "failed to write data to the file '" << name << "'";
    throw std::runtime_error(ss.str());
  }

  file_ << "Cd:,"     << problem.Cd     << "," << std::endl;
  file_ << "Cl:,"     << problem.Cl     << "," << std::endl;
  file_ << "Ch:,"     << problem.Ch     << "," << std::endl;
  file_ << "H:,"      << problem.H      << "," << std::endl;
  file_ << "Rho:,"    << problem.Rho    << "," << std::endl;
  file_ << "Gamma0:," << problem.Gamma0 << "," << std::endl;
  file_ << "M0:,"     << problem.M0     << "," << std::endl;
  file_ << std::endl;
  file_ << "Time," << "Mass," << "Velocity," << "Height," << "Distance," << "Gamma," << std::endl;

  t_next_ = (real)0.0f;

  return t_next_;
}

real CsvFromatter::Store(real t, real m, real v, real h, real l, real gamma)
{
  assert(file_.good());
  assert(file_.is_open());

  if (t >= t_next_)
  {
    t_next_ += dt_;
    file_ << t << ','
          << m << ','
          << v << ','
          << h << ','
          << l << ','
          << gamma << ',' << std::endl;
  }

  return t_next_;
}

void CsvFromatter::Finished(Reason reason, double accuracy)
{
  assert(file_.good());
  assert(file_.is_open());

  switch (reason)
  {
    case Reason::Burnt:
      file_ << "Burnt";
      break;

    case Reason::Collided:
      file_ << "Collided";
      break;

    case Reason::Timeouted:
      file_ << "Timeouted";
      break;

    default:
      assert(false);
  }

  file_ << ", accuracy " << accuracy << std::endl;
  file_.close();
  cur_ += 1;
}

//--------------------------
//--- BufferingFormatter ---
//--------------------------

BufferingFormatter::BufferingFormatter(real dt)
  : dt_(dt), t_next_((real)0.0)
{
  if (dt < (real)0.0)
  { throw std::runtime_error("'dt' must be declared as positive number or zero"); }
}

real BufferingFormatter::Started(const Case &problem)
{
  t_next_   = (real)0.0;
  logs_.emplace_back(Log(problem));
  return t_next_;
}

real BufferingFormatter::Store(real t, real m, real v, real h, real l, real gamma)
{
  if (t >= t_next_)
  {
    assert(!logs_.empty());
    auto &log = logs_[logs_.size() - 1];
    log.records.emplace_back(Record(t, m, v, h, l, gamma));
    t_next_ += dt_;
  }

  return t_next_;
}

void BufferingFormatter::Finished(Reason reason, double accuracy)
{
  assert(!logs_.empty());
  auto &log = logs_[logs_.size() - 1];
  log.reason = reason;
  log.accuracy = accuracy;
}

//---------------------
//--- MetaFormatter ---
//---------------------

MetaFormatter::MetaFormatter(size_t n_best, size_t buffer_size)
  : n_best_(n_best), buffer_size_(buffer_size),
    accuracy_threshold_(std::numeric_limits<double>::max())
{
  if (n_best == 0)
  { throw std::runtime_error("'n_best' must be a positive number"); }
  if (buffer_size < n_best)
  { throw std::runtime_error("'buffer_size' must be greater than 'n_best'"); }
}

real MetaFormatter::Started(const Case &problem)
{
  problems_.emplace_back(std::make_pair(problem, std::numeric_limits<double>::max()));
  return std::numeric_limits<real>::max();
}

real MetaFormatter::Store(real, real, real, real, real, real)
{
  return std::numeric_limits<real>::max();
}

namespace
{

void SelectBest(std::vector<std::pair<Case, double> > &problems, size_t n_best)
{
  std::sort(problems.begin(), problems.end(),
            [](const std::pair<Case, double> &el1,
                const std::pair<Case, double> &el2) -> bool
            {
              return el1.second < el2.second;
            });

  if (problems.size() > n_best)
  { problems.resize(n_best); }
}

} // unnamed namespace

void MetaFormatter::Finished(Reason reason, double accuracy)
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

void MetaFormatter::ExportAndReset(std::vector<std::pair<Case, double> > &results)
{
  if (problems_.size() >= n_best_)
  { SelectBest(problems_, n_best_); }

  results = std::move(problems_);
  accuracy_threshold_ = std::numeric_limits<double>::max();
}
