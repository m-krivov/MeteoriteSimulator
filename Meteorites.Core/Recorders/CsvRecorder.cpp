#include "CsvRecorder.h"


namespace
{

std::string FormatCsvName(const std::filesystem::path &directory,
                          const std::string &id, size_t cur)
{
  std::stringstream ss;
  ss << id << "_" << cur << ".csv";
  return (directory / ss.str()).string();
}

} // unnamed namespace

CsvRecorder::CsvRecorder(const std::filesystem::path &directory,
                         const std::string &id, real dt)
  : dt_(dt), t_next_((real)0.0), directory_(directory), id_(id), cur_(0)
{
  assert(std::filesystem::exists(directory));
  assert(dt >= (real)0.0);
}

CsvRecorder::~CsvRecorder()
{
  if (file_.good())
  { file_.close(); }
}

real CsvRecorder::Started(const VirtualMeteoroid &problem)
{
  assert(!file_.is_open());
  auto name = FormatCsvName(directory_, id_, cur_);
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

real CsvRecorder::Store(real t, real m, real v, real h, real l, real gamma)
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

void CsvRecorder::Finished(Reason reason, double accuracy)
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
