#include "BasicExporter.h"

void BasicExporter::OnProgress(const std::function<void(float)> &callback, float step)
{
  assert(step <= 1.0f);
  progress_callback_  = callback;
  progress_step_      = step;
  progress_threshold_ = step;
}

void BasicExporter::SetDirectory(const std::filesystem::path &directory)
{
  directory_ = directory;
  if (!std::filesystem::exists(directory_) &&
      !std::filesystem::create_directories(directory_))
  {
    std::ostringstream oss;
    oss << "Failed to create a directory ('" << directory.string() << "')";
    throw std::runtime_error(oss.str());
  }
}

void BasicExporter::SetMetaData(const std::string &date,
                                const std::shared_ptr<const IMeteorite> &meteorite)
{
  date_      = date;
  meteorite_ = meteorite;
}

const IMeteorite &BasicExporter::Meteorite() const
{
  assert(meteorite_ != nullptr);
  return *meteorite_;
}

void BasicExporter::UpdateProgress(size_t current, size_t total)
{
  // TODO: move this logic to a separate interface and update 'BasicMeteoroidGenerator'
  assert(current <= total);
  if (!progress_callback_)
  { return; }

  auto ratio = (float)current / total;
  while (ratio >= progress_threshold_)
  {
    progress_callback_(progress_threshold_);
    progress_threshold_ += progress_step_;
  }

  if (ratio < progress_threshold_ && current == total)
  { progress_callback_(1.0f); }
}
