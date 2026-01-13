#pragma once
#include "IExporter.h"

class BasicExporter : public IExporter
{
  public:
    BasicExporter(const BasicExporter &) = delete;
    BasicExporter &operator =(const BasicExporter &) = delete;

    // The member of 'IExporter'
    virtual void OnProgress(const std::function<void(float)> &callback, float step) override final;

    // The member of 'IExporter'
    virtual void SetDirectory(const std::filesystem::path &directory) override final;

    // The member of 'IExporter'
    virtual void SetMetaData(const std::string &date,
                             const std::shared_ptr<const IMeteorite> &meteorite) override final;

  protected:
    BasicExporter() = default;

    const std::filesystem::path &Directory() const { return directory_; }

    const std::string &Date() const { return date_; }

    const IMeteorite &Meteorite() const;

    void UpdateProgress(size_t current, size_t total);

  private:
    std::filesystem::path directory_;
    std::string date_;
    std::shared_ptr<const IMeteorite> meteorite_;
    std::function<void(float)> progress_callback_;
    float progress_step_{0.0f}, progress_threshold_{0.0f};
};
