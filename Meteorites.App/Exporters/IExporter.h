#pragma once
#include "Meteorites.Core/Defs.h"

#include "Meteorites.Core/IMeteorite.h"
#include "Meteorites.Core/Recorders/BufferingRecorder.h"

// Exports a set of meteorite trajectories to some container (text file, table, chart, etc)
class IExporter
{
  public:
    IExporter(const IExporter &) = delete;
    IExporter &operator =(const IExporter &) = delete;
    virtual ~IExporter() = default;

    // Sets a callback function that will be invoked after processing each portion of the data
    // For instance, 'step = 0.05f' means the callback will be called after every 5% progress
    virtual void OnProgress(const std::function<void(float)> &callback, float step) = 0;

    // Sets the output directory where the exported data will be stored
    // If the directory does not exist, it will be created
    // The actual names of files and their content depend on the exporter implementation
    virtual void SetDirectory(const std::filesystem::path &directory) = 0;

    // Sets optional metadata about the exporting data, such as the date of simulation and the meteorite information
    virtual void SetMetaData(const std::string &date,
                             const IMeteorite &meteorite) = 0;

    // Processes the provided meteoroid trajectories and stores them in the requested directory
    virtual void Export(const std::vector<MeteoroidTrajectory> &trajectories) = 0;

  protected:
    IExporter() = default;
};
