#pragma once

#include "Meteorites.CudaSolvers/CudaDefs.h"
#include "Meteorites.Core/Meteoroids/VirtualMeteoroid.h"

class MeteoroidsManager
{
  public:
    MeteoroidsManager() = delete;
    MeteoroidsManager(const MeteoroidsManager &) = delete;
    MeteoroidsManager &operator =(const MeteoroidsManager &) = delete;

    MeteoroidsManager(size_t m_bests,
                      size_t total_threads,
                      size_t best_meteoroids_per_thread);
    ~MeteoroidsManager();

    MeteoroidDeviation* GetDeviceDeviationsBuffer();
    uint64_t*           GetDeviceSeeds();

    void GenerateSeeds();

    void UpdateTopMeteoroids();

    real GetActualThresholdDeviation();

    std::vector<std::pair<VirtualMeteoroid, real> >
    GetTopMeteoroids(const TrajectoryPoint &ref_point0);
  
  private:
    // Top m_bests meteoroids share one memory block with deviations buffer
    // on GPU to simplify sorting between batch launches
    CudaPtr<MeteoroidDeviation> coalesced_deviations_;
    CudaPtr<uint64_t>           coalesced_seeds_;
    CudaPtr<uint64_t>     device_seeds_;
    std::vector<uint64_t> host_seeds_;
    size_t m_bests_;
    size_t total_threads_;
    size_t best_meteoroids_per_thread_;
};