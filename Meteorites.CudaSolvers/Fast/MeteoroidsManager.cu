#include "MeteoroidsManager.h"

#include "Meteorites.CudaSolvers/Random.h"
#include "CudaMeteoroidsGenerator.h"

MeteoroidsManager::MeteoroidsManager(size_t m_bests,
                                     size_t total_threads,
                                     size_t best_meteoroids_per_thread)
      : m_bests_(m_bests), total_threads_(total_threads),
        best_meteoroids_per_thread_(best_meteoroids_per_thread)
{
  HANDLE_ERROR(CudaAlloc(coalesced_deviations_, m_bests_ + best_meteoroids_per_thread_ * total_threads_));
  MeteoroidDeviation initial_m_bests[m_bests];
  for (size_t i = 0; i < m_bests; i++) { initial_m_bests[i] = { std::numeric_limits<real>::max(), 0 }; }
  HANDLE_ERROR(cudaMemcpy(coalesced_deviations_.get(), initial_m_bests,
                          sizeof(MeteoroidDeviation) * m_bests, cudaMemcpyHostToDevice));

  HANDLE_ERROR(CudaAlloc(coalesced_seeds_, m_bests + best_meteoroids_per_thread_ * total_threads_));
  HANDLE_ERROR(cudaMemset(coalesced_seeds_.get(), 0, sizeof(uint64_t) * m_bests_));

  HANDLE_ERROR(CudaAlloc(device_seeds_, total_threads_));

  host_seeds_.resize(total_threads_);
}

MeteoroidsManager::~MeteoroidsManager()
{
  try
  {
    coalesced_deviations_.reset();
    coalesced_seeds_.reset();
    device_seeds_.reset();
  }
  catch(std::exception &)
  {
  }
}

MeteoroidDeviation* MeteoroidsManager::GetDeviceDeviationsBuffer()
{
  return coalesced_deviations_.get() + m_bests_;
}

uint64_t* MeteoroidsManager::GetDeviceSeeds()
{
  return device_seeds_.get();
}

void MeteoroidsManager::GenerateSeeds()
{
  Random::State rootState;
  Random::Initialize(rootState);

  std::vector<Random::State> tempStates(host_seeds_.size());

  Random::Multiply(rootState, tempStates);

  for (size_t i = 0; i < host_seeds_.size(); i++)
  {
    uint32_t part1 = Random::Next(tempStates[i]);
    uint32_t part2 = Random::Next(tempStates[i]);
    host_seeds_[i] = (static_cast<uint64_t>(part1) << 32) | part2;
  }

  HANDLE_ERROR(cudaMemcpy(device_seeds_.get(), host_seeds_.data(),
                          sizeof(uint64_t) * total_threads_, cudaMemcpyHostToDevice));
}

void MeteoroidsManager::UpdateTopMeteoroids()
{
  // Restore seeds corresponding to actual deviations before global sorting
  std::vector<uint64_t> restored_seeds(best_meteoroids_per_thread_ * total_threads_);
  for (size_t tid = 0; tid < total_threads_; tid++)
  {
    for (size_t i = 0; i < best_meteoroids_per_thread_; i++)
    {
      restored_seeds[tid * best_meteoroids_per_thread_ + i] = host_seeds_[tid];
    }
  }
  HANDLE_ERROR(cudaMemcpy(coalesced_seeds_.get() + m_bests_, restored_seeds.data(),
                          sizeof(uint64_t) * best_meteoroids_per_thread_ * total_threads_, cudaMemcpyHostToDevice));

  // Sorting

  // Causes strange performance drops in the main kernel
  /*
  auto zip_begin = thrust::make_zip_iterator(thrust::make_tuple(coalesced_deviations_.get(),
                                                                coalesced_seeds_.get()));
  auto zip_end = zip_begin + (m_bests_ + best_meteoroids_per_thread_ * total_threads_);

  thrust::sort(thrust::device, zip_begin, zip_end,
    [] __device__ (const auto& a, const auto& b) { return thrust::get<0>(a).dev < thrust::get<0>(b).dev; });

  HANDLE_ERROR(cudaDeviceSynchronize());
  */

  // Temporary solution 
  size_t coalesced_buf_size = m_bests_ + best_meteoroids_per_thread_ * total_threads_;

  std::vector<MeteoroidDeviation> host_deviations(coalesced_buf_size);
  std::vector<uint64_t>           host_seeds(coalesced_buf_size);

  HANDLE_ERROR(cudaMemcpy(host_deviations.data(), coalesced_deviations_.get(),
                          sizeof(MeteoroidDeviation) * coalesced_buf_size, cudaMemcpyDeviceToHost));
  HANDLE_ERROR(cudaMemcpy(host_seeds.data(), coalesced_seeds_.get(),
                          sizeof(uint64_t) * coalesced_buf_size, cudaMemcpyDeviceToHost));

  std::vector<std::pair<MeteoroidDeviation, uint64_t>> zipped(coalesced_buf_size);
  for (size_t i = 0; i < coalesced_buf_size; i++)
  {
    zipped[i] = {host_deviations[i], host_seeds[i]};
  }

  std::sort(zipped.begin(), zipped.end(),
            [](const auto& a, const auto& b) { return a.first.dev < b.first.dev; });

  for (size_t i = 0; i < coalesced_buf_size; i++)
  {
    host_deviations[i] = zipped[i].first;
    host_seeds[i]      = zipped[i].second;
  }

  HANDLE_ERROR(cudaMemcpy(coalesced_deviations_.get(), host_deviations.data(),
                          sizeof(MeteoroidDeviation) * coalesced_buf_size, cudaMemcpyHostToDevice));
  HANDLE_ERROR(cudaMemcpy(coalesced_seeds_.get(), host_seeds.data(),
                          sizeof(uint64_t) * coalesced_buf_size, cudaMemcpyHostToDevice));
}

real MeteoroidsManager::GetActualThresholdDeviation()
{
  real threshold_deviation;
  HANDLE_ERROR(cudaMemcpy(&threshold_deviation, coalesced_deviations_.get() + m_bests_ - 1,
                          sizeof(real), cudaMemcpyDeviceToHost));
  return threshold_deviation;
}

__global__ void RestoreMeteoroidsKernel(VirtualMeteoroid *meteoroids,
                                        const MeteoroidDeviation *deviations,
                                        const uint64_t *seeds,
                                        const size_t size, const TrajectoryPoint ref_point0)
{
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid < size)
  {
    CudaMeteoroidsGenerator generator(seeds[tid], deviations[tid].curand_offset);
    generator.Next(meteoroids[tid], ref_point0);
  }
}

std::vector<std::pair<VirtualMeteoroid, real>>
MeteoroidsManager::GetTopMeteoroids(const TrajectoryPoint &ref_point0)
{
  CudaPtr<VirtualMeteoroid> best_meteoroids_dev;
  HANDLE_ERROR(CudaAlloc(best_meteoroids_dev, m_bests_));

  RestoreMeteoroidsKernel<<<(m_bests_ - 1) / 128 + 1, 128>>>
      (best_meteoroids_dev.get(),
       coalesced_deviations_.get(),
       coalesced_seeds_.get(), m_bests_, ref_point0);
  HANDLE_ERROR(cudaDeviceSynchronize());

  std::vector<VirtualMeteoroid> best_meteoroids(m_bests_);
  HANDLE_ERROR(cudaMemcpy(best_meteoroids.data(), best_meteoroids_dev.get(),
                          sizeof(VirtualMeteoroid) * m_bests_, cudaMemcpyDeviceToHost));

  std::vector<MeteoroidDeviation> best_deviations(m_bests_);
  HANDLE_ERROR(cudaMemcpy(best_deviations.data(), coalesced_deviations_.get(),
                          sizeof(MeteoroidDeviation) * m_bests_, cudaMemcpyDeviceToHost));

  // Something actually got calculated :)
  assert(best_deviations[0].dev != std::numeric_limits<real>::max());

  std::vector<std::pair<VirtualMeteoroid, real>> result(m_bests_);
  for (size_t i = 0; i < m_bests_; i++)
  {
    result[i] = { best_meteoroids[i], best_deviations[i].dev };
  }
  return result;
}