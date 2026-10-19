#pragma once

#include <cstdint>
#include <vector>

namespace jrc {
namespace cuda {

/// Hard limits baked into the JRCLUST kernels. Kept here so callers can check
/// before launching rather than reading garbage out of shared memory.
///   MAXDIM    - features_primary[MAXDIM][CHUNKSIZE] in both kernels
///   NTHREADS  - rho_chunk[NTHREADS][CHUNKSIZE], so blockDim.x may not exceed it
///   CHUNKSIZE - spikes handled per block
constexpr int kMaxFeatures = 45;
constexpr int kMaxThreads = 128;
constexpr int kChunkSize = 16;

/// Mirrors the `consts` int32 vector that computeRho.m / computeDelta.m pack
/// and hand to the kernels: [n1, n12, nFeatures, dnMax, fDcSpk].
struct SiteConstants {
    int n1 = 0;          ///< spikes on this site (primary)
    int n12 = 0;         ///< primary + secondary spikes
    int nFeatures = 0;   ///< feature dimensions per spike (<= kMaxFeatures)
    int dnMax = 0;       ///< max separation in spike-time order
    int fDcSpk = 0;      ///< 1 = per-spike distance cutoff
};

/// Computes rho for one site on the GPU.
///
/// Port of computeRhoSite() in +jrclust/+sort/computeRho.m. The kernel launched
/// is the unmodified jrc_cuda_rho.
///
/// @param siteFeatures  nFeatures x n12, column-major (feature varies fastest),
///                      exactly the layout MATLAB passed as a gpuArray.
/// @param spikeOrder    n12 time-order ranks.
/// @param rhoCut        squared distance cutoff (`distCut2`).
/// @param rhoOut        receives n1 values.
/// @param nThreadsGPU   block width; must be <= kMaxThreads.
void computeRhoSite(const float* siteFeatures,
                    const std::int32_t* spikeOrder,
                    const SiteConstants& consts,
                    float rhoCut,
                    float* rhoOut,
                    int nThreadsGPU = kMaxThreads);

/// Computes delta and the nearest higher-density neighbour for one site.
///
/// Port of computeDeltaSite() in +jrclust/+sort/computeDelta.m, launching the
/// unmodified jrc_cuda_delta kernel.
///
/// @param rhoOrder    n12 descending-rho ranks.
/// @param deltaOut    receives n1 normalised distances.
/// @param neighOut    receives n1 neighbour indices. NOTE: the kernel adds 1 to
///                    produce MATLAB-style indices; this function subtracts it
///                    again so the result is 0-based for C++ callers.
void computeDeltaSite(const float* siteFeatures,
                      const std::int32_t* spikeOrder,
                      const std::int32_t* rhoOrder,
                      const SiteConstants& consts,
                      float distCut2,
                      float* deltaOut,
                      std::int32_t* neighOut,
                      int nThreadsGPU = kMaxThreads);

/// CPU reference for computeRhoSite, used by `jrclust selftest` to prove the
/// kernel path agrees with plain arithmetic. Same semantics, no GPU needed.
void computeRhoSiteCpu(const float* siteFeatures,
                       const std::int32_t* spikeOrder,
                       const SiteConstants& consts,
                       float rhoCut,
                       float* rhoOut);

/// CPU reference for computeDeltaSite. Neighbour indices are 0-based.
void computeDeltaSiteCpu(const float* siteFeatures,
                         const std::int32_t* spikeOrder,
                         const std::int32_t* rhoOrder,
                         const SiteConstants& consts,
                         float distCut2,
                         float* deltaOut,
                         std::int32_t* neighOut);

}  // namespace cuda
}  // namespace jrc
