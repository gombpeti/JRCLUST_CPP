// DensityPeak.cu - host side of the density-peak clustering kernels.
//
// The kernels themselves are the original JRCLUST sources, copied byte for byte
// from +jrclust/+CUDA/ and #included below. Nothing in them is modified: this
// file only replaces what MATLAB's parallel.gpu.CUDAKernel used to do (allocate
// device memory, marshal arguments, configure the grid, launch).

#include "cuda/DensityPeak.h"
#include "cuda/CudaUtils.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <limits>

// Both kernel files define the same helper macros with identical replacement
// lists. Undefining between includes keeps the sources verbatim and avoids
// leaning on the benign-redefinition rule.
#include "kernels/jrc_cuda_rho.cu"

#undef ABS
#undef MIN
#undef MAX
#undef NTHREADS
#undef MAXDIM
#undef CHUNKSIZE
#undef SINGLE_INF

#include "kernels/jrc_cuda_delta.cu"

namespace jrc {
namespace cuda {

namespace {

/// The kernels index shared memory with blockDim.x, so an oversized block would
/// run off the end of rho_chunk[NTHREADS][...]. Validate before every launch.
void validate(const SiteConstants& c, int nThreadsGPU)
{
    if (c.n1 <= 0 || c.n12 <= 0)
        throw std::runtime_error("DensityPeak: n1 and n12 must be positive");
    if (c.n1 > c.n12)
        throw std::runtime_error("DensityPeak: n1 may not exceed n12");
    if (c.nFeatures <= 0 || c.nFeatures > kMaxFeatures)
        throw std::runtime_error("DensityPeak: nFeatures must be in 1.." +
                                 std::to_string(kMaxFeatures) +
                                 " (kernel shared-memory limit MAXDIM)");
    if (nThreadsGPU <= 0 || nThreadsGPU > kMaxThreads)
        throw std::runtime_error("DensityPeak: nThreadsGPU must be in 1.." +
                                 std::to_string(kMaxThreads) +
                                 " (kernel shared-memory limit NTHREADS)");
}

/// Packs the constants exactly as computeRho.m / computeDelta.m do:
///   consts = int32([n1, n12, nC, dn_max, fDc_spk])
std::array<std::int32_t, 5> packConstants(const SiteConstants& c)
{
    return {c.n1, c.n12, c.nFeatures, c.dnMax, c.fDcSpk};
}

/// Squared euclidean distance between spike a and spike b in feature space.
float featureDist2(const float* f, int nFeatures, int a, int b)
{
    float acc = 0.0f;
    const float* fa = f + static_cast<std::size_t>(a) * nFeatures;
    const float* fb = f + static_cast<std::size_t>(b) * nFeatures;
    for (int i = 0; i < nFeatures; ++i) {
        const float d = fb[i] - fa[i];
        acc += d * d;
    }
    return acc;
}

/// Per-spike cutoff used when fDc_spk == 1: dist_cut2 * ||features||^2.
float spikeCutoff(const float* f, int nFeatures, int i, float distCut2)
{
    float acc = 0.0f;
    const float* fi = f + static_cast<std::size_t>(i) * nFeatures;
    for (int k = 0; k < nFeatures; ++k) acc += fi[k] * fi[k];
    return acc * distCut2;
}

}  // namespace

// ---------------------------------------------------------------------- rho --

void computeRhoSite(const float* siteFeatures,
                    const std::int32_t* spikeOrder,
                    const SiteConstants& consts,
                    float rhoCut,
                    float* rhoOut,
                    int nThreadsGPU)
{
    validate(consts, nThreadsGPU);

    const std::size_t nFeatureValues =
        static_cast<std::size_t>(consts.nFeatures) * consts.n12;

    DeviceBuffer<float> dFeatures(nFeatureValues);
    DeviceBuffer<std::int32_t> dSpikeOrder(consts.n12);
    DeviceBuffer<std::int32_t> dConsts(5);
    DeviceBuffer<float> dRho(consts.n1);

    dFeatures.upload(siteFeatures, nFeatureValues);
    dSpikeOrder.upload(spikeOrder, consts.n12);

    const auto packed = packConstants(consts);
    dConsts.upload(packed.data(), packed.size());
    dRho.zero();

    // MATLAB: rhoCK.GridSize = [ceil(n1/chunkSize), 1]; ThreadBlockSize = [n, 1]
    const dim3 grid((consts.n1 + kChunkSize - 1) / kChunkSize, 1);
    const dim3 block(nThreadsGPU, 1);

    jrc_cuda_rho<<<grid, block>>>(dRho.get(), dFeatures.get(), dSpikeOrder.get(),
                                  dConsts.get(), rhoCut);

    JRC_CUDA_CHECK(cudaGetLastError());
    JRC_CUDA_CHECK(cudaDeviceSynchronize());

    dRho.download(rhoOut, consts.n1);
}

void computeRhoSiteCpu(const float* siteFeatures,
                       const std::int32_t* spikeOrder,
                       const SiteConstants& consts,
                       float rhoCut,
                       float* rhoOut)
{
    for (int i1 = 0; i1 < consts.n1; ++i1) {
        const float cutoff =
            consts.fDcSpk == 1
                ? spikeCutoff(siteFeatures, consts.nFeatures, i1, rhoCut)
                : rhoCut;

        long long nRho = 0;
        long long nComp = 0;

        for (int i12 = 0; i12 < consts.n12; ++i12) {
            const int timeDist = std::abs(spikeOrder[i1] - spikeOrder[i12]);
            if (timeDist > consts.dnMax) continue;

            ++nComp;
            const float d2 = featureDist2(siteFeatures, consts.nFeatures, i1, i12);

            // The kernel uses <= for the global cutoff and < for the per-spike
            // one; mirror that exactly.
            const bool within = (consts.fDcSpk == 0) ? (d2 <= cutoff) : (d2 < cutoff);
            if (within) ++nRho;
        }

        rhoOut[i1] = nComp > 0
                         ? static_cast<float>(static_cast<double>(nRho) /
                                              static_cast<double>(nComp))
                         : 0.0f;
    }
}

// -------------------------------------------------------------------- delta --

void computeDeltaSite(const float* siteFeatures,
                      const std::int32_t* spikeOrder,
                      const std::int32_t* rhoOrder,
                      const SiteConstants& consts,
                      float distCut2,
                      float* deltaOut,
                      std::int32_t* neighOut,
                      int nThreadsGPU)
{
    validate(consts, nThreadsGPU);

    const std::size_t nFeatureValues =
        static_cast<std::size_t>(consts.nFeatures) * consts.n12;

    DeviceBuffer<float> dFeatures(nFeatureValues);
    DeviceBuffer<std::int32_t> dSpikeOrder(consts.n12);
    DeviceBuffer<std::int32_t> dRhoOrder(consts.n12);
    DeviceBuffer<std::int32_t> dConsts(5);
    DeviceBuffer<float> dDelta(consts.n1);
    DeviceBuffer<unsigned int> dNeigh(consts.n1);

    dFeatures.upload(siteFeatures, nFeatureValues);
    dSpikeOrder.upload(spikeOrder, consts.n12);
    dRhoOrder.upload(rhoOrder, consts.n12);

    const auto packed = packConstants(consts);
    dConsts.upload(packed.data(), packed.size());
    dDelta.zero();
    dNeigh.zero();

    const dim3 grid((consts.n1 + kChunkSize - 1) / kChunkSize, 1);
    const dim3 block(nThreadsGPU, 1);

    jrc_cuda_delta<<<grid, block>>>(dDelta.get(), dNeigh.get(), dFeatures.get(),
                                    dSpikeOrder.get(), dRhoOrder.get(),
                                    dConsts.get(), distCut2);

    JRC_CUDA_CHECK(cudaGetLastError());
    JRC_CUDA_CHECK(cudaDeviceSynchronize());

    dDelta.download(deltaOut, consts.n1);

    std::vector<unsigned int> neigh(consts.n1);
    dNeigh.download(neigh.data(), consts.n1);

    // The kernel emits MATLAB 1-based indices; hand C++ callers 0-based ones.
    for (int i = 0; i < consts.n1; ++i)
        neighOut[i] = static_cast<std::int32_t>(neigh[i]) - 1;
}

void computeDeltaSiteCpu(const float* siteFeatures,
                         const std::int32_t* spikeOrder,
                         const std::int32_t* rhoOrder,
                         const SiteConstants& consts,
                         float distCut2,
                         float* deltaOut,
                         std::int32_t* neighOut)
{
    // Matches SINGLE_INF in the kernels.
    constexpr float kSingleInf = 3.402E+38f;

    for (int i1 = 0; i1 < consts.n1; ++i1) {
        float minDist = kSingleInf;
        int minIdx = i1;  // kernel seeds the neighbour with self

        for (int i12 = 0; i12 < consts.n12; ++i12) {
            const bool rhoIsLarger = rhoOrder[i12] < rhoOrder[i1];
            if (!rhoIsLarger) continue;

            const int timeDist = std::abs(spikeOrder[i1] - spikeOrder[i12]);
            if (timeDist > consts.dnMax) continue;

            const float d2 = featureDist2(siteFeatures, consts.nFeatures, i1, i12);
            if (d2 < minDist) {
                minDist = d2;
                minIdx = i12;
            }
        }

        const float cutoff =
            consts.fDcSpk == 0
                ? distCut2
                : spikeCutoff(siteFeatures, consts.nFeatures, i1, distCut2);

        deltaOut[i1] = std::sqrt(std::fabs(minDist) / cutoff);
        neighOut[i1] = minIdx;
    }
}

// --------------------------------------------------------------- device info --

std::vector<DeviceInfo> enumerateDevices()
{
    std::vector<DeviceInfo> devices;

    int count = 0;
    if (cudaGetDeviceCount(&count) != cudaSuccess) return devices;

    for (int i = 0; i < count; ++i) {
        cudaDeviceProp prop{};
        if (cudaGetDeviceProperties(&prop, i) != cudaSuccess) continue;

        DeviceInfo info;
        info.index = i;
        info.name = prop.name;
        info.major = prop.major;
        info.minor = prop.minor;
        info.totalMemory = prop.totalGlobalMem;
        info.multiProcessorCount = prop.multiProcessorCount;
        info.maxThreadsPerBlock = prop.maxThreadsPerBlock;
        info.sharedMemPerBlock = prop.sharedMemPerBlock;
        devices.push_back(std::move(info));
    }

    return devices;
}

}  // namespace cuda
}  // namespace jrc
