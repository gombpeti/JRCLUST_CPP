#include "features/FeatureExtractor.h"

#include "core/Parallel.h"
#include "core/SymEigen.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <numeric>
#include <stdexcept>

namespace jrc {
namespace features {

FeatureSettings FeatureSettings::fromConfig(const Config& config)
{
    FeatureSettings settings;

    const double sampleRate = config.sampleRate();
    const auto window = config.getNumbers("evtWindow");

    if (window.size() >= 2) {
        // evtWindow is in milliseconds, e.g. [-0.25, 0.75].
        const double preMs = window[0];
        const double postMs = window[1];
        settings.preSamples =
            static_cast<int>(std::round(std::fabs(preMs) * sampleRate / 1000.0));
        settings.postSamples =
            static_cast<int>(std::round(postMs * sampleRate / 1000.0));
    }

    settings.nSiteDir = config.getInt("nSiteDir", 4);
    settings.nSitesExcl = config.getInt("nSitesExcl", 2);
    settings.nPCsPerSite = config.getInt("nPCsPerSite", 3);
    settings.nPeaksFeatures = config.getInt("nPeaksFeatures", 1);

    if (settings.nPCsPerSite < 1 || settings.nPCsPerSite > 3)
        throw std::runtime_error("nPCsPerSite must be 1, 2 or 3");
    if (settings.nSitesEvt() < 1)
        throw std::runtime_error("nSitesEvt resolved to < 1; check nSiteDir/nSitesExcl");

    return settings;
}

SpikeWindows extractWindows(const std::vector<float>& samples,
                            int nSamples,
                            int nSites,
                            const std::vector<std::int64_t>& spikeTimes,
                            const std::vector<std::int32_t>& spikeSites,
                            const Probe& probe,
                            const FeatureSettings& settings)
{
    const int nSpikes = static_cast<int>(spikeTimes.size());
    const int nWin = settings.windowSamples();
    const int nNeighborhood =
        std::min(settings.nSitesNeighborhood(), probe.nSiteNeighbors());

    SpikeWindows windows;
    windows.nSites = nNeighborhood;
    windows.nSpikes = nSpikes;
    windows.nSamples = nWin;
    if (nSpikes == 0) return windows;

    windows.data.assign(static_cast<std::size_t>(nNeighborhood) * nSpikes * nWin, 0.0f);

    const int pre = settings.preSamples;

    parallelFor(static_cast<std::size_t>(nSpikes), [&](std::size_t spikeIndex) {
        const int spike = static_cast<int>(spikeIndex);
        const std::int64_t peak = spikeTimes[spike];
        const int centreSite = spikeSites[spike];
        const int* neighbors = probe.siteNeighbors(centreSite);

        for (int k = 0; k < nNeighborhood; ++k) {
            const int site = neighbors[k];
            const float* trace =
                samples.data() + static_cast<std::size_t>(site) * nSamples;
            float* out = windows.at(k, spike);

            for (int s = 0; s < nWin; ++s) {
                // Clamp to the chunk, matching min(max(timeRanges,1),nSamples)
                // in extractWindows.m.
                std::int64_t t = peak - pre + s;
                if (t < 0) t = 0;
                if (t >= nSamples) t = nSamples - 1;
                out[s] = trace[t];
            }
        }
    });

    return windows;
}

void applyLocalCar(SpikeWindows& windows, const FeatureSettings& settings)
{
    const int nSitesEvt = std::min(settings.nSitesEvt(), windows.nSites);
    const int nWin = windows.nSamples;
    const int nSpikes = windows.nSpikes;
    if (nSpikes == 0 || nSitesEvt < 1) return;

    if (nSitesEvt == 1) {
        // localCAR.m returns early with no reference for a single site.
        windows.nSites = 1;
        windows.data.resize(static_cast<std::size_t>(nSpikes) * nWin);
        return;
    }

    // Reference is the mean over the far half of the neighbourhood: MATLAB's
    // farHalf = ceil(nSites/2):nSites, 1-based, so 0-based this starts at
    // ceil(nSites/2) - 1.
    const int farBegin = (windows.nSites + 1) / 2 - 1;
    const int nFar = windows.nSites - farBegin;

    std::vector<float> reference(static_cast<std::size_t>(nSpikes) * nWin, 0.0f);

    parallelFor(static_cast<std::size_t>(nSpikes), [&](std::size_t spikeIndex) {
        const int spike = static_cast<int>(spikeIndex);
        float* ref = reference.data() + static_cast<std::size_t>(spike) * nWin;

        for (int k = farBegin; k < windows.nSites; ++k) {
            const float* row = windows.at(k, spike);
            for (int s = 0; s < nWin; ++s) ref[s] += row[s];
        }
        const float scale = 1.0f / static_cast<float>(nFar);
        for (int s = 0; s < nWin; ++s) ref[s] *= scale;
    });

    // Trim to the near sites, then subtract the reference and each window's own
    // mean. Compacting in place is safe because nSitesEvt <= nSites, so the
    // destination for site k never runs past the source.
    for (int k = 0; k < nSitesEvt; ++k) {
        float* dst = windows.data.data() +
                     (static_cast<std::size_t>(k) * nSpikes) * nWin;
        const float* src = windows.at(k, 0);
        if (dst != src)
            std::memmove(dst, src, static_cast<std::size_t>(nSpikes) * nWin * sizeof(float));
    }

    windows.nSites = nSitesEvt;
    windows.data.resize(static_cast<std::size_t>(nSitesEvt) * nSpikes * nWin);

    parallelFor(static_cast<std::size_t>(nSitesEvt), [&](std::size_t siteIndex) {
        const int k = static_cast<int>(siteIndex);
        for (int spike = 0; spike < nSpikes; ++spike) {
            float* row = windows.at(k, spike);
            const float* ref =
                reference.data() + static_cast<std::size_t>(spike) * nWin;

            double sum = 0.0;
            for (int s = 0; s < nWin; ++s) {
                row[s] -= ref[s];
                sum += row[s];
            }
            // meanSubtract.m removes the mean along the sample axis.
            const float mean = static_cast<float>(sum / nWin);
            for (int s = 0; s < nWin; ++s) row[s] -= mean;
        }
    });
}

PrincipalVectors computePrincipalVectors(const SpikeWindows& windows,
                                         const FeatureSettings& settings)
{
    const int nWin = windows.nSamples;
    const int nSpikes = windows.nSpikes;

    PrincipalVectors result;
    result.nSamples = nWin;
    result.nComponents = 3;
    result.pc.assign(static_cast<std::size_t>(3) * nWin, 0.0f);
    if (nSpikes == 0) return result;

    // Subsample the primary site's windows, evenly rather than at random:
    // spikePrVecs.m uses randsample(), which makes its vectors differ between
    // runs on identical input.
    const int stride = std::max(1, nSpikes / std::max(1, settings.maxPcaSpikes));
    int nUsed = 0;

    std::vector<double> cov(static_cast<std::size_t>(nWin) * nWin, 0.0);

    for (int spike = 0; spike < nSpikes; spike += stride) {
        const float* w = windows.at(0, spike);
        for (int i = 0; i < nWin; ++i) {
            const double wi = w[i];
            double* row = cov.data() + static_cast<std::size_t>(i) * nWin;
            for (int j = 0; j <= i; ++j) row[j] += wi * w[j];
        }
        ++nUsed;
    }

    if (nUsed < 2) return result;

    // covMat = (X*X')/(n-1)
    const double scale = 1.0 / static_cast<double>(nUsed - 1);
    for (auto& value : cov) value *= scale;
    for (int i = 0; i < nWin; ++i)
        for (int j = i + 1; j < nWin; ++j)
            cov[static_cast<std::size_t>(i) * nWin + j] =
                cov[static_cast<std::size_t>(j) * nWin + i];

    std::vector<double> values;
    std::vector<double> vectors;
    symmetricEigen(cov, nWin, values, vectors);

    result.values = values;

    // The spike is at index preSamples within the window; spikePrVecs.m flips
    // each vector so that its value there is negative.
    const int centre = std::min(settings.preSamples, nWin - 1);

    for (int k = 0; k < 3 && k < nWin; ++k) {
        // zscore each vector: subtract mean, divide by sample std (N-1).
        double sum = 0.0;
        for (int i = 0; i < nWin; ++i)
            sum += vectors[static_cast<std::size_t>(i) * nWin + k];
        const double mean = sum / nWin;

        double variance = 0.0;
        for (int i = 0; i < nWin; ++i) {
            const double d = vectors[static_cast<std::size_t>(i) * nWin + k] - mean;
            variance += d * d;
        }
        double sd = std::sqrt(variance / std::max(1, nWin - 1));
        if (sd == 0.0) sd = 1.0;  // no dispersion; do not divide by zero

        float* out = result.pc.data() + static_cast<std::size_t>(k) * nWin;
        for (int i = 0; i < nWin; ++i)
            out[i] = static_cast<float>(
                (vectors[static_cast<std::size_t>(i) * nWin + k] - mean) / sd);

        const float sign = (out[centre] < 0.0f) ? 1.0f : -1.0f;
        if (sign < 0.0f)
            for (int i = 0; i < nWin; ++i) out[i] = -out[i];
    }

    return result;
}

SpikeFeatures projectFeatures(const SpikeWindows& windows,
                              const PrincipalVectors& vectors,
                              const FeatureSettings& settings)
{
    const int nWin = windows.nSamples;
    const int nSpikes = windows.nSpikes;
    const int nSites = windows.nSites;
    const int nPCs = settings.nPCsPerSite;

    SpikeFeatures out;
    out.nSpikes = nSpikes;
    out.nFeatures = nPCs * nSites;
    if (nSpikes == 0) return out;

    out.data.assign(static_cast<std::size_t>(nSpikes) * out.nFeatures, 0.0f);

    parallelFor(static_cast<std::size_t>(nSpikes), [&](std::size_t spikeIndex) {
        const int spike = static_cast<int>(spikeIndex);
        float* feature = out.at(spike);

        for (int pc = 0; pc < nPCs; ++pc) {
            const float* vector = vectors.component(pc);
            for (int site = 0; site < nSites; ++site) {
                const float* w = windows.at(site, spike);

                float acc = 0.0f;
                for (int s = 0; s < nWin; ++s) acc += vector[s] * w[s];

                // MATLAB's cat(1, f1, f2, f3): all sites for PC1, then PC2, PC3.
                feature[pc * nSites + site] = acc;
            }
        }
    });

    return out;
}

std::vector<std::int32_t> findSecondarySites(const SpikeWindows& windows,
                                             const std::vector<std::int32_t>& spikeSites,
                                             const Probe& probe,
                                             const FeatureSettings& settings)
{
    const int nSpikes = windows.nSpikes;
    std::vector<std::int32_t> result(nSpikes, 0);
    if (nSpikes == 0) return result;

    // findSecondaryPeaks.m searches neighbourhood entries 2..nSitesEvt, i.e.
    // everything except the spike's own site.
    const int last = std::min(settings.nSitesEvt(), windows.nSites);

    parallelFor(static_cast<std::size_t>(nSpikes), [&](std::size_t spikeIndex) {
        const int spike = static_cast<int>(spikeIndex);
        const int* neighbors = probe.siteNeighbors(spikeSites[spike]);

        float bestVpp = -1.0f;
        int bestSlot = 1;

        for (int k = 1; k < last; ++k) {
            const float* w = windows.at(k, spike);
            float lo = w[0], hi = w[0];
            for (int s = 1; s < windows.nSamples; ++s) {
                lo = std::min(lo, w[s]);
                hi = std::max(hi, w[s]);
            }
            const float vpp = hi - lo;
            if (vpp > bestVpp) {
                bestVpp = vpp;
                bestSlot = k;
            }
        }

        result[spike] = static_cast<std::int32_t>(neighbors[bestSlot]);
    });

    return result;
}

}  // namespace features
}  // namespace jrc
