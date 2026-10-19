#pragma once

#include "core/Config.h"
#include "core/Probe.h"
#include "detect/PeakDetector.h"

#include <cstdint>
#include <vector>

namespace jrc {
namespace features {

/// Geometry of the spike windows and the feature vector built from them.
struct FeatureSettings {
    int preSamples = 5;    ///< samples before the peak (|evtWindowSamp[0]|)
    int postSamples = 15;  ///< samples after the peak
    int nSiteDir = 4;      ///< neighbours either side
    int nSitesExcl = 2;    ///< of those, how many are reference-only
    int nPCsPerSite = 3;
    int nPeaksFeatures = 1;  ///< 1 = primary site only, 2 = plus secondary
    int maxPcaSpikes = 10000;

    /// Samples per window: pre + peak + post.
    int windowSamples() const { return preSamples + 1 + postSamples; }

    /// Sites a window is cut on, before local CAR trims it: 2*nSiteDir + 1.
    int nSitesNeighborhood() const { return 2 * nSiteDir + 1; }

    /// Sites kept after local CAR: 1 + 2*nSiteDir - nSitesExcl.
    int nSitesEvt() const { return 1 + 2 * nSiteDir - nSitesExcl; }

    /// Length of one spike's feature vector: nPCsPerSite * nSitesEvt.
    int nFeatures() const { return nPCsPerSite * nSitesEvt(); }

    static FeatureSettings fromConfig(const Config& config);
};

/// Spike windows for one peak assignment.
///
/// Layout is [site][spike][sample] - sample contiguous - so projecting a window
/// onto a principal vector is a dot product over adjacent floats.
struct SpikeWindows {
    std::vector<float> data;
    int nSites = 0;
    int nSpikes = 0;
    int nSamples = 0;

    float* at(int site, int spike)
    {
        return data.data() +
               ((static_cast<std::size_t>(site) * nSpikes) + spike) * nSamples;
    }
    const float* at(int site, int spike) const
    {
        return data.data() +
               ((static_cast<std::size_t>(site) * nSpikes) + spike) * nSamples;
    }
};

/// Feature vectors for every spike: [spike][feature], feature contiguous.
struct SpikeFeatures {
    std::vector<float> data;
    int nSpikes = 0;
    int nFeatures = 0;

    const float* at(int spike) const
    {
        return data.data() + static_cast<std::size_t>(spike) * nFeatures;
    }
    float* at(int spike)
    {
        return data.data() + static_cast<std::size_t>(spike) * nFeatures;
    }
};

/// The three principal vectors, each `nSamples` long.
struct PrincipalVectors {
    std::vector<float> pc;      ///< 3 x nSamples, row-major
    std::vector<double> values; ///< eigenvalues, descending
    int nSamples = 0;
    int nComponents = 0;

    const float* component(int k) const
    {
        return pc.data() + static_cast<std::size_t>(k) * nSamples;
    }
};

/// Cuts windows around each spike on its site neighbourhood.
///
/// Port of extractWindows.m. `spikeSites` gives the centre site per spike, so
/// the same routine serves both the primary and the secondary peak. Windows
/// that run past the ends of the chunk are clamped, as MATLAB does.
SpikeWindows extractWindows(const std::vector<float>& samples,
                            int nSamples,
                            int nSites,
                            const std::vector<std::int64_t>& spikeTimes,
                            const std::vector<std::int32_t>& spikeSites,
                            const Probe& probe,
                            const FeatureSettings& settings);

/// Subtracts the local common-average reference and trims to nSitesEvt.
///
/// Port of localCAR.m: the reference is the mean over the far half of the
/// neighbourhood, and afterwards each window has its own mean removed
/// (meanSubtract.m). Operates in place, shrinking `windows.nSites`.
void applyLocalCar(SpikeWindows& windows, const FeatureSettings& settings);

/// Computes principal vectors from a subsample of primary-site windows.
///
/// Port of spikePrVecs() in spikePCA.m, including the z-scoring of each vector
/// and the sign convention that forces the spike centre negative.
PrincipalVectors computePrincipalVectors(const SpikeWindows& windows,
                                         const FeatureSettings& settings);

/// Projects every window onto the principal vectors.
///
/// Port of projectInterp() without the interpPC refinement. Produces
/// nPCsPerSite * nSitesEvt features per spike, ordered as MATLAB's
/// cat(1, features1, features2, features3): all sites for PC1, then PC2, PC3.
SpikeFeatures projectFeatures(const SpikeWindows& windows,
                              const PrincipalVectors& vectors,
                              const FeatureSettings& settings);

/// Picks the second-strongest site for each spike, by peak-to-peak amplitude.
///
/// Port of findSecondaryPeaks.m. Returns the site index per spike, chosen from
/// the neighbourhood entries 1..nSitesEvt-1 (excluding the primary site).
std::vector<std::int32_t> findSecondarySites(const SpikeWindows& windows,
                                             const std::vector<std::int32_t>& spikeSites,
                                             const Probe& probe,
                                             const FeatureSettings& settings);

}  // namespace features
}  // namespace jrc
