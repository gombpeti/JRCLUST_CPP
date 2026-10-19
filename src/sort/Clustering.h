#pragma once

#include "core/Config.h"
#include "core/Probe.h"
#include "detect/DetectController.h"

#include <cstdint>
#include <functional>
#include <vector>

namespace jrc {
namespace sort {

/// Density-peak clustering parameters, from the .prm.
struct ClusterSettings {
    double distCut = 2.0;       ///< percentile of pairwise distances (distCut)
    int nClusterIntervals = 4;  ///< nClusterIntervals; sets the time window
    double log10RhoCut = -2.5;
    double log10DeltaCut = 0.6;
    int minClusterSize = 30;
    int maxClustersSite = 20;
    bool useGpu = true;

    static ClusterSettings fromConfig(const Config& config);
};

struct ClusterTimings {
    double rho = 0;
    double delta = 0;
    double assign = 0;
    double total = 0;
};

struct ClusterResult {
    std::vector<float> spikeRho;
    std::vector<float> spikeDelta;
    std::vector<std::int32_t> spikeNeigh;     ///< nearest denser spike, global index
    std::vector<std::int32_t> spikeClusters;  ///< 1-based; 0 = unassigned/noise
    std::vector<std::int32_t> clusterCenters; ///< global spike index per cluster
    std::vector<std::int32_t> clusterSites;   ///< modal site per cluster
    std::vector<std::int32_t> unitCount;
    std::vector<float> siteRhoCut;

    int nClusters = 0;
    ClusterTimings timings;
};

/// Progress callback: (site index, site count, clusters so far).
using ClusterProgressFn = std::function<void(int, int, int)>;

/// Runs density-peak clustering over detected spikes.
///
/// Ports +jrclust/+sort: computeRho, computeDelta, detrendRhoDelta and
/// assignClusters. Per site it gathers that site's spikes (primary, plus
/// secondary when nPeaksFeatures >= 2), launches the original CUDA kernels for
/// rho and delta, then assigns clusters globally from the rho-delta plot.
ClusterResult cluster(const detect::DetectResult& detected,
                      const Probe& probe,
                      const ClusterSettings& settings,
                      const ClusterProgressFn& progress = nullptr);

}  // namespace sort
}  // namespace jrc
