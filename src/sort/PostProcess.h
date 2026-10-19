#pragma once

#include "core/Config.h"
#include "core/Probe.h"
#include "detect/DetectController.h"
#include "sort/Clustering.h"

#include <cstdint>
#include <functional>
#include <string>
#include <vector>

namespace jrc {
namespace sort {

/// Post-processing parameters, from the .prm.
struct PostProcessSettings {
    double maxUnitSim = 0.98;        ///< merge above this waveform correlation
    int nPassesMerge = 10;           ///< repeat merging up to this many times
    std::string autoMergeBy = "pearson";
    double refracInt = 0.25;         ///< ms; also sets the shift search range
    double evtMergeRad = 35.0;       ///< um
    int meanInterpFactor = 1;
    int nClusterIntervals = 4;
    double outlierThresh = 7.5;      ///< MADs, for outlier removal
    double qqFactor = 5.0;

    /// Average only the middle third of the cluster by spike y-position,
    /// realigning each spike onto the cluster's own site neighbourhood
    /// (driftMerge in the .prm). Off averages the spikes nearest the middle of
    /// the recording instead. This materially changes which units merge.
    bool driftMerge = true;

    /// L-ratio degrees of freedom. JRCLUST passes the SPIKE COUNT where the
    /// chi-squared dimensionality belongs, which inflates the statistic by
    /// orders of magnitude; the correct value is the feature count. Set true
    /// only to reproduce MATLAB's numbers.
    bool matlabLRatioDof = false;

    /// Veto waveform-similarity merges that the spike trains contradict.
    ///
    /// Two units that are really one neuron cannot fire within each other's
    /// refractory period, so their cross-correlogram shows a dip at zero lag.
    /// Two genuinely different neurons that happen to share a waveform shape
    /// fire independently, giving a flat correlogram. Comparing observed
    /// coincidences against the count expected from independent trains
    /// separates the two cases; waveform shape alone cannot.
    bool refractoryMergeCheck = false;

    /// Maximum coincidence ratio (observed / expected-if-independent) for a
    /// merge to be allowed. 1.0 is what independent trains give; well below
    /// that is a refractory dip. 0.25 keeps clear same-neuron pairs and
    /// rejects independent ones.
    double maxMergeContamination = 0.25;

    /// Minimum expected coincidences before the veto is allowed to fire.
    ///
    /// For two small or short-lived units the expected count can be a fraction
    /// of one spike, and then the ratio is Poisson noise rather than evidence.
    /// Below this the test has no power, so the merge falls back to the
    /// waveform decision rather than being rejected on a coin flip.
    double minExpectedCoincidences = 5.0;
    double bitScaling = 0.195;
    double sampleRate = 20000.0;

    /// ndiff and sgdiff scale the filtered amplitude, so uV conversion has to
    /// divide it back out. Port of jrclust.utils.bit2uV.
    double filterNorm = 2.0;

    static PostProcessSettings fromConfig(const Config& config);
};

/// Mean waveforms per cluster, in microvolts.
///
/// `global` is [nSamples, nSites, nClusters] column-major - the cluster's mean
/// scattered onto the full probe, zero away from its neighbourhood. `local` is
/// [nSamples, nWaveformSites, nClusters], the neighbourhood alone.
struct MeanWaveforms {
    std::vector<float> global;
    std::vector<float> local;
    std::vector<float> globalRaw;
    std::vector<float> localRaw;

    int nSamples = 0;
    int nSamplesRaw = 0;
    int nSites = 0;
    int nLocalSites = 0;
    int nClusters = 0;

    float* globalAt(int cluster, int site)
    {
        return global.data() +
               ((static_cast<std::size_t>(cluster) * nSites) + site) * nSamples;
    }
    const float* globalAt(int cluster, int site) const
    {
        return global.data() +
               ((static_cast<std::size_t>(cluster) * nSites) + site) * nSamples;
    }
    const float* globalRawAt(int cluster, int site) const
    {
        return globalRaw.data() +
               ((static_cast<std::size_t>(cluster) * nSites) + site) * nSamplesRaw;
    }
};

/// Per-unit quality metrics, matching JRCLUST's names.
struct UnitMetrics {
    std::vector<float> unitPeaks;       ///< |min| of the mean waveform, uV
    std::vector<std::int32_t> unitPeakSites;
    std::vector<float> unitVpp;
    std::vector<float> unitPeaksRaw;
    std::vector<float> unitVppRaw;
    std::vector<float> unitSNR;
    std::vector<std::int32_t> nSitesOverThresh;
    std::vector<float> siteRMS;
    std::vector<double> unitISIRatio;   ///< #(ISI<=2ms) / #(ISI<=20ms)
    std::vector<double> unitIsoDist;
    std::vector<double> unitLRatio;
};

/// Everything post-processing produces.
struct PostProcessResult {
    MeanWaveforms waveforms;
    std::vector<float> waveformSim;  ///< nClusters x nClusters, column-major
    UnitMetrics metrics;
    int nMergePasses = 0;
    int nMergedAway = 0;
    int nMergesVetoed = 0;         ///< blocked by the refractory check
    std::int64_t nOutliersRemoved = 0;
    std::int64_t nRefractoryRemoved = 0;
    double seconds = 0;
};

/// Zero-lag coincidence statistics between two spike trains.
struct CoincidenceTest {
    std::int64_t observed = 0;
    double expected = 0;
    double ratio = 1.0;   ///< observed / expected; ~1 independent, <<1 same unit
    bool hasPower = false; ///< expected is large enough to draw a conclusion
};

/// Counts spike pairs within +/- `refracSamples` and compares them against the
/// count expected from independent trains. Both inputs must be sorted.
CoincidenceTest refractoryCoincidence(const std::vector<std::int64_t>& timesA,
                                      const std::vector<std::int64_t>& timesB,
                                      std::int64_t refracSamples,
                                      std::int64_t recordingSamples,
                                      double minExpected);

using PostProgressFn = std::function<void(const char*, int, int)>;

/// Computes the mean filtered and raw waveform for each cluster.
///
/// Port of computeMeanWaveforms.m taking the non-drift branch of
/// getClusterMean: average the cluster's spikes nearest the middle of the
/// recording (getMidmostSpikes), then subtract each channel's mean.
MeanWaveforms computeMeanWaveforms(const detect::DetectResult& detected,
                                   const ClusterResult& clusters,
                                   const Probe& probe,
                                   const PostProcessSettings& settings);

/// Pairwise mean-waveform similarity.
///
/// Port of computeWaveformSim.m + unitWaveformSim.m: Pearson correlation over
/// the shared site neighbourhood, maximised over sub-refractory time shifts.
/// Returns an nClusters x nClusters symmetric matrix.
std::vector<float> computeWaveformSim(const MeanWaveforms& waveforms,
                                      const ClusterResult& clusters,
                                      const Probe& probe,
                                      const PostProcessSettings& settings);

/// Merges clusters whose similarity exceeds maxUnitSim, repeatedly.
///
/// Port of autoMerge.m / mergeBySim.m. Mutates `clusters` in place: spike
/// labels, centres, sites and counts are all renumbered. Returns the number of
/// clusters merged away.
int autoMerge(const detect::DetectResult& detected, ClusterResult& clusters,
              const Probe& probe, const PostProcessSettings& settings,
              PostProcessResult& out, const PostProgressFn& progress = nullptr);

/// Computes the per-unit quality metrics.
///
/// Port of computeQualityScores.m (both the interface and DensityPeak
/// overrides). Isolation distance and L-ratio use Mahalanobis distances in
/// feature space, as MATLAB does.
UnitMetrics computeQualityScores(const detect::DetectResult& detected,
                                 const ClusterResult& clusters,
                                 const MeanWaveforms& waveforms,
                                 const Probe& probe,
                                 const PostProcessSettings& settings);

/// Renumbers clusters into ascending primary-site order.
/// Port of orderClusters('clusterSites').
void orderClustersBySite(ClusterResult& clusters);

/// Removes individual spikes that do not belong to their unit.
///
/// Port of rmOutlierSpikes.m. Whole-cluster tests cannot catch this: a unit
/// that is 98% real neuron and 2% artifact still has a clean mean waveform, so
/// nothing flags it, yet those 2% are what show up as wild outliers in the
/// amplitude view. This works per spike instead, on the self-Mahalanobis
/// distance of each spike's features from its own cluster - log-transformed to
/// symmetrise the distribution, then scored in MAD units so the threshold is
/// scale-free.
///
/// Excluded spikes are set to cluster 0, matching how MATLAB marks them: they
/// stay in the file but belong to no unit.
std::int64_t removeOutlierSpikes(const detect::DetectResult& detected,
                                 ClusterResult& clusters,
                                 const PostProcessSettings& settings);

/// Drops spikes that violate the refractory period within their own unit,
/// keeping the larger-amplitude spike of each offending pair.
/// Port of rmRefracSpikes.m. Returns the number removed.
std::int64_t removeRefractorySpikes(const detect::DetectResult& detected,
                                    ClusterResult& clusters,
                                    const PostProcessSettings& settings);

/// Runs the standard sequence: mean waveforms, similarity, auto-merge, then
/// quality scores - the same order SortController.m uses.
PostProcessResult runPostProcessing(const detect::DetectResult& detected,
                                    ClusterResult& clusters,
                                    const Probe& probe,
                                    const PostProcessSettings& settings,
                                    bool doAutoMerge,
                                    bool doRemoveRefractory,
                                    bool doRemoveOutliers = false,
                                    const PostProgressFn& progress = nullptr);

}  // namespace sort
}  // namespace jrc
