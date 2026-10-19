#pragma once

#include "core/Probe.h"
#include "detect/DetectController.h"
#include "sort/Clustering.h"
#include "sort/PostProcess.h"

#include <cstdint>
#include <string>
#include <vector>

namespace jrc {
namespace sort {

/// Unit quality, matching JRCLUST's annotations and their keyboard shortcuts:
///   0 noise  - the waveform is not an action potential
///   1 single - clean, well isolated
///   2 multi  - a real spike shape, but contaminated
enum class UnitQuality { Noise = 0, Single = 1, Multi = 2 };

const char* qualityName(UnitQuality quality);  ///< "noise" / "single" / "multi"

/// Thresholds for classifying units.
///
/// The design principle is that SHAPE is amplitude-invariant. A large unit
/// sitting right against the electrode still has the stereotyped shape of an
/// action potential; a digital glitch does not, however large or small it is.
/// So the decisive noise tests are on shape, and amplitude only ever
/// contributes supporting evidence. A plain microvolt ceiling would eventually
/// discard a real close-by unit while missing every low-amplitude artifact.
struct QualityCriteria {
    // -- noise: decisive, essentially never true of a neuron -------------
    double maxSNR = 50.0;              ///< peak / noise RMS; real units are 3-20
    double glitchSpikiness = 3.0;      ///< trough vs its immediate neighbours
    double glitchHalfWidthMs = 0.075;  ///< with the above: a one-sample excursion

    /// Amplitude outliers are judged against this recording's own units rather
    /// than an absolute microvolt figure, because the sensible ceiling depends
    /// on probe, prep and filter. A unit must be this many times the MEDIAN on
    /// both SNR and Vpp before amplitude alone condemns it - medians are used
    /// instead of means so a handful of artifacts cannot drag the reference up.
    double snrOutlierFactor = 5.0;
    double vppOutlierFactor = 5.0;

    /// Refractory violations so frequent that no neuron could produce them.
    /// A clean unit sits near 0.02 and even a strong burster stays under ~0.2.
    /// Measured on real contaminated units: 0.48, 0.72, 0.99, 1.00. Decisive
    /// on its own - a spike train this self-overlapping is not one cell,
    /// whatever its mean waveform looks like.
    double decisiveISIRatio = 0.45;

    /// Zero-lag autocorrelogram: near-lag pair count against mid-lag, scaled
    /// so that 1.0 is "no refractory period at all".
    ///
    /// CALIBRATED, not derived. The theoretical value of 1.0 is wrong in
    /// practice and was deleting real neurons: measured over a 2.6 h, 512-site
    /// recording, units that are clean by every independent measure (ISI ratio
    /// <= 0.2, SNR >= 4, >= 5000 spikes, n = 239) reach 1.59 at p95 and 2.05 at
    /// p99, while units that are unarguably contaminated (ISI ratio >= 0.6,
    /// n = 57) start at 3.19 at p1. The two populations are cleanly separated
    /// with an empty gap between 2.05 and 3.19, so the threshold belongs in
    /// that gap rather than on top of the good units.
    ///
    /// At 1.0 this test condemned 342 of 1200 clusters and 23% of all spikes,
    /// among them a 197,474-spike unit firing steadily at 21 Hz for the whole
    /// recording that scored 1.074. At 3.0 it condemns 173 clusters and 4.6% of
    /// spikes, keeping every unit in the clean population and still catching
    /// 99% of the contaminated one.
    double maxZeroLagRatio = 3.0;

    /// Sites carrying the unit, as an ABSOLUTE count. A point source reaches a
    /// handful of neighbours; an event on dozens of sites is common-mode.
    /// Counted in sites rather than as a probe fraction, because a fraction
    /// scales with probe size and lets a 50-site artifact pass on 512 sites.
    int maxSpreadSites = 20;

    // -- noise: the same event on a whole shank at one amplitude ----------
    /// A neuron is a point source. Its field peaks at one site and falls off,
    /// so it reaches two to four neighbouring contacts - being on several
    /// channels at once is normal and is NOT evidence of noise. What no neuron
    /// can do is appear on every contact of a shank at the SAME amplitude:
    /// that is the whole shank moving together, which means common-mode
    /// pickup, a movement artifact or a gap edge.
    ///
    /// Measured per shank rather than across the probe, because a 64-site
    /// shank fully covered is only 12.5% of a 512-site probe and slips past
    /// any probe-wide count or fraction - `maxSpreadSites` above never fired
    /// once on a real 512-site recording.
    ///
    /// "Same amplitude" is the operative part. Counting sites above a low
    /// fraction of the peak would catch a large neuron whose tail is still
    /// visible far away; the cut is set at half the peak so that only sites
    /// carrying a comparable deflection count.
    double shankUniformFraction = 0.5;  ///< >= this much of the peak = "same amplitude"
    double maxShankCoverage = 0.5;      ///< above this share of a shank = common mode
    int minShankSitesForTest = 8;       ///< skip shanks too small to judge

    /// Fraction of a unit's spikes allowed to sit in windows where its own
    /// rate far exceeds its own baseline. Catches artifacts confined to one
    /// unit, which a population-rate test cannot see.
    double maxBurstFraction = 0.30;
    double burstBinSeconds = 10.0;
    double burstMadThresh = 10.0;
    double burstMinRatio = 8.0;

    // -- noise: supporting, two or more together are enough --------------
    /// Deliberately set well above the population p95: a bursting neuron can
    /// carry a high ISI ratio legitimately, so only extreme values count.
    double noiseISIRatio = 0.50;
    double maxSpreadFraction = 0.5;    ///< on over half the probe = common mode
    double maxVpp = 0.0;               ///< uV; 0 disables the absolute ceiling

    // NOTE: spatial spread of one site and a weak positive rebound were tried
    // as noise evidence and dropped. Measured across a real recording the
    // median unit sits on 2 sites and the 5th percentile of rebound is 0, so
    // both fire on ordinary neurons and produced false positives.

    // -- single vs multi unit --------------------------------------------
    double singleMinSNR = 5.0;
    double singleMaxISIRatio = 0.05;
    double singleMinIsoDist = 20.0;    ///< NaN counts as unknown, not failing

    /// Sites counted as carrying the unit, as a fraction of its peak.
    double spreadThreshold = 0.25;

    /// Annotation written for noise units.
    ///
    /// 'noise' is the semantically correct label and is what the GUI's own
    /// Noise annotation writes. JRCLUST's Edit > Delete annotated, however,
    /// removes units annotated exactly 'to_delete' - so choose that instead to
    /// arm one-click bulk removal.
    std::string noiseLabel = "noise";

    /// Whether to annotate the units that pass as 'single' / 'multi' too.
    ///
    /// On by default, so every unit carries a 0/1/2 verdict. The trade-off is
    /// that the GUI's "Display annotated units only" filters on *any*
    /// annotation, so with every unit labelled it matches all of them; turn
    /// this off to make that filter show only the units the classifier doubts.
    bool labelGoodUnits = true;
};

/// Interpretable shape measurements for one unit, from its mean raw waveform.
struct UnitShape {
    double troughDepth = 0;      ///< uV, positive magnitude
    double vpp = 0;              ///< uV
    double halfWidthMs = 0;      ///< trough width at half depth
    double troughToPeakMs = 0;   ///< trough to the following positive peak
    double spikiness = 0;        ///< |trough| / mean(|neighbours|)
    double recoveryRatio = 0;    ///< rebound height / trough depth
    int spatialSpread = 0;       ///< sites carrying >= spreadThreshold of peak
    int peakSite = 0;

    // Spread measured within the peak site's own shank. A neuron reaches a few
    // neighbours here; an artifact covers the shank.
    int shankSites = 0;          ///< active (non-ignored) sites on that shank
    int shankAbove = 0;          ///< of those, at >= shankUniformFraction of peak
    double shankCoverage = 0;    ///< shankAbove / shankSites
};

/// Counts how much of the peak site's shank carries the unit at comparable
/// amplitude. Fills the shank* fields of `shape`.
void measureShankSpread(const MeanWaveforms& waveforms, int cluster,
                        const Probe& probe, double uniformFraction,
                        UnitShape& shape);

struct UnitVerdict {
    UnitQuality quality = UnitQuality::Multi;
    UnitShape shape;
    std::vector<std::string> reasons;  ///< why, in plain words
};

struct QualityClassification {
    std::vector<UnitVerdict> verdicts;
    int nNoise = 0;
    int nSingle = 0;
    int nMulti = 0;
};

/// Measures the shape of one unit's mean raw waveform.
UnitShape measureShape(const MeanWaveforms& waveforms, int cluster, int peakSite,
                       double sampleRate, double spreadThreshold);

/// Classifies every unit into noise / single / multi.
///
/// Runs after merging, so it judges final units rather than fragments. Nothing
/// is deleted - units are labelled and the caller decides.
QualityClassification classifyUnits(const detect::DetectResult& detected,
                                    const ClusterResult& clusters,
                                    const MeanWaveforms& waveforms,
                                    const UnitMetrics& metrics,
                                    const Probe& probe,
                                    const PostProcessSettings& settings,
                                    const QualityCriteria& criteria);

/// Per-cluster annotation strings for MATLAB, indexed by cluster - 1.
std::vector<std::string> annotationStrings(const QualityClassification& classification,
                                           const QualityCriteria& criteria);

/// Outcome of splitting bimodal clusters.
struct SplitReport {
    int nClustersTested = 0;
    int nSplit = 0;                  ///< clusters divided in two
    std::int64_t nSpikesMovedOut = 0; ///< spikes that ended up in the new units
    std::vector<int> splitParents;   ///< original cluster ids that were split
};

/// Settings for the bimodal split.
struct SplitCriteria {
    /// How separated the two halves must be before a split is accepted,
    /// measured as the distance between their centres in units of their
    /// pooled spread. 2.0 means the modes are two standard deviations apart -
    /// visibly distinct, not a smooth tail being cut arbitrarily.
    double minSeparation = 2.0;

    /// Smallest half worth creating; below this the "second mode" is just a
    /// handful of stragglers.
    int minHalfSize = 15;

    /// Only split when the two halves disagree about being a neuron. Both
    /// halves looking real means the cluster holds two units, which is a
    /// merging question rather than a denoising one, and is left alone.
    bool requireOneHalfNoise = true;

    int maxIterations = 50;  ///< k-means iterations
};

/// Splits clusters that contain two separable populations, then relabels.
///
/// This automates what curation does by hand: divide a suspect cluster in two,
/// see that one half is clean spikes and the other is artifact, keep the first
/// and discard the second. Whole-cluster tests miss these entirely, because a
/// cluster that is mostly neuron has a mostly-neuronal mean waveform however
/// bad the contaminating minority is.
///
/// Mutates `clusters`: the noise half of each split becomes a new cluster at
/// the end, which `classifyUnits` then labels like any other.
SplitReport splitBimodalClusters(const detect::DetectResult& detected,
                                 ClusterResult& clusters,
                                 const MeanWaveforms& waveforms,
                                 const UnitMetrics& metrics,
                                 const Probe& probe,
                                 const PostProcessSettings& settings,
                                 const QualityCriteria& criteria,
                                 const SplitCriteria& splitCriteria);

}  // namespace sort
}  // namespace jrc
