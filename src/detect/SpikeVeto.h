#pragma once

#include "features/FeatureExtractor.h"

#include <cstdint>
#include <vector>

namespace jrc {
namespace detect {

/// Per-spike rejection at detection time.
///
/// Threshold crossing is the only test a conventional sorter applies, so any
/// large transient becomes a "spike". Clustering cannot then refuse it:
/// density-peak assigns every spike by following a chain to its nearest denser
/// neighbour, and that chain always terminates in some cluster. The artifact
/// gets adopted, and no amount of merging removes it.
///
/// These tests ask whether an individual detected event could physically be an
/// action potential, before anything can adopt it. They are deliberately
/// conservative - the cost of discarding a real spike is higher than the cost
/// of passing a bad one to the later stages, which have their own defences.
struct SpikeVetoParams {
    bool enabled = false;

    /// Reject spikes that fall in a zero-filled recording gap. This is
    /// separate from `enabled`: it costs almost nothing, needs no thresholds,
    /// and the samples it rejects are not signal by definition, so it is on
    /// by default.
    bool rejectGaps = true;
    double gapZeroFraction = 0.9;  ///< sites that must read zero
    int gapPadSamples = 0;         ///< 0 = derive from the waveform window

    /// A membrane potential cannot move and return within one sample. The
    /// trough of a real spike spans several; a digital glitch spans one.
    /// Measured as the trough's width at half its depth.
    double minHalfWidthMs = 0.05;

    /// Trough depth relative to its immediate neighbours. A smooth waveform
    /// scores near 1; a one-sample excursion scores far above it.
    double maxSpikiness = 3.0;

    /// Spatial decay. A neuron is a point source, so its extracellular field
    /// falls off with distance: amplitude peaks at one site and decays across
    /// neighbours. An artifact is common-mode and appears at roughly equal
    /// amplitude everywhere. This is the ratio of amplitude on the far half of
    /// the neighbourhood to amplitude on the spike's own site - around 0.1-0.4
    /// for a real spike, near 1.0 for an artifact.
    ///
    /// Site neighbourhoods are ordered by distance, so "far half" is simply
    /// the outer slots.
    double maxDecayRatio = 0.65;

    /// Ignore the decay test when the spike's own site is this quiet, since
    /// the ratio becomes meaningless for a near-zero denominator.
    double minAmplitudeForDecay = 1.0;
};

/// Why spikes were rejected, for reporting.
struct VetoCounts {
    std::int64_t tested = 0;
    std::int64_t tooNarrow = 0;   ///< trough shorter than a real spike
    std::int64_t tooSpiky = 0;    ///< single-sample excursion
    std::int64_t tooFlat = 0;     ///< no spatial decay: common-mode artifact
    std::int64_t rejected = 0;    ///< union of the above

    std::int64_t inGap = 0;       ///< detected inside a zero-filled gap
    std::int64_t gapSamples = 0;  ///< samples blanked as gap

    void add(const VetoCounts& other)
    {
        tested += other.tested;
        tooNarrow += other.tooNarrow;
        tooSpiky += other.tooSpiky;
        tooFlat += other.tooFlat;
        rejected += other.rejected;
        inGap += other.inGap;
        gapSamples += other.gapSamples;
    }
};

/// Marks which spikes survive. `keep` is resized to the spike count.
///
/// The two families of test need DIFFERENT inputs, and using one for both
/// fails badly in opposite directions:
///
///   rawWindows  - for trough width and spikiness. These measure how fast the
///                 membrane moved, so they belong on the voltage trace. Run on
///                 ndiff-filtered data they reject real spikes wholesale,
///                 because a differentiator sharpens every waveform.
///
///   filtWindows - for spatial decay. Raw traces carry low-frequency activity
///                 shared across the whole probe, which makes distant sites
///                 look as large as near ones and hides the falloff. Filtering
///                 and common-average referencing remove it, leaving the
///                 spike's own spatial footprint.
///
/// Both must still hold the full site neighbourhood, before local CAR trims
/// it, so the far sites are available.
VetoCounts vetoSpikes(const features::SpikeWindows& rawWindows,
                      const features::SpikeWindows& filtWindows,
                      const SpikeVetoParams& params, double sampleRate,
                      std::vector<std::uint8_t>& keep);

/// Drops the rejected spikes from a window set, in place.
void compactWindows(features::SpikeWindows& windows,
                    const std::vector<std::uint8_t>& keep);

/// Marks samples that fall in a recording gap, plus a margin either side.
///
/// Acquisition systems fill dropped frames with zeros. Each edge of such a
/// block is a step from full signal to exactly zero on EVERY channel at once -
/// several hundred microvolts, synchronised across the probe - which the
/// detector sees as an enormous spike on every site simultaneously. Those
/// detections then cluster together into units that look like a black bar
/// spanning the whole probe.
///
/// Unlike every other artifact test here, this one needs no statistics: the
/// samples are literally zero. Detecting them from the data rather than from
/// the acquisition log means it also works when no log is available.
///
/// Must be given RAW samples, before filtering - ndiff and common-average
/// referencing turn the zero block's edges into non-zero values.
///
/// @param samples    site-major raw chunk
/// @param zeroFrac   fraction of sites that must be zero to call it a gap
/// @param padSamples margin added either side, to cover every spike whose
///                   waveform would overlap the step
/// @param mask       receives 1 for samples to reject, length nSamples
/// @return number of samples marked
std::int64_t findRecordingGaps(const std::vector<float>& samples, int nSamples,
                               int nSites, double zeroFrac, int padSamples,
                               std::vector<std::uint8_t>& mask);

}  // namespace detect
}  // namespace jrc
