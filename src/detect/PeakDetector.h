#pragma once

#include "core/Probe.h"

#include <cstdint>
#include <vector>

namespace jrc {
namespace detect {

/// Spikes found in one chunk, sorted by time then site.
struct SpikeTable {
    std::vector<std::int64_t> times;  ///< sample index within the chunk
    std::vector<float> amplitudes;    ///< filtered value at the peak (negative)
    std::vector<std::int32_t> sites;  ///< 0-based site index

    std::size_t size() const { return times.size(); }
    void clear();
    void reserve(std::size_t n);
};

/// Knobs taken straight from the .prm.
struct DetectParams {
    double qqFactor = 4.5;
    int minNeighborsDetect = 0;   ///< 0, 1 or 2 neighbours must exceed threshold
    bool detectBipolar = false;   ///< also detect positive peaks
    int refracIntSamples = 5;     ///< refracInt converted to samples
    double evtDetectRad = 200.0;  ///< micrometres
    double spikeThreshMaxRaw = 0; ///< 0 disables the overlarge-spike filter

    // -- shank-wide coincidence removal -----------------------------------
    /// The duplicate resolution above only compares sites within
    /// `evtDetectRad`, which is how JRCLUST does it: one event should survive
    /// on its largest site, and anything further away is assumed to be a
    /// different event. That assumption fails. Measured on a 512-site probe,
    /// one unit held 21,206 pairs of detections on two sites 292 um apart -
    /// far outside the radius - at a median lag of 0.000 ms and 94% identical
    /// amplitude. Two neurons cannot fire on the same sample 21,206 times,
    /// and one neuron cannot be equally strong 292 um away, so both copies
    /// are the same event arriving through crosstalk or a shared reference.
    /// They survived into clustering and produced refractory violations that
    /// then got whole units condemned as noise.
    ///
    /// This pass sweeps each shank as a whole and keeps only the largest of
    /// any group of near-simultaneous detections of comparable amplitude.
    /// Simply widening `evtDetectRad` would not do: that radius also governs
    /// which sites form a neighbourhood for features, and widening it would
    /// merge genuinely distinct nearby neurons.
    bool removeCoincident = true;

    /// Lag within which two detections count as the same event, in samples.
    /// Kept tight, because real synchrony between distinct neurons carries
    /// millisecond jitter while a duplicate sits at zero lag. Set by the
    /// caller from coincidenceMs and the sample rate.
    int coincidenceSamples = 6;

    /// Smaller amplitude over larger, above which the two are "comparable".
    /// This is what protects a large neuron spanning several sites: its field
    /// decays with distance, so the far copy is much smaller and is kept. A
    /// duplicate has the same amplitude on both sites.
    double coincidenceAmpRatio = 0.7;
};

/// Per-site detection threshold, in the same units as the filtered samples.
///
/// Port of computeThreshold.m + estimateRMS.m:
///   thresh = qqFactor * median(|x|) / 0.6745
/// evaluated on at most `maxSamples` evenly spaced samples per site. Ignored
/// sites get a threshold of 0, which suppresses detection there.
std::vector<float> computeThresholds(const std::vector<float>& samples,
                                     int nSamples,
                                     int nSites,
                                     double qqFactor,
                                     const std::vector<std::uint8_t>& ignoredMask,
                                     int maxSamples = 100000);

/// Finds threshold crossings that are local minima, per site, then resolves
/// duplicates across neighbouring sites within the refractory window.
///
/// Combines detectPeaksSite.m, detectPeaks.m and mergePeaks.m. Times are chunk
/// relative; the caller offsets them into the recording.
SpikeTable detectSpikes(const std::vector<float>& samples,
                        int nSamples,
                        int nSites,
                        const std::vector<float>& thresholds,
                        const Probe& probe,
                        const DetectParams& params);

}  // namespace detect
}  // namespace jrc
