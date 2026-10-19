#pragma once

#include "core/Config.h"
#include "core/Probe.h"
#include "core/Recording.h"
#include "detect/Filter.h"
#include "detect/PeakDetector.h"
#include "detect/SpikeVeto.h"
#include "features/FeatureExtractor.h"
#include "io/ResultWriter.h"

#include <cstdint>
#include <functional>
#include <string>
#include <vector>

namespace jrc {
namespace detect {

/// Everything detection needs, resolved from the .prm in one place.
struct DetectSettings {
    FilterType filterType = FilterType::NDiff;
    int nDiffOrder = 2;
    CarMode carMode = CarMode::Mean;

    DetectParams peaks;
    features::FeatureSettings featureSettings;
    SpikeVetoParams veto;

    /// Extract PCA features alongside detection. Off makes `detect` a pure
    /// spike finder, which is useful for benchmarking the front end.
    bool computeFeatures = true;

    /// Raw-waveform window, from evtWindowRaw. Only used when saving, since
    /// only the curation GUI needs unfiltered traces.
    int rawPreSamples = 10;
    int rawPostSamples = 30;

    /// Keep spike waveforms in memory. Post-processing (mean waveforms,
    /// waveform similarity, quality scores) needs them; plain detection does
    /// not, and they run to gigabytes.
    bool retainWaveforms = false;

    int nSamplesPad = 100;   ///< chunk overlap, from nSamplesPad

    /// Detection chunk size. Purely a performance knob - it must not affect
    /// results, which is why threshold estimation below is defined
    /// independently of it. Measured on a 512-channel 20 kHz recording, 8 s is
    /// the sweet spot: smaller chunks pay per-chunk overhead too often, larger
    /// ones spill further out of cache.
    double chunkSeconds = 8;

    /// Threshold estimation samples fixed windows spread across the recording.
    /// MATLAB averages thresholds over every chunk of the whole file; sampling
    /// a spread is statistically equivalent and avoids a second full pass over
    /// 50 GB. These are deliberately independent of chunkSeconds so that
    /// detection output is reproducible whatever the chunk size.
    double thresholdWindowSeconds = 4.0;
    int thresholdWindows = 12;  ///< 0 or negative means "the whole recording"

    /// Builds settings from a parsed .prm.
    static DetectSettings fromConfig(const Config& config);
};

/// Wall-clock spent in each stage, for finding the next bottleneck.
struct DetectTimings {
    double read = 0;       ///< mapping and converting raw samples
    double transpose = 0;  ///< channel-major to site-major
    double filter = 0;     ///< ndiff
    double car = 0;        ///< common average reference
    double threshold = 0;  ///< RMS estimation
    double peaks = 0;      ///< per-site peak finding and cross-site merge
    double features = 0;   ///< window extraction, local CAR, PCA projection
    double total = 0;
};

/// Detection results for a whole recording.
struct DetectResult {
    std::vector<std::int64_t> spikeTimes;  ///< sample index in the recording
    std::vector<float> spikeAmps;
    std::vector<std::int32_t> spikeSites;  ///< 0-based site index
    std::vector<float> siteThresholds;

    /// How many detected events the per-spike shape tests rejected.
    VetoCounts veto;

    /// Secondary site per spike, when nPeaksFeatures >= 2. Empty otherwise.
    std::vector<std::int32_t> spikeSites2;

    /// Feature vectors, [spike][peak][feature] with feature contiguous.
    /// `nFeatures` is per peak; a spike's row is nPeaks * nFeatures long.
    std::vector<float> spikeFeatures;
    int nFeatures = 0;
    int nPeaks = 0;

    /// The principal vectors features were projected onto - one basis for the
    /// whole recording, so features from different chunks are comparable.
    features::PrincipalVectors principalVectors;

    /// Feature-weighted position of each spike, [x0,y0, x1,y1, ...].
    /// Port of spikePos.m: the centroid of the spike's neighbourhood weighted
    /// by squared PC1 amplitude. Drift-aware mean waveforms need this.
    std::vector<float> spikePositions;

    /// Spike waveforms, kept only when settings.retainWaveforms is set.
    /// Stored exactly as the .jrc files hold them: int16, MATLAB column-major
    /// [nSamples, nWaveformSites, nSpikes], so the sample axis varies fastest.
    std::vector<std::int16_t> filtWaveforms;
    std::vector<std::int16_t> rawWaveforms;
    int filtSamples = 0;
    int rawSamples = 0;
    int nWaveformSites = 0;

    /// Waveform of one spike on one neighbourhood site.
    const std::int16_t* filtWaveform(std::size_t spike, int site) const
    {
        return filtWaveforms.data() +
               (spike * static_cast<std::size_t>(nWaveformSites) + site) * filtSamples;
    }
    const std::int16_t* rawWaveform(std::size_t spike, int site) const
    {
        return rawWaveforms.data() +
               (spike * static_cast<std::size_t>(nWaveformSites) + site) * rawSamples;
    }

    DetectTimings timings;

    std::size_t size() const { return spikeTimes.size(); }

    const float* featuresFor(std::size_t spike, int peak) const
    {
        return spikeFeatures.data() +
               (spike * static_cast<std::size_t>(nPeaks) + peak) * nFeatures;
    }
};

/// Progress callback: (chunk index, chunk count, spikes so far).
using ProgressFn = std::function<void(int, int, std::size_t)>;

/// Appends `part` to `combined`, shifting its spike times by `sampleOffset`.
///
/// Concatenated sorting detects each file on its own - separate thresholds,
/// separate filtering, so a file boundary never lands inside a filter window -
/// then places the spikes on one shared timeline before clustering. Port of the
/// per-recording accumulation in detect.m.
///
/// Both parts must already share a principal-vector basis; see setPrecomputed.
void appendDetectResult(DetectResult& combined, const DetectResult& part,
                        std::int64_t sampleOffset);

/// Runs filter -> CAR -> threshold -> peak detection over a whole recording.
///
/// Two passes, mirroring detectOneRecording.m: the first estimates per-site
/// thresholds, the second detects with those thresholds held fixed so results
/// do not drift between chunks.
class DetectController {
public:
    DetectController(const Config& config, const Probe& probe,
                     const DetectSettings& settings);

    /// Streams waveforms and features to .jrc files as detection proceeds.
    /// Must be called before run(); pass nullptr to disable.
    void setWriter(io::ResultWriter* writer) { writer_ = writer; }

    /// Detects over the whole recording.
    DetectResult run(Recording& recording, const ProgressFn& progress = nullptr);

    /// Detects over samples [begin, end). Used to process a leading slice of a
    /// long recording without reading all of it.
    DetectResult runRange(Recording& recording, std::int64_t begin, std::int64_t end,
                          const ProgressFn& progress = nullptr);

    /// Thresholds only, without detecting. Exposed for diagnostics.
    std::vector<float> estimateThresholds(Recording& recording, std::int64_t begin,
                                          std::int64_t end);

    /// Derives one set of principal vectors for the whole recording.
    ///
    /// MATLAB recomputes these per load, so features from different loads sit
    /// in different bases. Here a single basis is built from spikes sampled
    /// across the whole range, which keeps every spike's features comparable
    /// no matter which chunk found it.
    features::PrincipalVectors estimatePrincipalVectors(
        Recording& recording, std::int64_t begin, std::int64_t end,
        const std::vector<float>& thresholds);

    /// Appends primary-site spike windows sampled across [begin, end) to
    /// `pooled`, for building a principal-vector basis spanning several files.
    void collectPrincipalVectorSamples(Recording& recording, std::int64_t begin,
                                       std::int64_t end,
                                       const std::vector<float>& thresholds,
                                       features::SpikeWindows& pooled);

    /// Supplies thresholds and a principal-vector basis computed elsewhere,
    /// instead of deriving them from the recording being processed.
    ///
    /// Concatenated sorting needs this: every file must project onto the SAME
    /// basis or its features are not comparable, and clustering across files
    /// would be meaningless. Pass nullptr to go back to deriving them.
    void setPrecomputed(const std::vector<float>* thresholds,
                        const features::PrincipalVectors* principalVectors)
    {
        fixedThresholds_ = thresholds;
        fixedPrincipalVectors_ = principalVectors;
    }

private:
    /// Reads [begin,end) for every site, remapping channels to sites, and
    /// applies filter + CAR. Returns the filtered chunk, site-major.
    void loadAndFilter(Recording& recording, std::int64_t begin, std::int64_t end,
                       std::vector<float>& out, int& nSamplesOut);

    /// The read half of loadAndFilter, split out so it can run ahead on a
    /// prefetch thread while the previous chunk is still being processed.
    void loadChunk(Recording& recording, std::int64_t begin, std::int64_t end,
                   std::vector<float>& out, int& nSamplesOut);

    /// The compute half: filter then common-average reference, in place.
    void filterChunk(std::vector<float>& samples, int nSamples);

    /// Number of chunks and their size, given the recording length.
    void planChunks(std::int64_t nSamples, int& nChunks, std::int64_t& chunkSamples) const;

    /// Evenly spaced sampling windows over [begin, end), used by both the
    /// threshold pass and the principal-vector pass so the two agree.
    std::vector<std::int64_t> samplingWindows(std::int64_t begin, std::int64_t end,
                                              std::int64_t& windowSamples) const;

    const Config& config_;
    const Probe& probe_;
    DetectSettings settings_;

    std::vector<float> rawScratch_;
    DetectTimings timings_;
    io::ResultWriter* writer_ = nullptr;

    const std::vector<float>* fixedThresholds_ = nullptr;
    const features::PrincipalVectors* fixedPrincipalVectors_ = nullptr;
};

}  // namespace detect
}  // namespace jrc
