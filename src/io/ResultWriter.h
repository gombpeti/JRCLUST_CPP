#pragma once

#include "features/FeatureExtractor.h"

#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

namespace jrc {
namespace io {

/// Streams the large per-spike arrays to JRCLUST's `.jrc` binary files.
///
/// JRCLUST keeps waveforms and features out of the `_res.mat` and in flat
/// binaries beside it, recording only their shape in the .mat. This writer
/// produces byte-compatible files:
///
///   <session>_filt.jrc      int16, [nFiltSamples, nSiteNeighbors, nSpikes]
///   <session>_raw.jrc       int16, [nRawSamples,  nSiteNeighbors, nSpikes]
///   <session>_features.jrc  single, [nFeatures, nPeaks, nSpikes]
///
/// All three are MATLAB column-major, so the sample axis varies fastest.
/// Data is appended chunk by chunk during detection rather than held in RAM;
/// the filtered waveforms alone run to gigabytes on a long recording.
class ResultWriter {
public:
    ~ResultWriter();

    ResultWriter() = default;
    ResultWriter(const ResultWriter&) = delete;
    ResultWriter& operator=(const ResultWriter&) = delete;

    /// Opens the three binaries. `prefix` is the full path without suffix,
    /// e.g. "H:/data/amplifier_denoised".
    void open(const std::string& prefix, bool writeRaw);

    bool isOpen() const { return filtFile_ != nullptr; }

    /// Appends one chunk's filtered waveforms, converting float to int16.
    /// `windows` must still hold the full neighbourhood (pre-local-CAR).
    void appendFilt(const features::SpikeWindows& windows);

    /// Appends one chunk's raw (unfiltered) waveforms.
    void appendRaw(const features::SpikeWindows& windows);

    /// Appends one chunk's features, already laid out [spike][peak][feature].
    void appendFeatures(const float* data, std::size_t count);

    void close();

    // Shapes accumulated so far, in MATLAB order.
    int filtSamples() const { return filtSamples_; }
    int rawSamples() const { return rawSamples_; }
    int nSiteNeighbors() const { return nSiteNeighbors_; }
    std::int64_t nSpikesWritten() const { return nSpikesFilt_; }

private:
    void appendWindows(std::FILE* file, const features::SpikeWindows& windows,
                       int& samplesOut, std::int64_t& counter);

    std::FILE* filtFile_ = nullptr;
    std::FILE* rawFile_ = nullptr;
    std::FILE* featuresFile_ = nullptr;

    int filtSamples_ = 0;
    int rawSamples_ = 0;
    int nSiteNeighbors_ = 0;

    std::int64_t nSpikesFilt_ = 0;
    std::int64_t nSpikesRaw_ = 0;

    std::vector<std::int16_t> scratch_;
};

/// Converts one chunk of windows to int16 in MATLAB's column-major
/// [nSamples, nSites, nSpikes] order and APPENDS them to `out`.
///
/// Shared by the .jrc writer and by in-memory retention so that what
/// post-processing sees is bit-identical to what lands on disk.
void packWindowsInt16(const features::SpikeWindows& windows,
                      std::vector<std::int16_t>& out);

/// Everything the MATLAB importer needs that is not already in a .jrc file.
///
/// Written as `<session>_cpp.bin`: a small header followed by flat arrays.
/// Indices are converted to MATLAB's 1-based convention on the way out, so the
/// importer can use them directly.
struct SidecarData {
    std::int64_t nSpikes = 0;
    int nSites = 0;
    int nClusters = 0;
    int nFeatures = 0;
    int nPeaks = 0;
    int filtSamples = 0;
    int rawSamples = 0;
    int nSiteNeighbors = 0;
    double sampleRate = 0;
    double detectSeconds = 0;
    double sortSeconds = 0;

    const std::vector<std::int64_t>* spikeTimes = nullptr;
    const std::vector<float>* spikeAmps = nullptr;
    const std::vector<std::int32_t>* spikeSites = nullptr;   // 0-based in
    const std::vector<std::int32_t>* spikeSites2 = nullptr;  // 0-based in
    const std::vector<float>* spikeRho = nullptr;
    const std::vector<float>* spikeDelta = nullptr;
    const std::vector<std::int32_t>* spikeNeigh = nullptr;    // 0-based in
    const std::vector<std::int32_t>* spikeClusters = nullptr; // already 1-based
    const std::vector<std::int32_t>* clusterCenters = nullptr;// 0-based in
    const std::vector<std::int32_t>* clusterSites = nullptr;  // 0-based in
    const std::vector<std::int32_t>* unitCount = nullptr;
    const std::vector<float>* siteThresh = nullptr;
    const std::vector<float>* rhoCutSite = nullptr;

    /// Feature-weighted [x,y] per spike, interleaved. MATLAB's drift-aware
    /// averaging depends on these, so shipping the real centroids rather than
    /// letting the importer guess keeps both sides consistent.
    const std::vector<float>* spikePositions = nullptr;

    /// Per-cluster annotation ("noise" / "single" / "multi" / "to_delete" / "").
    /// These land in clusterNotes, so the curation GUI shows them and
    /// Edit > Delete annotated can act on them.
    const std::vector<std::string>* clusterNotes = nullptr;
};

/// Writes `<prefix>_cpp.bin`. Throws std::runtime_error on failure.
void writeSidecar(const std::string& prefix, const SidecarData& data);

}  // namespace io
}  // namespace jrc
