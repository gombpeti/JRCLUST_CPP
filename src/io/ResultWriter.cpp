#include "io/ResultWriter.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <stdexcept>

namespace jrc {
namespace io {

namespace {

std::FILE* openBinary(const std::string& path)
{
    std::FILE* file = std::fopen(path.c_str(), "wb");
    if (!file) throw std::runtime_error("cannot open for writing: " + path);
    return file;
}

/// Saturating cast to int16, matching MATLAB's int16() on out-of-range values.
inline std::int16_t toInt16(float value)
{
    const float rounded = std::round(value);
    if (rounded <= -32768.0f) return -32768;
    if (rounded >= 32767.0f) return 32767;
    return static_cast<std::int16_t>(rounded);
}

template <typename T>
void writeAll(std::FILE* file, const T* data, std::size_t count, const char* what)
{
    if (count == 0) return;
    if (std::fwrite(data, sizeof(T), count, file) != count)
        throw std::runtime_error(std::string("short write to ") + what);
}

/// Writes a vector as int32, adding `offset` to every element - used to turn
/// 0-based C++ indices into MATLAB's 1-based ones.
void writeInt32(std::FILE* file, const std::vector<std::int32_t>* values,
                std::int32_t offset, std::size_t expected, const char* what)
{
    std::vector<std::int32_t> buffer(expected, 0);
    if (values) {
        const std::size_t n = std::min(expected, values->size());
        for (std::size_t i = 0; i < n; ++i) buffer[i] = (*values)[i] + offset;
    }
    writeAll(file, buffer.data(), buffer.size(), what);
}

void writeFloat(std::FILE* file, const std::vector<float>* values,
                std::size_t expected, const char* what)
{
    std::vector<float> buffer(expected, 0.0f);
    if (values) {
        const std::size_t n = std::min(expected, values->size());
        std::copy_n(values->begin(), n, buffer.begin());
    }
    writeAll(file, buffer.data(), buffer.size(), what);
}

}  // namespace

void packWindowsInt16(const features::SpikeWindows& windows,
                      std::vector<std::int16_t>& out)
{
    if (windows.nSpikes == 0) return;

    const int nSamples = windows.nSamples;
    const int nSites = windows.nSites;
    const std::size_t perSpike = static_cast<std::size_t>(nSamples) * nSites;

    const std::size_t base = out.size();
    out.resize(base + perSpike * windows.nSpikes);

    // MATLAB order is [sample, site, spike] column-major, so for each spike the
    // sample axis varies fastest, then site. Our buffer is [site][spike][sample].
    for (int spike = 0; spike < windows.nSpikes; ++spike) {
        std::int16_t* dst = out.data() + base + static_cast<std::size_t>(spike) * perSpike;
        for (int site = 0; site < nSites; ++site) {
            const float* in = windows.at(site, spike);
            std::int16_t* row = dst + static_cast<std::size_t>(site) * nSamples;
            for (int s = 0; s < nSamples; ++s) row[s] = toInt16(in[s]);
        }
    }
}

ResultWriter::~ResultWriter() { close(); }

void ResultWriter::open(const std::string& prefix, bool writeRaw)
{
    close();
    filtFile_ = openBinary(prefix + "_filt.jrc");
    featuresFile_ = openBinary(prefix + "_features.jrc");
    if (writeRaw) rawFile_ = openBinary(prefix + "_raw.jrc");
}

void ResultWriter::appendWindows(std::FILE* file,
                                 const features::SpikeWindows& windows,
                                 int& samplesOut, std::int64_t& counter)
{
    if (!file || windows.nSpikes == 0) return;

    samplesOut = windows.nSamples;
    nSiteNeighbors_ = windows.nSites;

    scratch_.clear();
    packWindowsInt16(windows, scratch_);

    writeAll(file, scratch_.data(), scratch_.size(), "waveform file");
    counter += windows.nSpikes;
}

void ResultWriter::appendFilt(const features::SpikeWindows& windows)
{
    appendWindows(filtFile_, windows, filtSamples_, nSpikesFilt_);
}

void ResultWriter::appendRaw(const features::SpikeWindows& windows)
{
    appendWindows(rawFile_, windows, rawSamples_, nSpikesRaw_);
}

void ResultWriter::appendFeatures(const float* data, std::size_t count)
{
    if (!featuresFile_ || count == 0) return;
    writeAll(featuresFile_, data, count, "features file");
}

void ResultWriter::close()
{
    if (filtFile_) std::fclose(filtFile_);
    if (rawFile_) std::fclose(rawFile_);
    if (featuresFile_) std::fclose(featuresFile_);
    filtFile_ = rawFile_ = featuresFile_ = nullptr;
}

void writeSidecar(const std::string& prefix, const SidecarData& data)
{
    const std::string path = prefix + "_cpp.bin";
    std::FILE* file = openBinary(path);

    try {
        // Header: magic, then the shapes the importer needs to reshape the
        // .jrc files. Version byte lets the importer refuse a mismatch.
        const char magic[8] = {'J', 'R', 'C', 'P', 'P', '0', '0', '4'};
        writeAll(file, magic, 8, "sidecar");

        const std::int64_t nSpikes = data.nSpikes;
        writeAll(file, &nSpikes, 1, "sidecar");

        const std::int32_t header[] = {
            data.nSites,     data.nClusters,       data.nFeatures,
            data.nPeaks,     data.filtSamples,     data.rawSamples,
            data.nSiteNeighbors};
        writeAll(file, header, sizeof(header) / sizeof(header[0]), "sidecar");

        const double reals[] = {data.sampleRate, data.detectSeconds, data.sortSeconds};
        writeAll(file, reals, 3, "sidecar");

        const std::size_t n = static_cast<std::size_t>(nSpikes);

        // spikeTimes: MATLAB indexes samples from 1.
        {
            std::vector<std::int32_t> times(n, 0);
            if (data.spikeTimes) {
                const std::size_t m = std::min(n, data.spikeTimes->size());
                for (std::size_t i = 0; i < m; ++i)
                    times[i] = static_cast<std::int32_t>((*data.spikeTimes)[i] + 1);
            }
            writeAll(file, times.data(), times.size(), "sidecar");
        }

        // spikeAmps as int16, the type JRCLUST stores.
        {
            std::vector<std::int16_t> amps(n, 0);
            if (data.spikeAmps) {
                const std::size_t m = std::min(n, data.spikeAmps->size());
                for (std::size_t i = 0; i < m; ++i) amps[i] = toInt16((*data.spikeAmps)[i]);
            }
            writeAll(file, amps.data(), amps.size(), "sidecar");
        }

        writeInt32(file, data.spikeSites, 1, n, "sidecar");
        writeInt32(file, data.spikeSites2, 1, n, "sidecar");
        writeFloat(file, data.spikeRho, n, "sidecar");
        writeFloat(file, data.spikeDelta, n, "sidecar");
        writeInt32(file, data.spikeNeigh, 1, n, "sidecar");
        writeInt32(file, data.spikeClusters, 0, n, "sidecar");  // already 1-based

        const std::size_t c = static_cast<std::size_t>(data.nClusters);
        writeInt32(file, data.clusterCenters, 1, c, "sidecar");
        writeInt32(file, data.clusterSites, 1, c, "sidecar");
        writeInt32(file, data.unitCount, 0, c, "sidecar");

        const std::size_t s = static_cast<std::size_t>(data.nSites);
        writeFloat(file, data.siteThresh, s, "sidecar");
        writeFloat(file, data.rhoCutSite, s, "sidecar");

        // spikePositions, [x,y] interleaved per spike.
        writeFloat(file, data.spikePositions, n * 2, "sidecar");

        // Cluster annotations, length-prefixed so MATLAB can read them back as
        // a cell array of char.
        for (std::size_t k = 0; k < c; ++k) {
            const std::string note =
                (data.clusterNotes && k < data.clusterNotes->size())
                    ? (*data.clusterNotes)[k]
                    : std::string();
            const std::int32_t length = static_cast<std::int32_t>(note.size());
            writeAll(file, &length, 1, "sidecar");
            if (length > 0) writeAll(file, note.data(), note.size(), "sidecar");
        }
    } catch (...) {
        std::fclose(file);
        throw;
    }

    std::fclose(file);
}

}  // namespace io
}  // namespace jrc
