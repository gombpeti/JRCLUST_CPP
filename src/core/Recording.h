#pragma once

#include "core/Config.h"
#include "core/MemoryMappedFile.h"

#include <cstdint>
#include <string>
#include <vector>

namespace jrc {

/// A raw extracellular recording held in a flat binary file.
///
/// Port of +jrclust/+detect/@Recording. Samples are stored channel-interleaved:
/// element (channel c, sample s) sits at `headerOffset + (c + s*nChans)*bytes`,
/// which is what MATLAB's column-major memmapfile with dshape [nChans nSamples]
/// produced.
class Recording {
public:
    /// Maps `path` and derives the sample count from the file size.
    /// Throws std::runtime_error if the size is not a whole number of samples.
    void open(const std::string& path, const Config& config);
    void close();

    bool isOpen() const { return file_.isOpen(); }

    int nChans() const { return nChans_; }
    std::int64_t nSamples() const { return nSamples_; }
    DataType dataType() const { return dataType_; }
    double sampleRate() const { return sampleRate_; }
    double bitScaling() const { return bitScaling_; }
    std::size_t fileBytes() const { return file_.size(); }
    const std::string& path() const { return file_.path(); }

    /// Recording length in seconds.
    double duration() const
    {
        return sampleRate_ > 0 ? static_cast<double>(nSamples_) / sampleRate_ : 0.0;
    }

    /// Reads samples [sampleBegin, sampleEnd) for every channel.
    ///
    /// Output is channel-interleaved, matching the file: out[c + s*nChans].
    /// Values are raw ADC counts converted to float, without bitScaling applied.
    void readSamples(std::int64_t sampleBegin,
                     std::int64_t sampleEnd,
                     std::vector<float>& out) const;

    /// As readSamples, but multiplied by bitScaling so the result is in µV.
    void readSamplesMicrovolts(std::int64_t sampleBegin,
                               std::int64_t sampleEnd,
                               std::vector<float>& out) const;

    /// Reads one channel over [sampleBegin, sampleEnd) into `out`.
    void readChannel(int channel,
                     std::int64_t sampleBegin,
                     std::int64_t sampleEnd,
                     std::vector<float>& out) const;

    /// Reads a set of channels and transposes them in a single pass.
    ///
    /// `out` is filled site-major - out[site * nSamples + t] - where site `i`
    /// takes its data from channel `channels[i]`. This fuses the raw read, the
    /// integer-to-float conversion and the transpose, so the samples are only
    /// touched once.
    ///
    /// The loop is tiled over time: a tile's worth of source data stays in L2
    /// while every site reads its column out of it, which turns what would
    /// otherwise be one cache miss per sample into near-sequential access.
    void readSitesTransposed(std::int64_t sampleBegin,
                             std::int64_t sampleEnd,
                             const std::vector<int>& channels,
                             std::vector<float>& out) const;

private:
    /// Converts `count` samples starting at `byteOffset` into floats.
    void convert(std::size_t byteOffset, std::size_t count, std::ptrdiff_t stride,
                 float* out) const;

    MemoryMappedFile file_;
    int nChans_ = 0;
    std::int64_t nSamples_ = 0;
    DataType dataType_ = DataType::Int16;
    std::size_t headerOffset_ = 0;
    double sampleRate_ = 0.0;
    double bitScaling_ = 1.0;
};

}  // namespace jrc
