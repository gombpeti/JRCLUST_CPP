#include "core/Recording.h"

#include "core/Parallel.h"

#include <algorithm>
#include <cstring>
#include <stdexcept>
#include <string>

namespace jrc {

namespace {

/// Reads one sample of type T from a possibly unaligned address.
template <typename T>
inline T loadUnaligned(const std::uint8_t* p)
{
    T value;
    std::memcpy(&value, p, sizeof(T));
    return value;
}

template <typename T>
void convertTyped(const std::uint8_t* base, std::size_t count, std::ptrdiff_t stride,
                  float* out)
{
    for (std::size_t i = 0; i < count; ++i) {
        out[i] = static_cast<float>(
            loadUnaligned<T>(base + static_cast<std::ptrdiff_t>(i) * stride));
    }
}

}  // namespace

void Recording::open(const std::string& path, const Config& config)
{
    close();

    nChans_ = config.nChans();
    dataType_ = config.dataType();
    headerOffset_ = config.headerOffset();
    sampleRate_ = config.sampleRate();
    bitScaling_ = config.bitScaling();

    if (nChans_ <= 0)
        throw std::runtime_error("nChans must be positive (check your .prm)");

    file_.open(path);

    const std::size_t bytes = file_.size();
    if (bytes <= headerOffset_)
        throw std::runtime_error("file is smaller than headerOffset: " + path);

    const std::size_t payload = bytes - headerOffset_;
    const std::size_t frameBytes =
        static_cast<std::size_t>(nChans_) * dataTypeBytes(dataType_);

    if (payload % frameBytes != 0) {
        // Same check as Recording.m: a partial frame means the geometry is wrong.
        throw std::runtime_error(
            "recording size is not a whole number of samples - check nChans (" +
            std::to_string(nChans_) + "), dataType (" +
            std::string(dataTypeName(dataType_)) + ") or headerOffset (" +
            std::to_string(headerOffset_) + ")");
    }

    nSamples_ = static_cast<std::int64_t>(payload / frameBytes);
}

void Recording::close()
{
    file_.close();
    nSamples_ = 0;
}

void Recording::convert(std::size_t byteOffset, std::size_t count,
                        std::ptrdiff_t stride, float* out) const
{
    const std::uint8_t* base = file_.data() + byteOffset;

    switch (dataType_) {
        case DataType::Int16: convertTyped<std::int16_t>(base, count, stride, out); break;
        case DataType::Uint16: convertTyped<std::uint16_t>(base, count, stride, out); break;
        case DataType::Int32: convertTyped<std::int32_t>(base, count, stride, out); break;
        case DataType::Uint32: convertTyped<std::uint32_t>(base, count, stride, out); break;
        case DataType::Single: convertTyped<float>(base, count, stride, out); break;
        case DataType::Double: convertTyped<double>(base, count, stride, out); break;
    }
}

void Recording::readSamples(std::int64_t sampleBegin, std::int64_t sampleEnd,
                            std::vector<float>& out) const
{
    if (!isOpen()) throw std::runtime_error("recording is not open");
    if (sampleBegin < 0 || sampleEnd > nSamples_ || sampleBegin >= sampleEnd)
        throw std::runtime_error("sample range out of bounds");

    const std::size_t sampleBytes = dataTypeBytes(dataType_);
    const std::int64_t nRequested = sampleEnd - sampleBegin;
    const std::size_t count =
        static_cast<std::size_t>(nRequested) * static_cast<std::size_t>(nChans_);

    out.resize(count);

    // Channel-interleaved and contiguous, so this is one straight run.
    const std::size_t offset =
        headerOffset_ + static_cast<std::size_t>(sampleBegin) *
                            static_cast<std::size_t>(nChans_) * sampleBytes;

    convert(offset, count, static_cast<std::ptrdiff_t>(sampleBytes), out.data());
}

void Recording::readSamplesMicrovolts(std::int64_t sampleBegin, std::int64_t sampleEnd,
                                      std::vector<float>& out) const
{
    readSamples(sampleBegin, sampleEnd, out);

    const float scale = static_cast<float>(bitScaling_);
    for (auto& value : out) value *= scale;
}

namespace {

/// Converts and transposes one time tile for every requested channel.
template <typename T>
void transposeTile(const std::uint8_t* base, int nChans, const std::vector<int>& channels,
                   float* out, std::int64_t nSamplesOut, std::int64_t tileBegin,
                   std::int64_t tileEnd)
{
    const auto* src = reinterpret_cast<const T*>(base);
    const std::size_t nSites = channels.size();

    for (std::size_t site = 0; site < nSites; ++site) {
        const int channel = channels[site];
        float* row = out + site * static_cast<std::size_t>(nSamplesOut);

        for (std::int64_t s = tileBegin; s < tileEnd; ++s) {
            T value;
            std::memcpy(&value, src + s * nChans + channel, sizeof(T));
            row[s] = static_cast<float>(value);
        }
    }
}

}  // namespace

void Recording::readSitesTransposed(std::int64_t sampleBegin, std::int64_t sampleEnd,
                                    const std::vector<int>& channels,
                                    std::vector<float>& out) const
{
    if (!isOpen()) throw std::runtime_error("recording is not open");
    if (sampleBegin < 0 || sampleEnd > nSamples_ || sampleBegin >= sampleEnd)
        throw std::runtime_error("sample range out of bounds");

    for (int channel : channels) {
        if (channel < 0 || channel >= nChans_)
            throw std::runtime_error("channel out of range in site map");
    }

    const std::int64_t nSamplesOut = sampleEnd - sampleBegin;
    const std::size_t nSites = channels.size();
    out.resize(static_cast<std::size_t>(nSamplesOut) * nSites);

    const std::size_t sampleBytes = dataTypeBytes(dataType_);
    const std::uint8_t* base =
        file_.data() + headerOffset_ +
        static_cast<std::size_t>(sampleBegin) * nChans_ * sampleBytes;

    // 4096 samples of 512 int16 channels is 4 MB - comfortably inside L2 on a
    // modern core, so each site's strided column read hits cache.
    constexpr std::int64_t kTile = 4096;
    const std::int64_t nTiles = (nSamplesOut + kTile - 1) / kTile;

    float* dst = out.data();
    const int nChans = nChans_;
    const DataType type = dataType_;

    parallelFor(static_cast<std::size_t>(nTiles), [&](std::size_t tileIndex) {
        const std::int64_t tileBegin = static_cast<std::int64_t>(tileIndex) * kTile;
        const std::int64_t tileEnd = std::min(tileBegin + kTile, nSamplesOut);

        switch (type) {
            case DataType::Int16:
                transposeTile<std::int16_t>(base, nChans, channels, dst, nSamplesOut,
                                            tileBegin, tileEnd);
                break;
            case DataType::Uint16:
                transposeTile<std::uint16_t>(base, nChans, channels, dst, nSamplesOut,
                                             tileBegin, tileEnd);
                break;
            case DataType::Int32:
                transposeTile<std::int32_t>(base, nChans, channels, dst, nSamplesOut,
                                            tileBegin, tileEnd);
                break;
            case DataType::Uint32:
                transposeTile<std::uint32_t>(base, nChans, channels, dst, nSamplesOut,
                                             tileBegin, tileEnd);
                break;
            case DataType::Single:
                transposeTile<float>(base, nChans, channels, dst, nSamplesOut,
                                     tileBegin, tileEnd);
                break;
            case DataType::Double:
                transposeTile<double>(base, nChans, channels, dst, nSamplesOut,
                                      tileBegin, tileEnd);
                break;
        }
    });
}

void Recording::readChannel(int channel, std::int64_t sampleBegin,
                            std::int64_t sampleEnd, std::vector<float>& out) const
{
    if (!isOpen()) throw std::runtime_error("recording is not open");
    if (channel < 0 || channel >= nChans_)
        throw std::runtime_error("channel out of range");
    if (sampleBegin < 0 || sampleEnd > nSamples_ || sampleBegin >= sampleEnd)
        throw std::runtime_error("sample range out of bounds");

    const std::size_t sampleBytes = dataTypeBytes(dataType_);
    const std::size_t count = static_cast<std::size_t>(sampleEnd - sampleBegin);

    out.resize(count);

    const std::size_t offset =
        headerOffset_ +
        (static_cast<std::size_t>(sampleBegin) * static_cast<std::size_t>(nChans_) +
         static_cast<std::size_t>(channel)) *
            sampleBytes;

    // Stride one whole frame to walk a single channel down the file.
    const std::ptrdiff_t stride =
        static_cast<std::ptrdiff_t>(sampleBytes) * static_cast<std::ptrdiff_t>(nChans_);

    convert(offset, count, stride, out.data());
}

}  // namespace jrc
