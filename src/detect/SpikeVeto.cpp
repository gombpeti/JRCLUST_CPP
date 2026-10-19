#include "detect/SpikeVeto.h"

#include "core/Parallel.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace jrc {
namespace detect {

namespace {

/// Trough depth of one trace, as a positive magnitude, plus its index.
inline float troughOf(const float* w, int n, int& index)
{
    index = 0;
    for (int s = 1; s < n; ++s)
        if (w[s] < w[index]) index = s;
    return -w[index];
}

}  // namespace

VetoCounts vetoSpikes(const features::SpikeWindows& rawWindows,
                      const features::SpikeWindows& filtWindows,
                      const SpikeVetoParams& params, double sampleRate,
                      std::vector<std::uint8_t>& keep)
{
    const features::SpikeWindows& windows = rawWindows;

    VetoCounts counts;
    keep.assign(static_cast<std::size_t>(std::max(0, windows.nSpikes)), 1);

    if (!params.enabled || windows.nSpikes == 0 || windows.nSamples < 3)
        return counts;

    const int nSamples = windows.nSamples;
    const int nSites = windows.nSites;
    const double msPerSample = 1000.0 / sampleRate;

    // Far half of the neighbourhood. siteNeighbors is ordered by distance, so
    // the outer slots are the distant sites.
    const int farBegin = (nSites + 1) / 2;
    const bool canTestDecay = nSites >= 4 && farBegin < nSites &&
                              filtWindows.nSpikes == windows.nSpikes &&
                              filtWindows.nSites == nSites &&
                              filtWindows.nSamples >= 3;

    // Per-thread tallies, combined afterwards, so the counting does not need
    // synchronisation on the hot path.
    const unsigned workers =
        static_cast<unsigned>(std::min<std::size_t>(threadCount(), windows.nSpikes));
    std::vector<VetoCounts> partial(std::max(1u, workers));

    parallelFor(std::max(1u, workers), [&](std::size_t worker) {
        VetoCounts& tally = partial[worker];

        for (int spike = static_cast<int>(worker); spike < windows.nSpikes;
             spike += static_cast<int>(std::max(1u, workers))) {
            ++tally.tested;

            const float* w = windows.at(0, spike);  // the spike's own site

            int troughIndex = 0;
            const float depth = troughOf(w, nSamples, troughIndex);

            bool reject = false;

            if (depth > 0) {
                // -- width at half depth ------------------------------------
                const float half = -0.5f * depth;
                int left = troughIndex;
                while (left > 0 && w[left - 1] <= half) --left;
                int right = troughIndex;
                while (right + 1 < nSamples && w[right + 1] <= half) ++right;

                const double halfWidthMs = (right - left + 1) * msPerSample;
                if (halfWidthMs < params.minHalfWidthMs) {
                    ++tally.tooNarrow;
                    reject = true;
                }

                // -- trough versus its immediate neighbours ------------------
                const double before =
                    std::fabs(troughIndex > 0 ? w[troughIndex - 1] : w[troughIndex]);
                const double after = std::fabs(
                    troughIndex + 1 < nSamples ? w[troughIndex + 1] : w[troughIndex]);
                const double neighbourhood = 0.5 * (before + after);

                const double spikiness =
                    neighbourhood > 1e-9 ? depth / neighbourhood
                                         : std::numeric_limits<double>::infinity();
                if (spikiness > params.maxSpikiness) {
                    ++tally.tooSpiky;
                    reject = true;
                }

                // -- spatial decay -------------------------------------------
                // A point source falls off with distance. Equal amplitude on
                // near and far sites means the event has no spatial origin.
                //
                // Measured on the FILTERED windows: the raw trace carries
                // low-frequency activity common to the whole probe, which
                // would make every site look equally large.
                if (canTestDecay) {
                    const int fSamples = filtWindows.nSamples;
                    int unused = 0;
                    const double near =
                        troughOf(filtWindows.at(0, spike), fSamples, unused);

                    if (near >= params.minAmplitudeForDecay) {
                        double farSum = 0.0;
                        int farCount = 0;
                        for (int k = farBegin; k < nSites; ++k) {
                            const double d = troughOf(filtWindows.at(k, spike),
                                                      fSamples, unused);
                            // Ignored sites are zeroed; they carry no
                            // information either way, so leave them out.
                            if (d > 0) {
                                farSum += d;
                                ++farCount;
                            }
                        }

                        if (farCount > 0) {
                            const double decayRatio = (farSum / farCount) / near;
                            if (decayRatio > params.maxDecayRatio) {
                                ++tally.tooFlat;
                                reject = true;
                            }
                        }
                    }
                }
            }

            if (reject) {
                keep[spike] = 0;
                ++tally.rejected;
            }
        }
    });

    for (const auto& p : partial) counts.add(p);
    return counts;
}

std::int64_t findRecordingGaps(const std::vector<float>& samples, int nSamples,
                               int nSites, double zeroFrac, int padSamples,
                               std::vector<std::uint8_t>& mask)
{
    mask.assign(static_cast<std::size_t>(std::max(0, nSamples)), 0);
    if (nSamples <= 0 || nSites <= 0) return 0;

    // A gap sample is one where essentially every site reads exactly zero.
    // Counting per sample is a strided read, so the sites loop is outer and
    // the counts accumulate across a contiguous pass per site.
    std::vector<std::int32_t> zeroCount(nSamples, 0);

    const float* base = samples.data();
    const unsigned workers =
        static_cast<unsigned>(std::min<std::size_t>(threadCount(), nSites));
    std::vector<std::vector<std::int32_t>> partial(
        std::max(1u, workers), std::vector<std::int32_t>(nSamples, 0));

    parallelFor(std::max(1u, workers), [&](std::size_t worker) {
        std::vector<std::int32_t>& acc = partial[worker];
        for (int site = static_cast<int>(worker); site < nSites;
             site += static_cast<int>(std::max(1u, workers))) {
            const float* row = base + static_cast<std::size_t>(site) * nSamples;
            for (int s = 0; s < nSamples; ++s)
                if (row[s] == 0.0f) ++acc[s];
        }
    });

    for (const auto& acc : partial)
        for (int s = 0; s < nSamples; ++s) zeroCount[s] += acc[s];

    const int needed = static_cast<int>(std::ceil(zeroFrac * nSites));

    std::vector<std::uint8_t> raw(nSamples, 0);
    for (int s = 0; s < nSamples; ++s)
        if (zeroCount[s] >= needed) raw[s] = 1;

    // Widen by the margin so that no spike whose waveform merely touches the
    // step survives. Done as a single pass over the run boundaries rather
    // than by dilating every sample.
    std::int64_t marked = 0;
    int s = 0;
    while (s < nSamples) {
        if (!raw[s]) { ++s; continue; }

        int end = s;
        while (end + 1 < nSamples && raw[end + 1]) ++end;

        const int from = std::max(0, s - padSamples);
        const int to = std::min(nSamples - 1, end + padSamples);
        for (int k = from; k <= to; ++k) {
            if (!mask[k]) ++marked;
            mask[k] = 1;
        }

        s = end + 1;
    }

    return marked;
}

void compactWindows(features::SpikeWindows& windows,
                    const std::vector<std::uint8_t>& keep)
{
    const int nSpikes = windows.nSpikes;
    if (nSpikes == 0) return;

    int kept = 0;
    for (int i = 0; i < nSpikes; ++i)
        if (keep[i]) ++kept;

    if (kept == nSpikes) return;

    const int nSamples = windows.nSamples;
    const int nSites = windows.nSites;

    std::vector<float> packed(static_cast<std::size_t>(nSites) * kept * nSamples);

    // Layout is [site][spike][sample], so each site's surviving spikes are
    // gathered into a contiguous run.
    for (int site = 0; site < nSites; ++site) {
        float* dst = packed.data() + static_cast<std::size_t>(site) * kept * nSamples;
        int out = 0;
        for (int spike = 0; spike < nSpikes; ++spike) {
            if (!keep[spike]) continue;
            std::copy_n(windows.at(site, spike), nSamples,
                        dst + static_cast<std::size_t>(out) * nSamples);
            ++out;
        }
    }

    windows.data = std::move(packed);
    windows.nSpikes = kept;
}

}  // namespace detect
}  // namespace jrc
