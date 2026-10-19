#include "detect/PeakDetector.h"

#include "core/Parallel.h"

#include <algorithm>
#include <cmath>
#include <numeric>

namespace jrc {
namespace detect {

void SpikeTable::clear()
{
    times.clear();
    amplitudes.clear();
    sites.clear();
}

void SpikeTable::reserve(std::size_t n)
{
    times.reserve(n);
    amplitudes.reserve(n);
    sites.reserve(n);
}

std::vector<float> computeThresholds(const std::vector<float>& samples,
                                     int nSamples,
                                     int nSites,
                                     double qqFactor,
                                     const std::vector<std::uint8_t>& ignoredMask,
                                     int maxSamples)
{
    std::vector<float> thresholds(nSites, 0.0f);
    if (nSamples <= 0 || nSites <= 0) return thresholds;

    // MATLAB's estimateRMS subsamples with randsample(), which makes its
    // thresholds differ from run to run. We take an evenly spaced subsample
    // instead: the same estimator for stationary noise, but reproducible.
    const int stride =
        (maxSamples > 0 && nSamples > maxSamples) ? (nSamples / maxSamples) : 1;

    const float* base = samples.data();

    parallelFor(static_cast<std::size_t>(nSites), [&](std::size_t site) {
        if (!ignoredMask.empty() && ignoredMask[site]) return;  // stays 0

        const float* row = base + site * static_cast<std::size_t>(nSamples);

        std::vector<float> scratch;
        scratch.reserve(static_cast<std::size_t>(nSamples / stride) + 1);
        for (int s = 0; s < nSamples; s += stride) scratch.push_back(std::fabs(row[s]));

        if (scratch.empty()) return;

        const std::size_t mid = scratch.size() / 2;
        std::nth_element(scratch.begin(), scratch.begin() + mid, scratch.end());

        float median = scratch[mid];
        if (scratch.size() % 2 == 0 && mid > 0) {
            const float lower =
                *std::max_element(scratch.begin(), scratch.begin() + mid);
            median = 0.5f * (lower + median);
        }

        // estimateRMS: median(|x|) / 0.6745, scaled by qqFactor. MATLAB casts
        // the result to int16, so round to match its thresholds.
        const double rms = static_cast<double>(median) / 0.6745;
        thresholds[site] = static_cast<float>(std::round(rms * qqFactor));
    });

    return thresholds;
}

namespace {

/// Threshold crossings on one site that are also local minima.
///
/// Port of the local findPeaks() in detectPeaksSite.m. Times are indices into
/// the chunk. Detection runs on the negative-going deflection; `sign` is -1 for
/// the normal pass and +1 for the bipolar pass, which flips the data.
void findPeaksOnSite(const std::vector<float>& samples,
                     int nSamples,
                     int nSites,
                     int site,
                     float threshold,
                     int minNeighborsDetect,
                     float sign,
                     std::vector<std::int64_t>& outTimes,
                     std::vector<float>& outAmps)
{
    (void)nSites;  // layout is site-major; each trace is contiguous
    if (threshold == 0.0f) return;  // bad or ignored site
    const float limit = -std::fabs(threshold);

    const float* row = samples.data() + static_cast<std::size_t>(site) * nSamples;

    auto at = [&](int s) { return sign * row[s]; };
    auto exceeds = [&](int s) { return at(s) < limit; };

    // Interior only: a peak needs both neighbours, and the first and last
    // samples cannot have them. MATLAB drops these too, though its version
    // also removes the final candidate unconditionally - a bug we do not
    // reproduce, since it discards one real spike per site per chunk.
    for (int s = 1; s < nSamples - 1; ++s) {
        if (!exceeds(s)) continue;

        const float centre = at(s);

        // Local minimum, ties included (matches MATLAB's <=).
        if (!(centre <= at(s + 1) && centre <= at(s - 1))) continue;

        if (minNeighborsDetect == 1) {
            if (!(exceeds(s - 1) || exceeds(s + 1))) continue;
        } else if (minNeighborsDetect == 2) {
            if (!(exceeds(s - 1) && exceeds(s + 1))) continue;
        }

        outTimes.push_back(s);
        // Report the amplitude as it appears in the data, not sign-flipped.
        outAmps.push_back(row[s]);
    }
}

}  // namespace

SpikeTable detectSpikes(const std::vector<float>& samples,
                        int nSamples,
                        int nSites,
                        const std::vector<float>& thresholds,
                        const Probe& probe,
                        const DetectParams& params)
{
    // -- per-site candidate peaks ------------------------------------------
    std::vector<std::vector<std::int64_t>> siteTimes(nSites);
    std::vector<std::vector<float>> siteAmps(nSites);

    parallelFor(static_cast<std::size_t>(nSites), [&](std::size_t siteIndex) {
        const int site = static_cast<int>(siteIndex);

        findPeaksOnSite(samples, nSamples, nSites, site, thresholds[site],
                        params.minNeighborsDetect, 1.0f, siteTimes[site],
                        siteAmps[site]);

        if (params.detectBipolar) {
            std::vector<std::int64_t> posTimes;
            std::vector<float> posAmps;
            findPeaksOnSite(samples, nSamples, nSites, site, thresholds[site],
                            params.minNeighborsDetect, -1.0f, posTimes, posAmps);

            // Merge the two passes back into ascending time order.
            std::vector<std::size_t> order(siteTimes[site].size() + posTimes.size());
            std::vector<std::int64_t> mergedTimes;
            std::vector<float> mergedAmps;
            mergedTimes.reserve(order.size());
            mergedAmps.reserve(order.size());

            std::size_t a = 0, b = 0;
            while (a < siteTimes[site].size() || b < posTimes.size()) {
                const bool takeA = b >= posTimes.size() ||
                                   (a < siteTimes[site].size() &&
                                    siteTimes[site][a] <= posTimes[b]);
                if (takeA) {
                    mergedTimes.push_back(siteTimes[site][a]);
                    mergedAmps.push_back(siteAmps[site][a]);
                    ++a;
                } else {
                    mergedTimes.push_back(posTimes[b]);
                    mergedAmps.push_back(posAmps[b]);
                    ++b;
                }
            }
            siteTimes[site] = std::move(mergedTimes);
            siteAmps[site] = std::move(mergedAmps);
        }

        // Reject implausibly large events, as detectPeaksSite.m does when
        // spikeThreshMax is set.
        if (params.spikeThreshMaxRaw > 0) {
            const float limit = static_cast<float>(params.spikeThreshMaxRaw);
            std::size_t kept = 0;
            for (std::size_t i = 0; i < siteTimes[site].size(); ++i) {
                if (std::fabs(siteAmps[site][i]) < limit) {
                    siteTimes[site][kept] = siteTimes[site][i];
                    siteAmps[site][kept] = siteAmps[site][i];
                    ++kept;
                }
            }
            siteTimes[site].resize(kept);
            siteAmps[site].resize(kept);
        }
    });

    // -- resolve duplicates across neighbouring sites ----------------------
    // Port of mergeSpikesSite() in mergePeaks.m. A spike on site i is dropped
    // when a spike on a nearby site, within the refractory window, has a larger
    // magnitude; equal magnitudes are broken by time, then by site index.
    const std::int64_t nLims = std::abs(params.refracIntSamples);

    std::vector<std::vector<std::uint8_t>> keep(nSites);
    for (int site = 0; site < nSites; ++site)
        keep[site].assign(siteTimes[site].size(), 1);

    // Each iteration writes only keep[i] and reads its neighbours' peak lists,
    // which are complete and immutable by now, so sites are independent.
    parallelFor(static_cast<std::size_t>(nSites), [&](std::size_t siteIndex) {
        const int i = static_cast<int>(siteIndex);

        const auto& iTimes = siteTimes[i];
        const auto& iAmps = siteAmps[i];
        if (iTimes.empty()) return;

        const std::vector<int>& neighbors =
            probe.hasNeighborCache() ? probe.cachedNeighbors(i)
                                     : probe.nearbySites(i, params.evtDetectRad);

        for (int j : neighbors) {
            const auto& jTimes = siteTimes[j];
            const auto& jAmps = siteAmps[j];
            if (jTimes.empty()) continue;

            // Both lists are time-sorted, so sweep them together.
            std::size_t lo = 0;
            for (std::size_t a = 0; a < iTimes.size(); ++a) {
                if (!keep[i][a]) continue;

                const std::int64_t ti = iTimes[a];
                const float ai = std::fabs(iAmps[a]);

                while (lo < jTimes.size() && jTimes[lo] < ti - nLims) ++lo;

                for (std::size_t b = lo; b < jTimes.size(); ++b) {
                    const std::int64_t delay = jTimes[b] - ti;
                    if (delay > nLims) break;

                    if (i == j && delay == 0) continue;  // self

                    const float aj = std::fabs(jAmps[b]);

                    if (aj > ai) {
                        keep[i][a] = 0;
                        break;
                    }
                    if (aj == ai) {
                        // Earlier spike wins; at equal times, lower site index.
                        if (delay < 0 || (delay == 0 && j < i)) {
                            keep[i][a] = 0;
                            break;
                        }
                    }
                }
            }
        }
    });

    // -- collect survivors, ordered by time --------------------------------
    std::size_t total = 0;
    for (int site = 0; site < nSites; ++site)
        total += static_cast<std::size_t>(
            std::count(keep[site].begin(), keep[site].end(), 1));

    SpikeTable table;
    table.reserve(total);

    for (int site = 0; site < nSites; ++site) {
        for (std::size_t k = 0; k < siteTimes[site].size(); ++k) {
            if (!keep[site][k]) continue;
            table.times.push_back(siteTimes[site][k]);
            table.amplitudes.push_back(siteAmps[site][k]);
            table.sites.push_back(site);
        }
    }

    std::vector<std::size_t> order(table.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        if (table.times[a] != table.times[b]) return table.times[a] < table.times[b];
        return table.sites[a] < table.sites[b];
    });

    SpikeTable sorted;
    sorted.reserve(table.size());
    for (std::size_t idx : order) {
        sorted.times.push_back(table.times[idx]);
        sorted.amplitudes.push_back(table.amplitudes[idx]);
        sorted.sites.push_back(table.sites[idx]);
    }

    if (!params.removeCoincident || params.coincidenceSamples <= 0) return sorted;

    // -- shank-wide coincidence removal ------------------------------------
    // Deliberately confined to pairs FURTHER APART than evtDetectRad. Inside
    // that radius the loop above has already applied JRCLUST's own rule, and
    // re-judging those pairs on a different criterion would silently diverge
    // from the MATLAB result. This pass only decides cases JRCLUST never
    // looks at: two detections too far apart to be neighbours, yet landing on
    // the same sample with the same amplitude, which no pair of neurons can
    // produce and no single neuron can produce either.
    const double radSq = params.evtDetectRad * params.evtDetectRad;
    const std::size_t n = sorted.size();
    std::vector<std::uint8_t> alive(n, 1);

    for (std::size_t a = 0; a < n; ++a) {
        if (!alive[a]) continue;
        const int siteA = sorted.sites[a];
        const float ampA = std::fabs(sorted.amplitudes[a]);

        for (std::size_t b = a + 1; b < n; ++b) {
            if (sorted.times[b] - sorted.times[a] > params.coincidenceSamples) break;
            if (!alive[b]) continue;

            const int siteB = sorted.sites[b];
            if (siteA == siteB) continue;
            if (probe.shank(siteA) != probe.shank(siteB)) continue;

            const double dx = probe.x(siteA) - probe.x(siteB);
            const double dy = probe.y(siteA) - probe.y(siteB);
            if (dx * dx + dy * dy <= radSq) continue;  // already resolved above

            const float ampB = std::fabs(sorted.amplitudes[b]);
            const float hi = std::max(ampA, ampB);
            const float lo = std::min(ampA, ampB);
            if (hi <= 0 || lo / hi < params.coincidenceAmpRatio) continue;

            // Keep the larger, as JRCLUST does within its radius. Ties go to
            // the earlier spike, then to the lower site index, matching the
            // tie-break used above so the two passes cannot disagree.
            if (ampB > ampA) {
                alive[a] = 0;
                break;
            }
            alive[b] = 0;
        }
    }

    std::size_t survivors = 0;
    for (std::size_t i = 0; i < n; ++i) survivors += alive[i];
    if (survivors == n) return sorted;

    SpikeTable deduped;
    deduped.reserve(survivors);
    for (std::size_t i = 0; i < n; ++i) {
        if (!alive[i]) continue;
        deduped.times.push_back(sorted.times[i]);
        deduped.amplitudes.push_back(sorted.amplitudes[i]);
        deduped.sites.push_back(sorted.sites[i]);
    }

    return deduped;
}

}  // namespace detect
}  // namespace jrc
