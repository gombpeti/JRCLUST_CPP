#include "sort/PostProcess.h"

#include "core/Parallel.h"
#include "core/Stats.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace jrc {
namespace sort {

PostProcessSettings PostProcessSettings::fromConfig(const Config& config)
{
    PostProcessSettings settings;
    settings.maxUnitSim = config.getNumber("maxUnitSim", 0.98);
    settings.nPassesMerge = config.getInt("nPassesMerge", 10);
    settings.autoMergeBy = config.getString("autoMergeBy", "pearson");
    settings.refracInt = config.getNumber("refracInt", 0.25);
    settings.evtMergeRad = config.getNumber("evtMergeRad", 35.0);
    settings.meanInterpFactor = config.getInt("meanInterpFactor", 1);
    settings.nClusterIntervals = config.getInt("nClusterIntervals", 4);
    settings.outlierThresh = config.getNumber("outlierThresh", 7.5);
    settings.qqFactor = config.getNumber("qqFactor", 5.0);
    settings.bitScaling = config.getNumber("bitScaling", 0.195);
    settings.sampleRate = config.sampleRate();
    settings.driftMerge = config.getInt("driftMerge", 1) != 0;

    // bit2uV.m: the differentiator filters scale amplitudes, so uV conversion
    // divides that factor back out.
    const std::string filterType = config.getString("filterType", "ndiff");
    const int nDiffOrder = config.getInt("nDiffOrder", 2);
    if (filterType == "sgdiff") {
        double sum = 0;
        for (int i = 1; i <= nDiffOrder; ++i) sum += double(i) * i;
        settings.filterNorm = sum * 2.0;
    } else if (filterType == "ndiff") {
        settings.filterNorm = 2.0;
    } else {
        settings.filterNorm = 1.0;
    }

    return settings;
}

namespace {

/// MATLAB's round(): halves go away from zero, unlike std::round's ties-away
/// which happens to match. Kept explicit so the intent is clear.
inline int matlabRound(double v)
{
    return static_cast<int>(v < 0 ? -std::floor(-v + 0.5) : std::floor(v + 0.5));
}

/// Spike indices of each cluster, 0-based, in ascending order.
std::vector<std::vector<std::int32_t>> spikesByCluster(const ClusterResult& clusters,
                                                       std::size_t nSpikes)
{
    std::vector<std::vector<std::int32_t>> out(clusters.nClusters);
    for (std::size_t i = 0; i < nSpikes; ++i) {
        const std::int32_t c = clusters.spikeClusters[i];
        if (c > 0 && c <= clusters.nClusters)
            out[c - 1].push_back(static_cast<std::int32_t>(i));
    }
    return out;
}

/// The cluster's spikes nearest the middle of the recording.
/// Port of getMidmostSpikes.m.
std::vector<std::int32_t> midmostSpikes(const std::vector<std::int32_t>& spikes,
                                        std::size_t nSpikesTotal, int nClusterIntervals)
{
    const int keep = std::max(1, matlabRound(static_cast<double>(spikes.size()) /
                                             std::max(1, nClusterIntervals)));
    if (static_cast<int>(spikes.size()) <= keep) return spikes;

    const double middle = matlabRound(static_cast<double>(nSpikesTotal) / 2.0);

    std::vector<std::int32_t> ordered = spikes;
    std::nth_element(ordered.begin(), ordered.begin() + keep, ordered.end(),
                     [&](std::int32_t a, std::int32_t b) {
                         return std::fabs(a - middle) < std::fabs(b - middle);
                     });
    ordered.resize(keep);
    return ordered;
}

/// Spikes whose y-position falls in [lo, hi), capped at 1000 and spread evenly.
/// Port of getSpikesInBounds.m, which falls back to everything if the band is
/// empty. MATLAB subsamples at random here; even spacing keeps runs repeatable.
std::vector<std::int32_t> spikesInBand(const std::vector<std::int32_t>& spikes,
                                       const std::vector<float>& yPos, float lo,
                                       float hi)
{
    std::vector<std::int32_t> inBand;
    for (std::size_t i = 0; i < spikes.size(); ++i)
        if (yPos[i] >= lo && yPos[i] < hi) inBand.push_back(spikes[i]);

    const std::vector<std::int32_t>& source = inBand.empty() ? spikes : inBand;

    constexpr std::size_t kMax = 1000;
    if (source.size() <= kMax) return source;

    std::vector<std::int32_t> out;
    out.reserve(kMax);
    const double step = static_cast<double>(source.size()) / kMax;
    for (std::size_t i = 0; i < kMax; ++i)
        out.push_back(source[static_cast<std::size_t>(i * step)]);
    return out;
}

/// Averages spikes onto `iSite`'s neighbourhood, matching each spike's own
/// neighbourhood by absolute site id and ignoring sites the two do not share.
///
/// Port of nanmeanInt16.m. `fetch(spike, slot)` returns that spike's waveform
/// on slot `slot` of its OWN neighbourhood. Output is [nLocalSites][nSamples].
template <typename Fetch>
void averageAligned(const Fetch& fetch, const std::vector<std::int32_t>& spikes,
                    const std::vector<std::int32_t>& spikeSites, int iSite,
                    const Probe& probe, int nSamples, int nLocalSites,
                    std::vector<double>& out)
{
    out.assign(static_cast<std::size_t>(nLocalSites) * nSamples, 0.0);
    std::vector<int> counts(nLocalSites, 0);

    const int* iNeighbors = probe.siteNeighbors(iSite);
    const int nNeighbors = probe.nSiteNeighbors();

    // Group by the spike's own site so the slot mapping is built once per site
    // rather than once per spike.
    std::vector<std::size_t> order(spikes.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        return spikeSites[spikes[a]] < spikeSites[spikes[b]];
    });

    std::size_t i = 0;
    std::vector<std::pair<int, int>> slotMap;  // (slot in iNeighbors, slot in jNeighbors)

    while (i < order.size()) {
        const int jSite = spikeSites[spikes[order[i]]];
        std::size_t j = i;
        while (j < order.size() && spikeSites[spikes[order[j]]] == jSite) ++j;

        const int* jNeighbors = probe.siteNeighbors(jSite);
        slotMap.clear();
        for (int a = 0; a < nLocalSites && a < nNeighbors; ++a)
            for (int b = 0; b < nNeighbors; ++b)
                if (iNeighbors[a] == jNeighbors[b]) {
                    slotMap.emplace_back(a, b);
                    break;
                }

        if (!slotMap.empty()) {
            for (std::size_t k = i; k < j; ++k) {
                const std::int32_t spike = spikes[order[k]];
                for (const auto& pair : slotMap) {
                    const auto* w = fetch(spike, pair.second);
                    double* acc = out.data() +
                                  static_cast<std::size_t>(pair.first) * nSamples;
                    for (int s = 0; s < nSamples; ++s) acc[s] += w[s];
                    ++counts[pair.first];
                }
            }
            // counts is per (site, spike); each spike contributed nSamples
            // values to one slot, so divide by the spike count for that slot.
        }

        i = j;
    }

    for (int a = 0; a < nLocalSites; ++a) {
        double* acc = out.data() + static_cast<std::size_t>(a) * nSamples;
        if (counts[a] == 0) continue;
        const double inv = 1.0 / counts[a];
        for (int s = 0; s < nSamples; ++s) acc[s] *= inv;
    }
}

/// Shift index pairs for the correlation search. Port of shiftRange.m.
struct ShiftPair {
    std::vector<int> back;     // 0-based sample indices
    std::vector<int> forward;
};

std::vector<ShiftPair> buildShiftRange(int nSamples, int nShifts)
{
    std::vector<ShiftPair> pairs;
    pairs.reserve(2 * nShifts + 1);

    for (int shift = -nShifts; shift <= nShifts; ++shift) {
        const int back = -matlabRound(shift / 2.0);
        const int fwd = shift + back;

        ShiftPair pair;
        for (int s = 0; s < nSamples; ++s) {
            const int b = s + back;
            const int f = s + fwd;
            if (b >= 0 && b < nSamples && f >= 0 && f < nSamples) {
                pair.back.push_back(b);
                pair.forward.push_back(f);
            }
        }
        if (!pair.back.empty()) pairs.push_back(std::move(pair));
    }

    return pairs;
}

/// Pearson correlation of two equal-length sequences.
double pearson(const std::vector<float>& x, const std::vector<float>& y)
{
    const std::size_t n = x.size();
    if (n < 2 || y.size() != n) return 0.0;

    double sumX = 0, sumY = 0;
    for (std::size_t i = 0; i < n; ++i) {
        sumX += x[i];
        sumY += y[i];
    }
    const double meanX = sumX / n;
    const double meanY = sumY / n;

    double sxx = 0, syy = 0, sxy = 0;
    for (std::size_t i = 0; i < n; ++i) {
        const double dx = x[i] - meanX;
        const double dy = y[i] - meanY;
        sxx += dx * dx;
        syy += dy * dy;
        sxy += dx * dy;
    }

    const double denom = std::sqrt(sxx * syy);
    return denom > 0 ? sxy / denom : 0.0;
}

/// Regularized lower incomplete gamma P(a, x), for chi2cdf.
double gammaP(double a, double x)
{
    if (x <= 0 || a <= 0) return 0.0;

    const double logGammaA = std::lgamma(a);

    if (x < a + 1.0) {
        // Series expansion.
        double term = 1.0 / a;
        double sum = term;
        for (int n = 1; n < 1000; ++n) {
            term *= x / (a + n);
            sum += term;
            if (std::fabs(term) < std::fabs(sum) * 1e-15) break;
        }
        return sum * std::exp(-x + a * std::log(x) - logGammaA);
    }

    // Continued fraction for Q(a,x), then P = 1 - Q.
    constexpr double tiny = 1e-300;
    double b = x + 1.0 - a;
    double c = 1.0 / tiny;
    double d = 1.0 / b;
    double h = d;

    for (int i = 1; i < 1000; ++i) {
        const double an = -static_cast<double>(i) * (i - a);
        b += 2.0;
        d = an * d + b;
        if (std::fabs(d) < tiny) d = tiny;
        c = b + an / c;
        if (std::fabs(c) < tiny) c = tiny;
        d = 1.0 / d;
        const double delta = d * c;
        h *= delta;
        if (std::fabs(delta - 1.0) < 1e-15) break;
    }

    const double q = std::exp(-x + a * std::log(x) - logGammaA) * h;
    return 1.0 - q;
}

/// chi-squared CDF with `dof` degrees of freedom.
inline double chi2cdf(double x, double dof) { return gammaP(dof / 2.0, x / 2.0); }

/// Cholesky decomposition of a symmetric positive-definite matrix (row-major).
/// Returns false if the matrix is not positive definite.
bool cholesky(std::vector<double>& a, int n)
{
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= i; ++j) {
            double sum = a[static_cast<std::size_t>(i) * n + j];
            for (int k = 0; k < j; ++k)
                sum -= a[static_cast<std::size_t>(i) * n + k] *
                       a[static_cast<std::size_t>(j) * n + k];

            if (i == j) {
                if (sum <= 1e-12) return false;  // singular or indefinite
                a[static_cast<std::size_t>(i) * n + i] = std::sqrt(sum);
            } else {
                a[static_cast<std::size_t>(i) * n + j] =
                    sum / a[static_cast<std::size_t>(j) * n + j];
            }
        }
        for (int j = i + 1; j < n; ++j) a[static_cast<std::size_t>(i) * n + j] = 0.0;
    }
    return true;
}

}  // namespace

// ------------------------------------------------------------ mean waveforms --

MeanWaveforms computeMeanWaveforms(const detect::DetectResult& detected,
                                   const ClusterResult& clusters, const Probe& probe,
                                   const PostProcessSettings& settings)
{
    MeanWaveforms out;
    out.nSamples = detected.filtSamples;
    out.nSamplesRaw = detected.rawSamples;
    out.nSites = probe.nSites();
    out.nLocalSites = detected.nWaveformSites;
    out.nClusters = clusters.nClusters;

    if (out.nClusters == 0 || out.nSamples == 0) return out;

    const std::size_t nSpikes = detected.size();
    const auto byCluster = spikesByCluster(clusters, nSpikes);

    const std::size_t globalSize = static_cast<std::size_t>(out.nClusters) *
                                   out.nSites * out.nSamples;
    const std::size_t localSize = static_cast<std::size_t>(out.nClusters) *
                                  out.nLocalSites * out.nSamples;

    out.global.assign(globalSize, 0.0f);
    out.local.assign(localSize, 0.0f);

    const bool haveRaw = !detected.rawWaveforms.empty() && out.nSamplesRaw > 0;
    if (haveRaw) {
        out.globalRaw.assign(static_cast<std::size_t>(out.nClusters) * out.nSites *
                                 out.nSamplesRaw,
                             0.0f);
        out.localRaw.assign(static_cast<std::size_t>(out.nClusters) * out.nLocalSites *
                                out.nSamplesRaw,
                            0.0f);
    }

    // bit2uV: the filtered traces carry the differentiator's gain.
    const float filtScale =
        static_cast<float>(settings.bitScaling / settings.filterNorm);
    const float rawScale = static_cast<float>(settings.bitScaling);

    const bool useDrift =
        settings.driftMerge && detected.spikePositions.size() >= nSpikes * 2;

    parallelFor(static_cast<std::size_t>(out.nClusters), [&](std::size_t clusterIndex) {
        const int c = static_cast<int>(clusterIndex);
        const auto& spikes = byCluster[c];
        if (spikes.empty()) return;

        const int site = clusters.clusterSites[c];
        const int* neighbors = probe.siteNeighbors(site);

        // Which spikes to average, and how.
        //
        // Drift path (driftMerge): take the middle third of the cluster by
        // spike y-position and realign each spike onto this cluster's site
        // neighbourhood, so a unit that drifted across sites still averages
        // coherently. Otherwise: average the spikes nearest the middle of the
        // recording on their own neighbourhood slots.
        std::vector<std::int32_t> chosen;
        std::vector<float> yPos;
        float yLow = 0, yMid1 = 0, yMid2 = 0, yHigh = 0;

        if (useDrift) {
            yPos.resize(spikes.size());
            for (std::size_t i = 0; i < spikes.size(); ++i)
                yPos[i] = detected.spikePositions[static_cast<std::size_t>(spikes[i]) * 2 + 1];

            std::vector<float> sorted = yPos;
            yLow = matlabQuantile(sorted, 0.0);
            yMid1 = matlabQuantile(sorted, 1.0 / 3.0);
            yMid2 = matlabQuantile(sorted, 2.0 / 3.0);
            yHigh = matlabQuantile(sorted, 1.0);

            chosen = spikesInBand(spikes, yPos, yMid1, yMid2);
        } else {
            chosen = midmostSpikes(spikes, nSpikes, settings.nClusterIntervals);
        }
        if (chosen.empty()) return;

        std::vector<double> acc;

        const auto writeOut = [&](const std::vector<double>& values, int nSamples,
                                  float scale, float* localBase, float* globalBase,
                                  int globalStride) {
            for (int k = 0; k < out.nLocalSites; ++k) {
                const double* a = values.data() + static_cast<std::size_t>(k) * nSamples;

                // meanSubtract along the sample axis, then convert to uV.
                double mean = 0.0;
                for (int s = 0; s < nSamples; ++s) mean += a[s];
                mean /= nSamples;

                float* localRow =
                    localBase + static_cast<std::size_t>(k) * nSamples;
                float* globalRow =
                    globalBase + static_cast<std::size_t>(neighbors[k]) * globalStride;

                for (int s = 0; s < nSamples; ++s) {
                    const float v = static_cast<float>((a[s] - mean) * scale);
                    localRow[s] = v;
                    globalRow[s] = v;
                }
            }
        };

        // -- filtered ------------------------------------------------------
        const auto fetchFilt = [&](std::int32_t spike, int slot) {
            return detected.filtWaveform(spike, slot);
        };

        if (useDrift) {
            averageAligned(fetchFilt, chosen, detected.spikeSites, site, probe,
                           out.nSamples, out.nLocalSites, acc);
        } else {
            acc.assign(static_cast<std::size_t>(out.nSamples) * out.nLocalSites, 0.0);
            for (std::int32_t spike : chosen)
                for (int k = 0; k < out.nLocalSites; ++k) {
                    const std::int16_t* w = detected.filtWaveform(spike, k);
                    double* a = acc.data() + static_cast<std::size_t>(k) * out.nSamples;
                    for (int s = 0; s < out.nSamples; ++s) a[s] += w[s];
                }
            const double inv = 1.0 / static_cast<double>(chosen.size());
            for (double& v : acc) v *= inv;
        }

        writeOut(acc, out.nSamples, filtScale,
                 out.local.data() + (static_cast<std::size_t>(c) * out.nLocalSites) *
                                        out.nSamples,
                 out.global.data() + (static_cast<std::size_t>(c) * out.nSites) *
                                         out.nSamples,
                 out.nSamples);

        if (!haveRaw) return;

        // -- raw -----------------------------------------------------------
        const auto fetchRaw = [&](std::int32_t spike, int slot) {
            return detected.rawWaveform(spike, slot);
        };

        if (useDrift) {
            averageAligned(fetchRaw, chosen, detected.spikeSites, site, probe,
                           out.nSamplesRaw, out.nLocalSites, acc);
        } else {
            acc.assign(static_cast<std::size_t>(out.nSamplesRaw) * out.nLocalSites, 0.0);
            for (std::int32_t spike : chosen)
                for (int k = 0; k < out.nLocalSites; ++k) {
                    const std::int16_t* w = detected.rawWaveform(spike, k);
                    double* a =
                        acc.data() + static_cast<std::size_t>(k) * out.nSamplesRaw;
                    for (int s = 0; s < out.nSamplesRaw; ++s) a[s] += w[s];
                }
            const double inv = 1.0 / static_cast<double>(chosen.size());
            for (double& v : acc) v *= inv;
        }

        writeOut(acc, out.nSamplesRaw, rawScale,
                 out.localRaw.data() + (static_cast<std::size_t>(c) * out.nLocalSites) *
                                           out.nSamplesRaw,
                 out.globalRaw.data() + (static_cast<std::size_t>(c) * out.nSites) *
                                            out.nSamplesRaw,
                 out.nSamplesRaw);

        (void)yLow;
        (void)yHigh;
    });

    return out;
}

// -------------------------------------------------------- waveform similarity --

namespace {

/// The three strongest sites per cluster. Port of getSecondaryPeaks.m: rank
/// sites by how far the mean waveform dips below its own first sample.
void secondaryPeaks(const MeanWaveforms& waveforms, std::vector<int>& sites1,
                    std::vector<int>& sites2, std::vector<int>& sites3)
{
    const int nClusters = waveforms.nClusters;
    const int nSites = waveforms.nSites;

    sites1.assign(nClusters, 0);
    sites2.assign(nClusters, 0);
    sites3.assign(nClusters, 0);

    std::vector<float> depth(nSites);

    for (int c = 0; c < nClusters; ++c) {
        for (int site = 0; site < nSites; ++site) {
            const float* w = waveforms.globalAt(c, site);
            float lo = w[0];
            for (int s = 1; s < waveforms.nSamples; ++s) lo = std::min(lo, w[s]);
            depth[site] = lo - w[0];
        }

        int best = 0;
        for (int site = 1; site < nSites; ++site)
            if (depth[site] < depth[best]) best = site;
        sites1[c] = best;
        depth[best] = 0.0f;

        int second = 0;
        for (int site = 1; site < nSites; ++site)
            if (depth[site] < depth[second]) second = site;
        sites2[c] = second;
        depth[second] = 0.0f;

        int third = 0;
        for (int site = 1; site < nSites; ++site)
            if (depth[site] < depth[third]) third = site;
        sites3[c] = third;
    }
}

}  // namespace

std::vector<float> computeWaveformSim(const MeanWaveforms& waveforms,
                                      const ClusterResult& clusters,
                                      const Probe& probe,
                                      const PostProcessSettings& settings)
{
    const int nClusters = waveforms.nClusters;
    std::vector<float> sim(static_cast<std::size_t>(nClusters) * nClusters, 0.0f);
    if (nClusters < 2) return sim;

    std::vector<int> sites1, sites2, sites3;
    secondaryPeaks(waveforms, sites1, sites2, sites3);

    // Shift up to one refractory period either way, so a unit whose waveform is
    // a few samples out of alignment still matches.
    const int nShifts = static_cast<int>(std::ceil(
        settings.meanInterpFactor * settings.refracInt * settings.sampleRate / 1000.0));
    const auto shifts = buildShiftRange(waveforms.nSamples, nShifts);

    const int nNeighbors = probe.nSiteNeighbors();

    parallelFor(static_cast<std::size_t>(nClusters), [&](std::size_t clusterIndex) {
        const int i = static_cast<int>(clusterIndex);
        const int iSite = sites1[i];
        const int* iNeighbors = probe.siteNeighbors(iSite);

        // Candidates: clusters peaking on this site (primary, secondary or
        // tertiary), or whose primary is one of this cluster's other peaks.
        // Only j < i, since the matrix is symmetrised afterwards.
        std::vector<int> candidates;
        for (int j = 0; j < i; ++j) {
            if (sites1[j] == iSite || sites2[j] == iSite || sites3[j] == iSite ||
                sites1[j] == sites2[i] || sites1[j] == sites3[i])
                candidates.push_back(j);
        }
        if (candidates.empty()) return;

        std::vector<float> x, y;

        for (int j : candidates) {
            const int jSite = sites1[j];

            // Sites shared by both neighbourhoods; identical sites means the
            // whole neighbourhood is comparable.
            std::vector<int> compareSites;
            if (jSite == iSite) {
                compareSites.assign(iNeighbors, iNeighbors + nNeighbors);
            } else {
                const int* jNeighbors = probe.siteNeighbors(jSite);
                for (int a = 0; a < nNeighbors; ++a) {
                    for (int b = 0; b < nNeighbors; ++b) {
                        if (iNeighbors[a] == jNeighbors[b]) {
                            compareSites.push_back(iNeighbors[a]);
                            break;
                        }
                    }
                }
            }
            if (compareSites.empty()) continue;

            double best = -1.0;
            for (const auto& shift : shifts) {
                x.clear();
                y.clear();
                for (int site : compareSites) {
                    const float* wi = waveforms.globalAt(i, site);
                    const float* wj = waveforms.globalAt(j, site);
                    for (std::size_t k = 0; k < shift.back.size(); ++k) {
                        x.push_back(wi[shift.back[k]]);
                        y.push_back(wj[shift.forward[k]]);
                    }
                }
                best = std::max(best, pearson(x, y));
            }

            if (best > 0)
                sim[static_cast<std::size_t>(i) * nClusters + j] =
                    static_cast<float>(best);
        }
    });

    // Symmetrise: waveformSim = max(S, S').
    for (int i = 0; i < nClusters; ++i) {
        for (int j = 0; j < i; ++j) {
            const float v = std::max(sim[static_cast<std::size_t>(i) * nClusters + j],
                                     sim[static_cast<std::size_t>(j) * nClusters + i]);
            sim[static_cast<std::size_t>(i) * nClusters + j] = v;
            sim[static_cast<std::size_t>(j) * nClusters + i] = v;
        }
    }

    return sim;
}

// ------------------------------------------------------------------- merging --

namespace {

/// Applies an old->new cluster mapping (1-based, 0 = drop) and compacts.
void relabelClusters(ClusterResult& clusters, const std::vector<std::int32_t>& mapping,
                     const detect::DetectResult& detected)
{
    const std::size_t nSpikes = clusters.spikeClusters.size();

    for (std::size_t i = 0; i < nSpikes; ++i) {
        const std::int32_t c = clusters.spikeClusters[i];
        clusters.spikeClusters[i] =
            (c > 0 && c < static_cast<std::int32_t>(mapping.size())) ? mapping[c] : 0;
    }

    std::int32_t nNew = 0;
    for (std::int32_t v : mapping) nNew = std::max(nNew, v);

    // Rebuild the per-cluster tables from the relabelled spikes.
    std::vector<std::int32_t> counts(nNew, 0);
    std::vector<std::vector<std::int32_t>> votes(nNew);
    for (std::size_t i = 0; i < nSpikes; ++i) {
        const std::int32_t c = clusters.spikeClusters[i];
        if (c <= 0) continue;
        ++counts[c - 1];
        votes[c - 1].push_back(detected.spikeSites[i]);
    }

    std::vector<std::int32_t> newCenters(nNew, 0);
    for (std::size_t old = 1; old < mapping.size(); ++old) {
        const std::int32_t neu = mapping[old];
        if (neu > 0 && newCenters[neu - 1] == 0 &&
            old - 1 < clusters.clusterCenters.size())
            newCenters[neu - 1] = clusters.clusterCenters[old - 1];
    }

    std::vector<std::int32_t> newSites(nNew, 0);
    for (std::int32_t c = 0; c < nNew; ++c) {
        if (votes[c].empty()) continue;
        std::sort(votes[c].begin(), votes[c].end());
        std::int32_t best = votes[c][0], bestRun = 0, current = votes[c][0], run = 0;
        for (std::int32_t s : votes[c]) {
            run = (s == current) ? run + 1 : 1;
            current = s;
            if (run > bestRun) {
                bestRun = run;
                best = current;
            }
        }
        newSites[c] = best;
    }

    clusters.nClusters = nNew;
    clusters.unitCount = std::move(counts);
    clusters.clusterSites = std::move(newSites);
    clusters.clusterCenters = std::move(newCenters);
}

}  // namespace

CoincidenceTest refractoryCoincidence(const std::vector<std::int64_t>& timesA,
                                      const std::vector<std::int64_t>& timesB,
                                      std::int64_t refracSamples,
                                      std::int64_t recordingSamples, double minExpected)
{
    CoincidenceTest test;
    if (timesA.empty() || timesB.empty() || recordingSamples <= 0) return test;

    // Sweep both sorted trains, counting B spikes inside each A spike's window.
    std::size_t lo = 0;
    for (std::int64_t t : timesA) {
        while (lo < timesB.size() && timesB[lo] < t - refracSamples) ++lo;
        for (std::size_t k = lo; k < timesB.size(); ++k) {
            if (timesB[k] > t + refracSamples) break;
            ++test.observed;
        }
    }

    // Independent trains would place this many pairs in the window.
    const double window = static_cast<double>(2 * refracSamples + 1);
    test.expected = static_cast<double>(timesA.size()) *
                    static_cast<double>(timesB.size()) * window /
                    static_cast<double>(recordingSamples);

    test.hasPower = test.expected >= minExpected;
    test.ratio = test.expected > 0
                     ? static_cast<double>(test.observed) / test.expected
                     : 1.0;
    return test;
}

int autoMerge(const detect::DetectResult& detected, ClusterResult& clusters,
              const Probe& probe, const PostProcessSettings& settings,
              PostProcessResult& out, const PostProgressFn& progress)
{
    int totalMerged = 0;

    const std::int64_t refracSamples = static_cast<std::int64_t>(
        std::lround(settings.refracInt * settings.sampleRate / 1000.0));
    std::int64_t recordingSamples = 1;
    if (!detected.spikeTimes.empty())
        recordingSamples = detected.spikeTimes.back() - detected.spikeTimes.front() + 1;

    // orderClusters('clusterSites') first, as autoMerge.m does.
    orderClustersBySite(clusters);

    for (int pass = 0; pass < settings.nPassesMerge; ++pass) {
        if (progress) progress("mean waveforms", pass + 1, settings.nPassesMerge);
        out.waveforms = computeMeanWaveforms(detected, clusters, probe, settings);

        if (progress) progress("waveform similarity", pass + 1, settings.nPassesMerge);
        out.waveformSim = computeWaveformSim(out.waveforms, clusters, probe, settings);

        const int n = clusters.nClusters;
        if (n < 2) break;

        // For each cluster j, the best-matching cluster i < j. Zeroing the
        // lower triangle (and diagonal) forces merges to run downwards, so a
        // unit always folds into a lower-numbered one.
        std::vector<float> bestScore(n, 0.0f);
        std::vector<int> bestTarget(n, -1);

        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < j; ++i) {
                const float v = out.waveformSim[static_cast<std::size_t>(j) * n + i];
                if (v > bestScore[j]) {
                    bestScore[j] = v;
                    bestTarget[j] = i;
                }
            }
        }

        std::vector<int> merging;
        std::vector<int> targets;
        for (int j = 0; j < n; ++j) {
            if (bestTarget[j] >= 0 && bestScore[j] >= settings.maxUnitSim) {
                merging.push_back(j);
                targets.push_back(bestTarget[j]);
            }
        }
        if (merging.empty()) break;

        // Waveform shape says these look alike; the spike trains say whether
        // they can be the same neuron. Drop the pairs that fire independently.
        if (settings.refractoryMergeCheck && refracSamples > 0) {
            const auto byCluster = spikesByCluster(clusters, detected.size());

            std::vector<int> keptMerging, keptTargets;
            for (std::size_t k = 0; k < merging.size(); ++k) {
                const auto collect = [&](int cluster) {
                    std::vector<std::int64_t> times;
                    times.reserve(byCluster[cluster].size());
                    for (std::int32_t spike : byCluster[cluster])
                        times.push_back(detected.spikeTimes[spike]);
                    return times;  // already ascending: spikes are time-ordered
                };

                const CoincidenceTest test = refractoryCoincidence(
                    collect(merging[k]), collect(targets[k]), refracSamples,
                    recordingSamples, settings.minExpectedCoincidences);

                // Only reject when the test actually has power; otherwise the
                // waveform decision stands.
                const bool veto =
                    test.hasPower && test.ratio > settings.maxMergeContamination;

                if (veto) {
                    ++out.nMergesVetoed;
                } else {
                    keptMerging.push_back(merging[k]);
                    keptTargets.push_back(targets[k]);
                }
            }

            merging = std::move(keptMerging);
            targets = std::move(keptTargets);
            if (merging.empty()) break;
        }

        // A unit that is itself a merge target must not also be merged away in
        // the same pass, or the chain would lose spikes.
        std::vector<char> isTarget(n, 0);
        for (int t : targets) isTarget[t] = 1;

        std::vector<int> keepMerging, keepTargets;
        for (std::size_t k = 0; k < merging.size(); ++k) {
            if (!isTarget[merging[k]]) {
                keepMerging.push_back(merging[k]);
                keepTargets.push_back(targets[k]);
            }
        }
        if (keepMerging.empty()) break;

        // Build the old(1-based) -> new(1-based) mapping.
        std::vector<std::int32_t> redirect(n);
        std::iota(redirect.begin(), redirect.end(), 0);
        for (std::size_t k = 0; k < keepMerging.size(); ++k)
            redirect[keepMerging[k]] = keepTargets[k];

        std::vector<char> survives(n, 1);
        for (int m : keepMerging) survives[m] = 0;

        std::vector<std::int32_t> compact(n, 0);
        std::int32_t next = 1;
        for (int c = 0; c < n; ++c)
            if (survives[c]) compact[c] = next++;

        std::vector<std::int32_t> mapping(n + 1, 0);
        for (int c = 0; c < n; ++c) mapping[c + 1] = compact[redirect[c]];

        relabelClusters(clusters, mapping, detected);

        totalMerged += static_cast<int>(keepMerging.size());
        out.nMergePasses = pass + 1;
    }

    out.nMergedAway = totalMerged;
    return totalMerged;
}

// ----------------------------------------------------------------- reordering --

void orderClustersBySite(ClusterResult& clusters)
{
    const int n = clusters.nClusters;
    if (n < 2) return;

    std::vector<int> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        if (clusters.clusterSites[a] != clusters.clusterSites[b])
            return clusters.clusterSites[a] < clusters.clusterSites[b];
        return a < b;
    });

    std::vector<std::int32_t> mapping(n + 1, 0);
    for (int newIndex = 0; newIndex < n; ++newIndex)
        mapping[order[newIndex] + 1] = newIndex + 1;

    for (auto& c : clusters.spikeClusters)
        if (c > 0 && c <= n) c = mapping[c];

    const auto reorder = [&](auto& vec) {
        using T = typename std::decay_t<decltype(vec)>::value_type;
        if (static_cast<int>(vec.size()) != n) return;
        std::vector<T> tmp(n);
        for (int newIndex = 0; newIndex < n; ++newIndex)
            tmp[newIndex] = vec[order[newIndex]];
        vec = std::move(tmp);
    };

    reorder(clusters.clusterSites);
    reorder(clusters.unitCount);
    reorder(clusters.clusterCenters);
}

// ------------------------------------------------------------- quality scores --

UnitMetrics computeQualityScores(const detect::DetectResult& detected,
                                 const ClusterResult& clusters,
                                 const MeanWaveforms& waveforms, const Probe& probe,
                                 const PostProcessSettings& settings)
{
    const int nClusters = clusters.nClusters;
    const int nSites = probe.nSites();

    UnitMetrics m;
    m.unitPeaks.assign(nClusters, 0.0f);
    m.unitPeakSites.assign(nClusters, 0);
    m.unitVpp.assign(nClusters, 0.0f);
    m.unitPeaksRaw.assign(nClusters, 0.0f);
    m.unitVppRaw.assign(nClusters, 0.0f);
    m.unitSNR.assign(nClusters, 0.0f);
    m.nSitesOverThresh.assign(nClusters, 0);
    m.unitISIRatio.assign(nClusters, std::numeric_limits<double>::quiet_NaN());
    m.unitIsoDist.assign(nClusters, std::numeric_limits<double>::quiet_NaN());
    m.unitLRatio.assign(nClusters, std::numeric_limits<double>::quiet_NaN());

    if (nClusters == 0) return m;

    // siteRMS = bit2uV(siteThresh / qqFactor)
    m.siteRMS.assign(nSites, 0.0f);
    const double uvScale = settings.bitScaling / settings.filterNorm;
    for (int s = 0; s < nSites && s < static_cast<int>(detected.siteThresholds.size());
         ++s)
        m.siteRMS[s] =
            static_cast<float>(detected.siteThresholds[s] / settings.qqFactor * uvScale);

    const bool haveRaw = !waveforms.globalRaw.empty();

    for (int c = 0; c < nClusters; ++c) {
        const int site = clusters.clusterSites[c];

        // Peak and Vpp on the cluster's own site.
        const float* w = waveforms.globalAt(c, site);
        float lo = w[0], hi = w[0];
        for (int s = 1; s < waveforms.nSamples; ++s) {
            lo = std::min(lo, w[s]);
            hi = std::max(hi, w[s]);
        }
        m.unitPeaks[c] = std::fabs(lo);
        m.unitVpp[c] = hi - lo;

        if (haveRaw) {
            const float* r = waveforms.globalRawAt(c, site);
            float rlo = r[0], rhi = r[0];
            for (int s = 1; s < waveforms.nSamplesRaw; ++s) {
                rlo = std::min(rlo, r[s]);
                rhi = std::max(rhi, r[s]);
            }
            m.unitPeaksRaw[c] = rlo;
            m.unitVppRaw[c] = rhi - rlo;
        }

        // Deepest local site, and how many sites the unit clears threshold on.
        int bestLocal = 0;
        float bestDepth = 0.0f;
        for (int k = 0; k < waveforms.nLocalSites; ++k) {
            const float* lw = waveforms.local.data() +
                              ((static_cast<std::size_t>(c) * waveforms.nLocalSites) + k) *
                                  waveforms.nSamples;
            float d = lw[0];
            for (int s = 1; s < waveforms.nSamples; ++s) d = std::min(d, lw[s]);
            if (d < bestDepth) {
                bestDepth = d;
                bestLocal = k;
            }
        }
        m.unitPeakSites[c] = bestLocal;

        int over = 0;
        for (int s = 0; s < nSites; ++s) {
            const float* gw = waveforms.globalAt(c, s);
            float d = gw[0];
            for (int k = 1; k < waveforms.nSamples; ++k) d = std::min(d, gw[k]);
            if (d < -m.siteRMS[s] * settings.qqFactor) ++over;
        }
        m.nSitesOverThresh[c] = over;

        if (site < static_cast<int>(m.siteRMS.size()) && m.siteRMS[site] > 0)
            m.unitSNR[c] = m.unitPeaks[c] / m.siteRMS[site];
    }

    // -- ISI ratio ---------------------------------------------------------
    const auto byCluster = spikesByCluster(clusters, detected.size());
    const std::int64_t samples2ms =
        static_cast<std::int64_t>(std::lround(settings.sampleRate * 0.002));
    const std::int64_t samples20ms =
        static_cast<std::int64_t>(std::lround(settings.sampleRate * 0.02));

    for (int c = 0; c < nClusters; ++c) {
        const auto& spikes = byCluster[c];
        if (spikes.size() < 2) continue;

        std::int64_t nShort = 0, nLong = 0;
        for (std::size_t i = 1; i < spikes.size(); ++i) {
            const std::int64_t d =
                detected.spikeTimes[spikes[i]] - detected.spikeTimes[spikes[i - 1]];
            if (d <= samples2ms) ++nShort;
            if (d <= samples20ms) ++nLong;
        }
        if (nLong > 0)
            m.unitISIRatio[c] = static_cast<double>(nShort) / static_cast<double>(nLong);
    }

    // -- isolation distance and L-ratio ------------------------------------
    // Mahalanobis distance of every spike on the unit's site from that unit's
    // feature distribution, exactly as computeQualityScores.m does.
    const int nFeatures = detected.nFeatures;
    if (nFeatures <= 0) return m;

    std::vector<std::vector<std::int32_t>> siteSpikes(nSites);
    std::vector<std::vector<int>> sitePeak(nSites);
    for (std::size_t i = 0; i < detected.size(); ++i) {
        siteSpikes[detected.spikeSites[i]].push_back(static_cast<std::int32_t>(i));
        sitePeak[detected.spikeSites[i]].push_back(0);
    }
    if (detected.nPeaks >= 2 && detected.spikeSites2.size() == detected.size()) {
        for (std::size_t i = 0; i < detected.size(); ++i) {
            siteSpikes[detected.spikeSites2[i]].push_back(static_cast<std::int32_t>(i));
            sitePeak[detected.spikeSites2[i]].push_back(1);
        }
    }

    parallelFor(static_cast<std::size_t>(nClusters), [&](std::size_t clusterIndex) {
        const int c = static_cast<int>(clusterIndex);
        const int site = clusters.clusterSites[c];

        const auto& spikes = siteSpikes[site];
        const auto& peaks = sitePeak[site];
        const int n = static_cast<int>(spikes.size());
        if (n < nFeatures + 2) return;

        std::vector<char> inCluster(n, 0);
        int nIn = 0;
        for (int i = 0; i < n; ++i) {
            if (clusters.spikeClusters[spikes[i]] == c + 1) {
                inCluster[i] = 1;
                ++nIn;
            }
        }
        if (nIn < nFeatures + 2 || nIn == n) return;

        // Mean and covariance of the cluster's own features.
        std::vector<double> mean(nFeatures, 0.0);
        for (int i = 0; i < n; ++i) {
            if (!inCluster[i]) continue;
            const float* f = detected.featuresFor(spikes[i], peaks[i]);
            for (int k = 0; k < nFeatures; ++k) mean[k] += f[k];
        }
        for (double& v : mean) v /= nIn;

        std::vector<double> cov(static_cast<std::size_t>(nFeatures) * nFeatures, 0.0);
        for (int i = 0; i < n; ++i) {
            if (!inCluster[i]) continue;
            const float* f = detected.featuresFor(spikes[i], peaks[i]);
            for (int a = 0; a < nFeatures; ++a) {
                const double da = f[a] - mean[a];
                for (int b = 0; b <= a; ++b)
                    cov[static_cast<std::size_t>(a) * nFeatures + b] +=
                        da * (f[b] - mean[b]);
            }
        }
        const double norm = 1.0 / (nIn - 1);
        for (int a = 0; a < nFeatures; ++a)
            for (int b = 0; b <= a; ++b) {
                cov[static_cast<std::size_t>(a) * nFeatures + b] *= norm;
                cov[static_cast<std::size_t>(b) * nFeatures + a] =
                    cov[static_cast<std::size_t>(a) * nFeatures + b];
            }

        // A singular covariance means the unit's features are degenerate;
        // MATLAB bails out here too rather than reporting a bogus distance.
        std::vector<double> chol = cov;
        if (!cholesky(chol, nFeatures)) return;

        std::vector<double> outside;
        outside.reserve(n - nIn);
        std::vector<double> z(nFeatures);

        for (int i = 0; i < n; ++i) {
            if (inCluster[i]) continue;
            const float* f = detected.featuresFor(spikes[i], peaks[i]);

            // Forward-substitute L z = (x - mu); then d^2 = z'z.
            for (int a = 0; a < nFeatures; ++a) {
                double sum = f[a] - mean[a];
                for (int b = 0; b < a; ++b)
                    sum -= chol[static_cast<std::size_t>(a) * nFeatures + b] * z[b];
                z[a] = sum / chol[static_cast<std::size_t>(a) * nFeatures + a];
            }
            double d2 = 0.0;
            for (int a = 0; a < nFeatures; ++a) d2 += z[a] * z[a];
            outside.push_back(d2);
        }

        if (outside.empty()) return;

        // L-ratio: how much of the chi-squared tail the outsiders occupy.
        //
        // The degrees of freedom are the feature dimensionality. MATLAB passes
        // size(siteFeatures,2) here, which is the SPIKE COUNT - with dof in the
        // thousands the tail probability saturates at 1 and the statistic
        // degenerates to nOutside/nIn. matlabLRatioDof reproduces that.
        const double dof = settings.matlabLRatioDof ? static_cast<double>(n)
                                                    : static_cast<double>(nFeatures);
        double lsum = 0.0;
        for (double d2 : outside) lsum += 1.0 - chi2cdf(d2, dof);
        m.unitLRatio[c] = lsum / nIn;

        // Isolation distance is undefined once the unit owns most of the site.
        if (static_cast<double>(nIn) / n > 0.5) return;
        if (static_cast<int>(outside.size()) < nIn) return;

        std::nth_element(outside.begin(), outside.begin() + (nIn - 1), outside.end());
        m.unitIsoDist[c] = outside[nIn - 1];
    });

    return m;
}

// ------------------------------------------------------- refractory cleanup --

std::int64_t removeOutlierSpikes(const detect::DetectResult& detected,
                                 ClusterResult& clusters,
                                 const PostProcessSettings& settings)
{
    if (settings.outlierThresh <= 0) return 0;

    const int nFeatures = detected.nFeatures;
    if (nFeatures <= 0 || detected.spikeFeatures.empty()) return 0;

    const auto byCluster = spikesByCluster(clusters, detected.size());
    const bool haveSecondary =
        detected.nPeaks >= 2 && detected.spikeSites2.size() == detected.size();

    std::vector<std::vector<std::int32_t>> toDrop(clusters.nClusters);

    parallelFor(static_cast<std::size_t>(clusters.nClusters), [&](std::size_t index) {
        const int c = static_cast<int>(index);
        const auto& spikes = byCluster[c];
        const int n = static_cast<int>(spikes.size());

        // Need more spikes than dimensions for the covariance to be invertible.
        if (n < nFeatures + 2) return;

        const int site = clusters.clusterSites[c];

        // Where this cluster's site is a spike's SECONDARY site, that peak's
        // features describe it better, so prefer them - as rmOutlierSpikes.m does.
        const auto featuresFor = [&](std::int32_t spike) {
            const int peak =
                (haveSecondary && detected.spikeSites2[spike] == site) ? 1 : 0;
            return detected.featuresFor(spike, peak);
        };

        std::vector<double> mean(nFeatures, 0.0);
        for (std::int32_t spike : spikes) {
            const float* f = featuresFor(spike);
            for (int k = 0; k < nFeatures; ++k) mean[k] += f[k];
        }
        for (double& v : mean) v /= n;

        std::vector<double> cov(static_cast<std::size_t>(nFeatures) * nFeatures, 0.0);
        for (std::int32_t spike : spikes) {
            const float* f = featuresFor(spike);
            for (int a = 0; a < nFeatures; ++a) {
                const double da = f[a] - mean[a];
                for (int b = 0; b <= a; ++b)
                    cov[static_cast<std::size_t>(a) * nFeatures + b] +=
                        da * (f[b] - mean[b]);
            }
        }
        const double norm = 1.0 / (n - 1);
        for (int a = 0; a < nFeatures; ++a)
            for (int b = 0; b <= a; ++b) {
                cov[static_cast<std::size_t>(a) * nFeatures + b] *= norm;
                cov[static_cast<std::size_t>(b) * nFeatures + a] =
                    cov[static_cast<std::size_t>(a) * nFeatures + b];
            }

        // A degenerate covariance means the distance is meaningless; MATLAB
        // skips the cluster rather than guessing, and so do we.
        std::vector<double> chol = cov;
        if (!cholesky(chol, nFeatures)) return;

        std::vector<double> logDist(n);
        std::vector<double> z(nFeatures);

        for (int i = 0; i < n; ++i) {
            const float* f = featuresFor(spikes[i]);
            for (int a = 0; a < nFeatures; ++a) {
                double sum = f[a] - mean[a];
                for (int b = 0; b < a; ++b)
                    sum -= chol[static_cast<std::size_t>(a) * nFeatures + b] * z[b];
                z[a] = sum / chol[static_cast<std::size_t>(a) * nFeatures + a];
            }
            double d2 = 0.0;
            for (int a = 0; a < nFeatures; ++a) d2 += z[a] * z[a];
            logDist[i] = std::log(std::max(d2, 1e-300));
        }

        // madScore: (x - median) / median(|x - median|). No 1.4826 scaling,
        // matching madScore.m, so the threshold means the same thing.
        std::vector<double> sorted = logDist;
        std::sort(sorted.begin(), sorted.end());
        const double median = sorted[sorted.size() / 2];

        std::vector<double> deviation(n);
        for (int i = 0; i < n; ++i) deviation[i] = std::fabs(logDist[i] - median);
        std::sort(deviation.begin(), deviation.end());
        const double mad = deviation[deviation.size() / 2];
        if (mad <= 0) return;

        for (int i = 0; i < n; ++i)
            if ((logDist[i] - median) / mad > settings.outlierThresh)
                toDrop[c].push_back(spikes[i]);
    });

    std::int64_t removed = 0;
    for (int c = 0; c < clusters.nClusters; ++c) {
        for (std::int32_t spike : toDrop[c]) {
            clusters.spikeClusters[spike] = 0;  // belongs to no unit
            ++removed;
        }
    }

    if (removed > 0) {
        std::vector<std::int32_t> identity(clusters.nClusters + 1, 0);
        for (int c = 1; c <= clusters.nClusters; ++c) identity[c] = c;
        relabelClusters(clusters, identity, detected);
    }

    return removed;
}

std::int64_t removeRefractorySpikes(const detect::DetectResult& detected,
                                    ClusterResult& clusters,
                                    const PostProcessSettings& settings)
{
    const std::int64_t refracSamples = static_cast<std::int64_t>(
        std::lround(settings.refracInt * settings.sampleRate / 1000.0));
    if (refracSamples <= 0) return 0;

    const auto byCluster = spikesByCluster(clusters, detected.size());
    std::int64_t removed = 0;

    for (int c = 0; c < clusters.nClusters; ++c) {
        const auto& spikes = byCluster[c];
        if (spikes.size() < 2) continue;

        // Spikes are time-ordered, so a single sweep finds every violation.
        // Of each offending pair the smaller-amplitude spike is dropped.
        std::int32_t previous = spikes[0];
        for (std::size_t i = 1; i < spikes.size(); ++i) {
            const std::int32_t current = spikes[i];
            if (clusters.spikeClusters[previous] != c + 1) {
                previous = current;
                continue;
            }

            const std::int64_t gap =
                detected.spikeTimes[current] - detected.spikeTimes[previous];
            if (gap > refracSamples) {
                previous = current;
                continue;
            }

            const float ampPrev = std::fabs(detected.spikeAmps[previous]);
            const float ampCurr = std::fabs(detected.spikeAmps[current]);

            if (ampCurr > ampPrev) {
                clusters.spikeClusters[previous] = 0;
                previous = current;
            } else {
                clusters.spikeClusters[current] = 0;
            }
            ++removed;
        }
    }

    if (removed > 0) {
        // Counts and sites shift when spikes are dropped.
        std::vector<std::int32_t> identity(clusters.nClusters + 1, 0);
        for (int c = 1; c <= clusters.nClusters; ++c) identity[c] = c;
        relabelClusters(clusters, identity, detected);
    }

    return removed;
}

// ------------------------------------------------------------------ pipeline --

PostProcessResult runPostProcessing(const detect::DetectResult& detected,
                                    ClusterResult& clusters, const Probe& probe,
                                    const PostProcessSettings& settings,
                                    bool doAutoMerge, bool doRemoveRefractory,
                                    bool doRemoveOutliers,
                                    const PostProgressFn& progress)
{
    const auto start = std::chrono::steady_clock::now();
    PostProcessResult out;

    if (detected.filtWaveforms.empty())
        throw std::runtime_error(
            "post-processing needs spike waveforms; enable retainWaveforms");

    // Before merging: a cluster carrying a burst of artifacts has a distorted
    // mean waveform, which would steer the merge decisions.
    if (doRemoveOutliers) {
        if (progress) progress("outlier spike removal", 0, 0);
        out.nOutliersRemoved = removeOutlierSpikes(detected, clusters, settings);
    }

    if (doRemoveRefractory) {
        if (progress) progress("refractory cleanup", 0, 0);
        out.nRefractoryRemoved =
            removeRefractorySpikes(detected, clusters, settings);
    }

    if (doAutoMerge) {
        autoMerge(detected, clusters, probe, settings, out, progress);
    } else {
        orderClustersBySite(clusters);
    }

    // Recompute against the final labelling, so what is reported matches what
    // is saved even when the last merge pass changed the clusters.
    if (progress) progress("mean waveforms", 0, 0);
    out.waveforms = computeMeanWaveforms(detected, clusters, probe, settings);

    if (progress) progress("waveform similarity", 0, 0);
    out.waveformSim = computeWaveformSim(out.waveforms, clusters, probe, settings);

    if (progress) progress("quality scores", 0, 0);
    out.metrics =
        computeQualityScores(detected, clusters, out.waveforms, probe, settings);

    out.seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    return out;
}

}  // namespace sort
}  // namespace jrc
