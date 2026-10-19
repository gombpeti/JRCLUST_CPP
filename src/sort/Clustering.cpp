#include "sort/Clustering.h"

#include "core/Parallel.h"
#include "core/Stats.h"
#include "cuda/DensityPeak.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace jrc {
namespace sort {

ClusterSettings ClusterSettings::fromConfig(const Config& config)
{
    ClusterSettings settings;
    settings.distCut = config.getNumber("distCut", 2.0);
    settings.nClusterIntervals = config.getInt("nClusterIntervals", 4);
    settings.log10RhoCut = config.getNumber("log10RhoCut", -2.5);
    settings.log10DeltaCut = config.getNumber("log10DeltaCut", 0.6);
    settings.minClusterSize = config.getInt("minClusterSize", 30);
    settings.maxClustersSite = config.getInt("maxClustersSite", 20);
    return settings;
}

namespace {

/// Rank of each element in ascending order; ties broken by position.
/// Port of jrclust.utils.rankorder.
std::vector<std::int32_t> rankOrderAscending(const std::vector<std::int64_t>& values)
{
    const std::size_t n = values.size();
    std::vector<std::size_t> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(),
              [&](std::size_t a, std::size_t b) { return values[a] < values[b]; });

    std::vector<std::int32_t> rank(n);
    for (std::size_t i = 0; i < n; ++i)
        rank[order[i]] = static_cast<std::int32_t>(i);
    return rank;
}

/// Squared euclidean distance between two feature vectors.
inline float dist2(const float* a, const float* b, int n)
{
    float acc = 0.0f;
    for (int i = 0; i < n; ++i) {
        const float d = a[i] - b[i];
        acc += d * d;
    }
    return acc;
}

/// Per-site distance cutoff.
///
/// Port of estRhoCutSite() in computeRho.m: subsample 1000 primary and 4000
/// comparison spikes, take each primary spike's distCut-percentile distance to
/// the comparison set (restricted to spikes nearby in time), then the median
/// across primaries.
float estimateRhoCut(const std::vector<float>& siteFeatures,
                     const std::vector<std::int32_t>& spikeOrder, int n1, int n2,
                     int nFeatures, const ClusterSettings& settings)
{
    if (n1 < 2) return 1.0f;

    constexpr int kSubsample = 1000;
    const int nPrimary = std::min(n1, kSubsample);
    const int nSecondary = std::min(n1, 4 * kSubsample);

    const int primaryStride = std::max(1, n1 / nPrimary);
    const int secondaryStride = std::max(1, n1 / nSecondary);

    const double timeWindow =
        static_cast<double>(n1 + n2) / std::max(1, settings.nClusterIntervals);

    std::vector<float> perPrimary;
    perPrimary.reserve(static_cast<std::size_t>(nPrimary));

    for (int i = 0; i < n1; i += primaryStride) {
        const float* fi = siteFeatures.data() + static_cast<std::size_t>(i) * nFeatures;

        std::vector<float> dists;
        dists.reserve(static_cast<std::size_t>(nSecondary));

        for (int j = 0; j < n1; j += secondaryStride) {
            if (std::abs(spikeOrder[j] - spikeOrder[i]) >= timeWindow) continue;

            const float d =
                dist2(siteFeatures.data() + static_cast<std::size_t>(j) * nFeatures,
                      fi, nFeatures);
            if (d > 0.0f) dists.push_back(d);  // MATLAB drops <= 0 as nan
        }

        if (dists.empty()) continue;
        const float q = matlabQuantile(dists, settings.distCut / 100.0);
        if (std::isfinite(q)) perPrimary.push_back(q);
    }

    if (perPrimary.empty()) return 1.0f;

    // nanmedian across the sampled primary spikes.
    const std::size_t mid = perPrimary.size() / 2;
    std::nth_element(perPrimary.begin(), perPrimary.begin() + mid, perPrimary.end());
    float cut = perPrimary[mid];
    if (cut <= 0.0f) cut = 1.0f;
    return cut;
}

/// Least-squares fit of y ~ a*x^2 + b*x + c over `indices`, returning the
/// z-scored residual for every element. Port of detrendQuadratic().
std::vector<double> detrendQuadratic(const std::vector<double>& x,
                                     const std::vector<double>& y,
                                     const std::vector<std::size_t>& indices)
{
    const std::size_t n = x.size();
    std::vector<double> z(n, std::numeric_limits<double>::quiet_NaN());
    if (indices.size() < 3) return z;

    // Normal equations for the 3-parameter fit: X'X b = X'y, X = [x^2, x, 1].
    double s[5] = {0, 0, 0, 0, 0};  // sums of x^0..x^4
    double t[3] = {0, 0, 0};        // sums of y*x^0..x^2

    for (std::size_t i : indices) {
        const double xi = x[i];
        const double x2 = xi * xi;
        s[0] += 1.0;
        s[1] += xi;
        s[2] += x2;
        s[3] += x2 * xi;
        s[4] += x2 * x2;
        t[0] += y[i];
        t[1] += y[i] * xi;
        t[2] += y[i] * x2;
    }

    // Solve the symmetric 3x3 system by Cramer's rule.
    const double m[3][3] = {{s[4], s[3], s[2]}, {s[3], s[2], s[1]}, {s[2], s[1], s[0]}};
    const double rhs[3] = {t[2], t[1], t[0]};

    const auto det3 = [](const double a[3][3]) {
        return a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) -
               a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]) +
               a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
    };

    const double base = det3(m);
    if (std::fabs(base) < 1e-30) return z;

    double b[3];
    for (int k = 0; k < 3; ++k) {
        double mk[3][3];
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c) mk[r][c] = (c == k) ? rhs[r] : m[r][c];
        b[k] = det3(mk) / base;
    }

    std::vector<double> residual(n);
    for (std::size_t i = 0; i < n; ++i)
        residual[i] = y[i] - (b[0] * x[i] * x[i] + b[1] * x[i] + b[2]);

    // z-score against the fitted subset only, as MATLAB does.
    double mean = 0.0;
    for (std::size_t i : indices) mean += residual[i];
    mean /= static_cast<double>(indices.size());

    double variance = 0.0;
    for (std::size_t i : indices) {
        const double d = residual[i] - mean;
        variance += d * d;
    }
    double sd = std::sqrt(variance / std::max<std::size_t>(1, indices.size() - 1));
    if (sd == 0.0) sd = 1.0;

    for (std::size_t i = 0; i < n; ++i) z[i] = (residual[i] - mean) / sd;
    return z;
}

}  // namespace

ClusterResult cluster(const detect::DetectResult& detected, const Probe& probe,
                      const ClusterSettings& settings,
                      const ClusterProgressFn& progress)
{
    const std::size_t nSpikes = detected.size();
    const int nSites = probe.nSites();
    const int nFeatures = detected.nFeatures;

    if (nFeatures <= 0)
        throw std::runtime_error("clustering needs features; run detect with them on");
    if (nFeatures > cuda::kMaxFeatures)
        throw std::runtime_error(
            "feature vector is " + std::to_string(nFeatures) +
            " long but the CUDA kernels cap at " + std::to_string(cuda::kMaxFeatures) +
            " (MAXDIM); reduce nPCsPerSite or nSiteDir");

    ClusterResult result;
    result.spikeRho.assign(nSpikes, 0.0f);
    result.spikeDelta.assign(nSpikes, 0.0f);
    result.spikeNeigh.assign(nSpikes, -1);
    result.siteRhoCut.assign(nSites, 0.0f);

    const auto tStart = std::chrono::steady_clock::now();

    // -- group spikes by site ----------------------------------------------
    std::vector<std::vector<std::int32_t>> primary(nSites);
    std::vector<std::vector<std::int32_t>> secondary(nSites);

    for (std::size_t i = 0; i < nSpikes; ++i)
        primary[detected.spikeSites[i]].push_back(static_cast<std::int32_t>(i));

    const bool hasSecondary =
        detected.nPeaks >= 2 && detected.spikeSites2.size() == nSpikes;
    if (hasSecondary) {
        for (std::size_t i = 0; i < nSpikes; ++i)
            secondary[detected.spikeSites2[i]].push_back(static_cast<std::int32_t>(i));
    }

    // -- per site: rho, then delta -----------------------------------------
    double rhoSeconds = 0.0;
    double deltaSeconds = 0.0;

    for (int site = 0; site < nSites; ++site) {
        const int n1 = static_cast<int>(primary[site].size());
        const int n2 = static_cast<int>(secondary[site].size());
        const int n12 = n1 + n2;
        if (n1 < 2) {
            if (progress) progress(site + 1, nSites, 0);
            continue;
        }

        // Feature matrix for this site: primary spikes' peak-1 features
        // followed by secondary spikes' peak-2 features, as getSiteFeatures.m
        // assembles them. Column-major (feature fastest) for the kernels.
        std::vector<float> siteFeatures(static_cast<std::size_t>(nFeatures) * n12);
        std::vector<std::int64_t> spikeIds(n12);

        for (int k = 0; k < n1; ++k) {
            const std::int32_t spike = primary[site][k];
            std::copy_n(detected.featuresFor(spike, 0), nFeatures,
                        siteFeatures.begin() + static_cast<std::size_t>(k) * nFeatures);
            spikeIds[k] = spike;
        }
        for (int k = 0; k < n2; ++k) {
            const std::int32_t spike = secondary[site][k];
            std::copy_n(detected.featuresFor(spike, 1), nFeatures,
                        siteFeatures.begin() +
                            static_cast<std::size_t>(n1 + k) * nFeatures);
            spikeIds[n1 + k] = spike;
        }

        const std::vector<std::int32_t> spikeOrder = rankOrderAscending(spikeIds);

        cuda::SiteConstants consts;
        consts.n1 = n1;
        consts.n12 = n12;
        consts.nFeatures = nFeatures;
        consts.dnMax = static_cast<int>(std::lround(
            static_cast<double>(n12) / std::max(1, settings.nClusterIntervals)));
        consts.fDcSpk = 0;

        const float rhoCut =
            estimateRhoCut(siteFeatures, spikeOrder, n1, n2, nFeatures, settings);
        result.siteRhoCut[site] = rhoCut;

        // -- rho ------------------------------------------------------------
        std::vector<float> rho(n1, 0.0f);
        {
            const auto t0 = std::chrono::steady_clock::now();
            if (settings.useGpu) {
                cuda::computeRhoSite(siteFeatures.data(), spikeOrder.data(), consts,
                                     rhoCut, rho.data());
            } else {
                cuda::computeRhoSiteCpu(siteFeatures.data(), spikeOrder.data(), consts,
                                        rhoCut, rho.data());
            }
            rhoSeconds +=
                std::chrono::duration<double>(std::chrono::steady_clock::now() - t0)
                    .count();
        }

        for (int k = 0; k < n1; ++k) result.spikeRho[primary[site][k]] = rho[k];

        // -- delta ----------------------------------------------------------
        // rhoOrder ranks all n12 spikes by rho descending; the kernel uses it
        // to decide which neighbours are denser.
        std::vector<float> rhoAll(n12);
        for (int k = 0; k < n1; ++k) rhoAll[k] = rho[k];
        for (int k = 0; k < n2; ++k)
            rhoAll[n1 + k] = result.spikeRho[secondary[site][k]];

        std::vector<std::size_t> order(n12);
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
            return rhoAll[a] > rhoAll[b];  // descending
        });
        std::vector<std::int32_t> rhoOrder(n12);
        for (std::size_t k = 0; k < order.size(); ++k)
            rhoOrder[order[k]] = static_cast<std::int32_t>(k);

        std::vector<float> delta(n1, 0.0f);
        std::vector<std::int32_t> neigh(n1, 0);
        {
            const auto t0 = std::chrono::steady_clock::now();
            if (settings.useGpu) {
                cuda::computeDeltaSite(siteFeatures.data(), spikeOrder.data(),
                                       rhoOrder.data(), consts, rhoCut, delta.data(),
                                       neigh.data());
            } else {
                cuda::computeDeltaSiteCpu(siteFeatures.data(), spikeOrder.data(),
                                          rhoOrder.data(), consts, rhoCut, delta.data(),
                                          neigh.data());
            }
            deltaSeconds +=
                std::chrono::duration<double>(std::chrono::steady_clock::now() - t0)
                    .count();
        }

        for (int k = 0; k < n1; ++k) {
            const std::int32_t spike = primary[site][k];
            result.spikeDelta[spike] = delta[k];

            // Kernel neighbours index into this site's local list; map back to
            // global spike ids.
            const std::int32_t local = neigh[k];
            result.spikeNeigh[spike] =
                (local >= 0 && local < n12) ? static_cast<std::int32_t>(spikeIds[local])
                                            : spike;
        }

        if (progress) progress(site + 1, nSites, 0);
    }

    result.timings.rho = rhoSeconds;
    result.timings.delta = deltaSeconds;

    // -- find cluster centers ----------------------------------------------
    const auto tAssign = std::chrono::steady_clock::now();

    std::vector<double> logRho(nSpikes);
    std::vector<double> delta(nSpikes);
    for (std::size_t i = 0; i < nSpikes; ++i) {
        logRho[i] = (result.spikeRho[i] > 0) ? std::log10(result.spikeRho[i])
                                             : -std::numeric_limits<double>::infinity();
        delta[i] = result.spikeDelta[i];
    }

    // Global detrend: fit delta against log-rho over the bulk of the cloud,
    // then treat the largest positive residuals as candidate centers.
    const double rhoFloor = std::pow(10.0, settings.log10RhoCut);
    std::vector<std::size_t> fitIndices;
    for (std::size_t i = 0; i < nSpikes; ++i) {
        if (delta[i] > 0 && delta[i] < 1 && result.spikeRho[i] > rhoFloor &&
            result.spikeRho[i] < 0.1 && std::isfinite(logRho[i]))
            fitIndices.push_back(i);
    }

    std::vector<double> z = detrendQuadratic(logRho, delta, fitIndices);
    for (std::size_t i = 0; i < nSpikes; ++i)
        if (delta[i] == 0) z[i] = std::numeric_limits<double>::quiet_NaN();

    std::vector<std::size_t> candidates;
    for (std::size_t i = 0; i < nSpikes; ++i)
        if (result.spikeRho[i] > rhoFloor && !std::isnan(z[i])) candidates.push_back(i);

    const std::size_t maxCenters =
        static_cast<std::size_t>(settings.maxClustersSite) * nSites;

    std::sort(candidates.begin(), candidates.end(),
              [&](std::size_t a, std::size_t b) { return z[a] > z[b]; });
    if (candidates.size() > maxCenters) candidates.resize(maxCenters);

    const double deltaFloor = std::pow(10.0, settings.log10DeltaCut);
    std::vector<std::int32_t> centers;
    for (std::size_t i : candidates)
        if (z[i] >= deltaFloor) centers.push_back(static_cast<std::int32_t>(i));

    // -- assign spikes to clusters -----------------------------------------
    // Walk spikes from densest to sparsest: a spike's nearest denser neighbour
    // is always assigned before it, so one pass suffices. MATLAB iterates the
    // same relation to a fixed point.
    result.spikeClusters.assign(nSpikes, 0);
    for (std::size_t k = 0; k < centers.size(); ++k)
        result.spikeClusters[centers[k]] = static_cast<std::int32_t>(k + 1);

    std::vector<std::size_t> byRho(nSpikes);
    std::iota(byRho.begin(), byRho.end(), 0);
    std::sort(byRho.begin(), byRho.end(), [&](std::size_t a, std::size_t b) {
        return result.spikeRho[a] > result.spikeRho[b];
    });

    for (std::size_t spike : byRho) {
        if (result.spikeClusters[spike] > 0) continue;
        const std::int32_t neighbor = result.spikeNeigh[spike];
        if (neighbor < 0 || static_cast<std::size_t>(neighbor) >= nSpikes) continue;
        if (static_cast<std::size_t>(neighbor) == spike) continue;
        result.spikeClusters[spike] = result.spikeClusters[neighbor];
    }

    // -- drop undersized clusters and compact labels ------------------------
    // minClusterSize is raised to twice the feature count, as removeSmall does.
    const int minSize = std::max(settings.minClusterSize, 2 * nFeatures);

    std::vector<std::int32_t> counts(centers.size() + 1, 0);
    for (std::size_t i = 0; i < nSpikes; ++i)
        if (result.spikeClusters[i] > 0) ++counts[result.spikeClusters[i]];

    std::vector<std::int32_t> remap(centers.size() + 1, 0);
    std::int32_t next = 1;
    for (std::size_t c = 1; c < counts.size(); ++c)
        if (counts[c] >= minSize) remap[c] = next++;

    for (std::size_t i = 0; i < nSpikes; ++i)
        result.spikeClusters[i] = remap[result.spikeClusters[i]];

    result.nClusters = next - 1;

    result.clusterCenters.assign(result.nClusters, 0);
    for (std::size_t c = 1; c < remap.size(); ++c)
        if (remap[c] > 0) result.clusterCenters[remap[c] - 1] = centers[c - 1];

    // Modal site and spike count per surviving cluster.
    result.unitCount.assign(result.nClusters, 0);
    result.clusterSites.assign(result.nClusters, 0);

    std::vector<std::vector<std::int32_t>> siteVotes(result.nClusters);
    for (std::size_t i = 0; i < nSpikes; ++i) {
        const std::int32_t c = result.spikeClusters[i];
        if (c <= 0) continue;
        ++result.unitCount[c - 1];
        siteVotes[c - 1].push_back(detected.spikeSites[i]);
    }
    for (int c = 0; c < result.nClusters; ++c) {
        if (siteVotes[c].empty()) continue;
        std::sort(siteVotes[c].begin(), siteVotes[c].end());
        std::int32_t best = siteVotes[c][0], bestRun = 0;
        std::int32_t current = siteVotes[c][0], run = 0;
        for (std::int32_t s : siteVotes[c]) {
            if (s == current) {
                ++run;
            } else {
                current = s;
                run = 1;
            }
            if (run > bestRun) {
                bestRun = run;
                best = current;
            }
        }
        result.clusterSites[c] = best;
    }

    result.timings.assign =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - tAssign).count();
    result.timings.total =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - tStart).count();

    return result;
}

}  // namespace sort
}  // namespace jrc
