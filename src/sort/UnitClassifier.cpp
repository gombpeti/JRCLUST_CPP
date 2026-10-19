#include "sort/UnitClassifier.h"

#include "core/Parallel.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

namespace jrc {
namespace sort {

const char* qualityName(UnitQuality quality)
{
    switch (quality) {
        case UnitQuality::Noise: return "noise";
        case UnitQuality::Single: return "single";
        case UnitQuality::Multi: return "multi";
    }
    return "";
}

UnitShape measureShape(const MeanWaveforms& waveforms, int cluster, int peakSite,
                       double sampleRate, double spreadThreshold)
{
    UnitShape shape;
    shape.peakSite = peakSite;

    // Shape is measured on the RAW mean waveform. The filtered trace has been
    // through a differentiator, which sharpens everything and would make real
    // spikes look as narrow as glitches.
    const bool haveRaw = !waveforms.globalRaw.empty() && waveforms.nSamplesRaw > 2;
    const float* w = haveRaw ? waveforms.globalRawAt(cluster, peakSite)
                             : waveforms.globalAt(cluster, peakSite);
    const int n = haveRaw ? waveforms.nSamplesRaw : waveforms.nSamples;
    if (n < 3) return shape;

    const double msPerSample = 1000.0 / sampleRate;

    int trough = 0;
    for (int s = 1; s < n; ++s)
        if (w[s] < w[trough]) trough = s;

    double hi = w[0];
    for (int s = 1; s < n; ++s) hi = std::max<double>(hi, w[s]);

    shape.troughDepth = -static_cast<double>(w[trough]);
    shape.vpp = hi - w[trough];

    if (shape.troughDepth <= 0) return shape;  // no negative deflection at all

    // Width of the trough at half its depth: how long the waveform actually
    // spends below half-minimum. A one-sample glitch scores 1 sample; a real
    // spike several.
    const double half = w[trough] / 2.0;
    int left = trough;
    while (left > 0 && w[left - 1] <= half) --left;
    int right = trough;
    while (right + 1 < n && w[right + 1] <= half) ++right;
    shape.halfWidthMs = (right - left + 1) * msPerSample;

    // How much the trough stands out from its immediate neighbours. A membrane
    // potential cannot move and return within one sample, so a large value
    // means a digital glitch rather than a spike.
    const double before = std::fabs(trough > 0 ? w[trough - 1] : w[trough]);
    const double after = std::fabs(trough + 1 < n ? w[trough + 1] : w[trough]);
    const double neighbourhood = 0.5 * (before + after);
    shape.spikiness = neighbourhood > 1e-9 ? shape.troughDepth / neighbourhood
                                           : std::numeric_limits<double>::infinity();

    // A real spike repolarises: the trough is followed by a positive rebound.
    int peakAfter = trough;
    for (int s = trough; s < n; ++s)
        if (w[s] > w[peakAfter]) peakAfter = s;
    shape.troughToPeakMs = (peakAfter - trough) * msPerSample;
    shape.recoveryRatio = w[peakAfter] > 0 ? w[peakAfter] / shape.troughDepth : 0.0;

    // Spatial footprint: a real spike appears on a handful of neighbouring
    // sites with a smooth decay. One site suggests a local glitch; most of the
    // probe suggests a common-mode artifact.
    const double cut = spreadThreshold * shape.troughDepth;
    shape.spatialSpread = 0;
    for (int site = 0; site < waveforms.nSites; ++site) {
        const float* sw = haveRaw ? waveforms.globalRawAt(cluster, site)
                                  : waveforms.globalAt(cluster, site);
        double depth = 0;
        for (int s = 0; s < n; ++s) depth = std::min<double>(depth, sw[s]);
        if (-depth >= cut) ++shape.spatialSpread;
    }

    return shape;
}

void measureShankSpread(const MeanWaveforms& waveforms, int cluster,
                        const Probe& probe, double uniformFraction,
                        UnitShape& shape)
{
    shape.shankSites = 0;
    shape.shankAbove = 0;
    shape.shankCoverage = 0;
    if (shape.troughDepth <= 0) return;

    const bool haveRaw = !waveforms.globalRaw.empty() && waveforms.nSamplesRaw > 2;
    const int n = haveRaw ? waveforms.nSamplesRaw : waveforms.nSamples;
    if (n < 3) return;

    const int shank = probe.shank(shape.peakSite);
    const double cut = uniformFraction * shape.troughDepth;

    for (int site = 0; site < waveforms.nSites; ++site) {
        if (probe.shank(site) != shank) continue;
        if (probe.isIgnored(site)) continue;  // zeroed, would dilute the share
        ++shape.shankSites;

        const float* sw = haveRaw ? waveforms.globalRawAt(cluster, site)
                                  : waveforms.globalAt(cluster, site);
        double depth = 0;
        for (int s = 0; s < n; ++s) depth = std::min<double>(depth, sw[s]);
        if (-depth >= cut) ++shape.shankAbove;
    }

    if (shape.shankSites > 0)
        shape.shankCoverage =
            static_cast<double>(shape.shankAbove) / shape.shankSites;
}

namespace {

/// Spike-train statistics that no waveform test can reach.
struct TrainStats {
    double isiRatio = std::numeric_limits<double>::quiet_NaN();
    double zeroLagRatio = std::numeric_limits<double>::quiet_NaN();
    double burstFraction = 0.0;
};

/// Measures one unit's spike train. `times` must be ascending, in seconds.
TrainStats measureTrain(const std::vector<double>& times,
                        const QualityCriteria& criteria)
{
    TrainStats out;
    const std::size_t n = times.size();
    if (n < 20) return out;

    // -- refractory violations -------------------------------------------
    std::int64_t short2 = 0, short20 = 0;
    for (std::size_t i = 1; i < n; ++i) {
        const double ms = (times[i] - times[i - 1]) * 1000.0;
        if (ms <= 2.0) ++short2;
        if (ms <= 20.0) ++short20;
    }
    if (short20 > 0)
        out.isiRatio = static_cast<double>(short2) / static_cast<double>(short20);

    // -- zero-lag autocorrelogram ----------------------------------------
    // Pairs within 2 ms against pairs in 10-25 ms, scaled by window width so
    // that 1.0 means "no refractory period". A neuron scores well below it.
    std::int64_t near = 0, far = 0;
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            const double ms = (times[j] - times[i]) * 1000.0;
            if (ms > 25.0) break;
            if (ms <= 2.0) ++near;
            else if (ms >= 10.0) ++far;
        }
    }
    if (far > 0)
        out.zeroLagRatio =
            (static_cast<double>(near) / static_cast<double>(far)) * (15.0 / 2.0);

    // -- rate bursts against the unit's own baseline ----------------------
    const double span = times.back() - times.front();
    if (span > criteria.burstBinSeconds) {
        const int nBins =
            static_cast<int>(std::ceil(span / criteria.burstBinSeconds)) + 1;
        std::vector<int> counts(nBins, 0);
        for (double t : times) {
            const int b = static_cast<int>((t - times.front()) /
                                           criteria.burstBinSeconds);
            if (b >= 0 && b < nBins) ++counts[b];
        }

        std::vector<int> nonEmpty;
        for (int c : counts)
            if (c > 0) nonEmpty.push_back(c);
        if (!nonEmpty.empty()) {
            std::sort(nonEmpty.begin(), nonEmpty.end());
            const double med = nonEmpty[nonEmpty.size() / 2];

            std::vector<double> dev;
            dev.reserve(nonEmpty.size());
            for (int c : nonEmpty) dev.push_back(std::fabs(c - med));
            std::sort(dev.begin(), dev.end());
            double mad = dev[dev.size() / 2] * 1.4826;
            if (mad <= 0) mad = 1.0;

            const double cutoff = std::max(med + criteria.burstMadThresh * mad,
                                           criteria.burstMinRatio * med);
            std::int64_t inBurst = 0;
            for (int c : counts)
                if (c > cutoff) inBurst += c;

            out.burstFraction = static_cast<double>(inBurst) / static_cast<double>(n);
        }
    }

    return out;
}

}  // namespace

QualityClassification classifyUnits(const detect::DetectResult& detected,
                                    const ClusterResult& clusters,
                                    const MeanWaveforms& waveforms,
                                    const UnitMetrics& metrics, const Probe& probe,
                                    const PostProcessSettings& settings,
                                    const QualityCriteria& criteria)
{
    QualityClassification result;
    result.verdicts.resize(clusters.nClusters);

    // Spike times per unit, for the train statistics. Spikes are already in
    // time order, so each unit's list comes out ascending.
    std::vector<std::vector<double>> timesByUnit(clusters.nClusters);
    {
        const double rate = settings.sampleRate > 0 ? settings.sampleRate : 20000.0;
        for (std::size_t i = 0; i < detected.size(); ++i) {
            const std::int32_t c = clusters.spikeClusters[i];
            if (c > 0 && c <= clusters.nClusters)
                timesByUnit[c - 1].push_back(detected.spikeTimes[i] / rate);
        }
    }

    std::vector<TrainStats> trains(clusters.nClusters);
    parallelFor(static_cast<std::size_t>(clusters.nClusters), [&](std::size_t c) {
        trains[c] = measureTrain(timesByUnit[c], criteria);
    });

    const double sampleRate = settings.sampleRate;
    const int nSites = probe.nSites();

    // Measure shapes first: the amplitude outlier test needs the population
    // median, so it has to see every unit before judging any of them.
    for (int c = 0; c < clusters.nClusters; ++c) {
        result.verdicts[c].shape = measureShape(waveforms, c, clusters.clusterSites[c],
                                                sampleRate, criteria.spreadThreshold);
        measureShankSpread(waveforms, c, probe, criteria.shankUniformFraction,
                           result.verdicts[c].shape);
    }

    const auto medianOf = [](std::vector<double> v) {
        v.erase(std::remove_if(v.begin(), v.end(),
                               [](double x) { return !std::isfinite(x) || x <= 0; }),
                v.end());
        if (v.empty()) return 0.0;
        std::sort(v.begin(), v.end());
        return v[v.size() / 2];
    };

    std::vector<double> snrValues, vppValues;
    for (int c = 0; c < clusters.nClusters; ++c) {
        if (c < static_cast<int>(metrics.unitSNR.size()))
            snrValues.push_back(metrics.unitSNR[c]);
        vppValues.push_back(result.verdicts[c].shape.vpp);
    }
    const double medianSNR = medianOf(snrValues);
    const double medianVpp = medianOf(vppValues);

    const double snrOutlier = medianSNR > 0
                                  ? medianSNR * criteria.snrOutlierFactor
                                  : std::numeric_limits<double>::infinity();
    const double vppOutlier = medianVpp > 0
                                  ? medianVpp * criteria.vppOutlierFactor
                                  : std::numeric_limits<double>::infinity();

    for (int c = 0; c < clusters.nClusters; ++c) {
        UnitVerdict& verdict = result.verdicts[c];

        const UnitShape& s = verdict.shape;

        const double snr =
            c < static_cast<int>(metrics.unitSNR.size()) ? metrics.unitSNR[c] : 0.0;
        const double isi = c < static_cast<int>(metrics.unitISIRatio.size())
                               ? metrics.unitISIRatio[c]
                               : std::numeric_limits<double>::quiet_NaN();
        const double isoDist = c < static_cast<int>(metrics.unitIsoDist.size())
                                   ? metrics.unitIsoDist[c]
                                   : std::numeric_limits<double>::quiet_NaN();

        // -- decisive noise tests -----------------------------------------
        bool decisive = false;

        if (s.troughDepth <= 0) {
            verdict.reasons.push_back("no negative deflection");
            decisive = true;
        }
        if (snr > criteria.maxSNR) {
            verdict.reasons.push_back("SNR " + std::to_string(static_cast<int>(snr)) +
                                      " exceeds anything biological");
            decisive = true;
        }
        // A membrane potential cannot swing and return inside one sample.
        if (s.spikiness > criteria.glitchSpikiness &&
            s.halfWidthMs <= criteria.glitchHalfWidthMs) {
            verdict.reasons.push_back("single-sample excursion, not a spike shape");
            decisive = true;
        }
        // Extreme on BOTH amplitude measures relative to this recording's own
        // units. Either alone would catch genuinely large neurons.
        if (snr > snrOutlier && s.vpp > vppOutlier) {
            verdict.reasons.push_back("amplitude far outside the population");
            decisive = true;
        }

        // A waveform on most of the probe at once is a common-mode artifact,
        // whatever its shape on any one site. Counted in sites, not as a probe
        // fraction: on a 512-site probe a 50-site artifact is only 10% and
        // would slip past a fractional test, yet a point source never reaches
        // anywhere near that many.
        if (s.spatialSpread > criteria.maxSpreadSites ||
            (nSites > 0 && s.spatialSpread >
                               criteria.maxSpreadFraction * static_cast<double>(nSites))) {
            verdict.reasons.push_back("present on " +
                                      std::to_string(s.spatialSpread) +
                                      " sites at once (common mode)");
            decisive = true;
        }

        // The same event on most of one shank at comparable amplitude. A
        // neuron reaching two to four neighbouring contacts is ordinary and
        // passes here; a deflection covering the shank is the shank moving as
        // a whole, which no point source can do. This is the test the
        // probe-wide count above cannot make, because a fully covered 64-site
        // shank is only 12.5% of a 512-site probe.
        if (s.shankSites >= criteria.minShankSitesForTest &&
            s.shankCoverage > criteria.maxShankCoverage) {
            verdict.reasons.push_back(
                "covers " + std::to_string(s.shankAbove) + " of " +
                std::to_string(s.shankSites) +
                " sites on its shank at the same amplitude (common mode)");
            decisive = true;
        }

        // -- spike-train tests, which no waveform can reach -----------------
        // A mean waveform is an average: a unit can look perfectly spike-like
        // and still hold a train that no single cell could produce.
        const TrainStats& train = trains[c];

        if (std::isfinite(train.isiRatio) &&
            train.isiRatio > criteria.decisiveISIRatio) {
            verdict.reasons.push_back(
                "refractory period violated by most short intervals");
            decisive = true;
        }

        if (std::isfinite(train.zeroLagRatio) &&
            train.zeroLagRatio > criteria.maxZeroLagRatio) {
            verdict.reasons.push_back(
                "autocorrelogram peaks at zero lag instead of dipping");
            decisive = true;
        }

        if (train.burstFraction > criteria.maxBurstFraction) {
            verdict.reasons.push_back(
                "most spikes arrive in bursts far above the unit's own rate");
            decisive = true;
        }

        // Recorded for context, but NOT grounds on its own. A large,
        // contaminated unit is multi-unit activity - real neurons, poorly
        // separated - and deleting it would throw away spikes.
        if (std::isfinite(isi) && isi > criteria.noiseISIRatio &&
            !(std::isfinite(train.isiRatio) &&
              train.isiRatio > criteria.decisiveISIRatio))
            verdict.reasons.push_back("refractory violations far too frequent");
        if (s.vpp > vppOutlier)
            verdict.reasons.push_back("amplitude far above the population median");
        if (criteria.maxVpp > 0 && s.vpp > criteria.maxVpp)
            verdict.reasons.push_back("amplitude above the configured ceiling");

        if (decisive) {
            verdict.quality = UnitQuality::Noise;
            ++result.nNoise;
            continue;
        }

        // -- single vs multi ------------------------------------------------
        // Everything here has a real spike shape; the question is only how
        // cleanly it is isolated.
        const bool cleanIsi = std::isfinite(isi) && isi <= criteria.singleMaxISIRatio;
        const bool goodSnr = snr >= criteria.singleMinSNR;
        // An unknown isolation distance is not held against a unit.
        const bool goodIso =
            !std::isfinite(isoDist) || isoDist >= criteria.singleMinIsoDist;

        if (cleanIsi && goodSnr && goodIso) {
            verdict.quality = UnitQuality::Single;
            ++result.nSingle;
        } else {
            verdict.quality = UnitQuality::Multi;
            ++result.nMulti;
            if (!goodSnr) verdict.reasons.push_back("SNR below single-unit cut");
            if (!cleanIsi) verdict.reasons.push_back("refractory violations present");
            if (!goodIso) verdict.reasons.push_back("poorly isolated in feature space");
        }
    }

    return result;
}

namespace {

/// Two-means on one cluster's feature vectors.
///
/// Seeded at the two most distant spikes rather than at random, so the result
/// is reproducible and the split lands on the dominant separation - which is
/// what a hierarchical split by hand finds too.
/// Returns the assignment, or an empty vector when the split is not usable.
std::vector<std::uint8_t> twoMeans(const std::vector<const float*>& features,
                                   int nFeatures, int maxIterations,
                                   double& separation)
{
    separation = 0.0;
    const int n = static_cast<int>(features.size());
    if (n < 4) return {};

    // Seed: the spike furthest from the centroid, and the one furthest from it.
    std::vector<double> centroid(nFeatures, 0.0);
    for (const float* f : features)
        for (int k = 0; k < nFeatures; ++k) centroid[k] += f[k];
    for (double& v : centroid) v /= n;

    const auto dist2ToVec = [&](const float* f, const std::vector<double>& v) {
        double d = 0;
        for (int k = 0; k < nFeatures; ++k) {
            const double e = f[k] - v[k];
            d += e * e;
        }
        return d;
    };

    int seedA = 0;
    double best = -1;
    for (int i = 0; i < n; ++i) {
        const double d = dist2ToVec(features[i], centroid);
        if (d > best) { best = d; seedA = i; }
    }

    std::vector<double> a(nFeatures), b(nFeatures);
    for (int k = 0; k < nFeatures; ++k) a[k] = features[seedA][k];

    int seedB = 0;
    best = -1;
    for (int i = 0; i < n; ++i) {
        const double d = dist2ToVec(features[i], a);
        if (d > best) { best = d; seedB = i; }
    }
    for (int k = 0; k < nFeatures; ++k) b[k] = features[seedB][k];

    std::vector<std::uint8_t> label(n, 0);

    for (int iter = 0; iter < maxIterations; ++iter) {
        bool changed = false;
        for (int i = 0; i < n; ++i) {
            const std::uint8_t which =
                dist2ToVec(features[i], a) <= dist2ToVec(features[i], b) ? 0 : 1;
            if (which != label[i]) { label[i] = which; changed = true; }
        }

        std::vector<double> sumA(nFeatures, 0.0), sumB(nFeatures, 0.0);
        int nA = 0, nB = 0;
        for (int i = 0; i < n; ++i) {
            auto& sum = label[i] == 0 ? sumA : sumB;
            (label[i] == 0 ? nA : nB) += 1;
            for (int k = 0; k < nFeatures; ++k) sum[k] += features[i][k];
        }
        if (nA == 0 || nB == 0) return {};

        for (int k = 0; k < nFeatures; ++k) {
            a[k] = sumA[k] / nA;
            b[k] = sumB[k] / nB;
        }
        if (!changed) break;
    }

    // Separation: centre distance over pooled within-group spread. This is what
    // distinguishes two real modes from an arbitrary cut through one blob.
    double centreDist2 = 0;
    for (int k = 0; k < nFeatures; ++k) {
        const double e = a[k] - b[k];
        centreDist2 += e * e;
    }

    double withinA = 0, withinB = 0;
    int nA = 0, nB = 0;
    for (int i = 0; i < n; ++i) {
        if (label[i] == 0) { withinA += dist2ToVec(features[i], a); ++nA; }
        else               { withinB += dist2ToVec(features[i], b); ++nB; }
    }
    if (nA < 2 || nB < 2) return {};

    const double spread =
        std::sqrt(withinA / nA) + std::sqrt(withinB / nB);
    separation = spread > 0 ? std::sqrt(centreDist2) / (0.5 * spread) : 0.0;

    return label;
}

/// Mean raw waveform of an arbitrary spike subset, on `site`'s neighbourhood.
MeanWaveforms subsetMeanWaveform(const detect::DetectResult& detected,
                                 const std::vector<std::int32_t>& spikes, int site,
                                 const Probe& probe, double bitScaling,
                                 double filterNorm, int nSites)
{
    MeanWaveforms w;
    w.nSamples = detected.filtSamples;
    w.nSamplesRaw = detected.rawSamples;
    w.nSites = nSites;
    w.nLocalSites = detected.nWaveformSites;
    w.nClusters = 1;

    w.global.assign(static_cast<std::size_t>(nSites) * w.nSamples, 0.0f);
    const bool haveRaw = !detected.rawWaveforms.empty() && w.nSamplesRaw > 0;
    if (haveRaw)
        w.globalRaw.assign(static_cast<std::size_t>(nSites) * w.nSamplesRaw, 0.0f);

    if (spikes.empty()) return w;

    const int* neighbors = probe.siteNeighbors(site);
    const double inv = 1.0 / spikes.size();

    for (int k = 0; k < w.nLocalSites; ++k) {
        std::vector<double> accFilt(w.nSamples, 0.0);
        std::vector<double> accRaw(haveRaw ? w.nSamplesRaw : 0, 0.0);

        for (std::int32_t spike : spikes) {
            const std::int16_t* f = detected.filtWaveform(spike, k);
            for (int s = 0; s < w.nSamples; ++s) accFilt[s] += f[s];
            if (haveRaw) {
                const std::int16_t* r = detected.rawWaveform(spike, k);
                for (int s = 0; s < w.nSamplesRaw; ++s) accRaw[s] += r[s];
            }
        }

        const auto finish = [&](std::vector<double>& acc, int nSamples, double scale,
                                std::vector<float>& out) {
            double mean = 0;
            for (int s = 0; s < nSamples; ++s) { acc[s] *= inv; mean += acc[s]; }
            mean /= nSamples;
            float* row = out.data() + static_cast<std::size_t>(neighbors[k]) * nSamples;
            for (int s = 0; s < nSamples; ++s)
                row[s] = static_cast<float>((acc[s] - mean) * scale);
        };

        finish(accFilt, w.nSamples, bitScaling / filterNorm, w.global);
        if (haveRaw) finish(accRaw, w.nSamplesRaw, bitScaling, w.globalRaw);
    }

    return w;
}

}  // namespace

SplitReport splitBimodalClusters(const detect::DetectResult& detected,
                                 ClusterResult& clusters,
                                 const MeanWaveforms& waveforms,
                                 const UnitMetrics& metrics, const Probe& probe,
                                 const PostProcessSettings& settings,
                                 const QualityCriteria& criteria,
                                 const SplitCriteria& splitCriteria)
{
    (void)waveforms;
    (void)metrics;

    SplitReport report;

    const int nFeatures = detected.nFeatures;
    if (nFeatures <= 0 || detected.spikeFeatures.empty()) return report;
    if (detected.filtWaveforms.empty()) return report;

    const std::size_t nSpikes = detected.size();
    const int nSites = probe.nSites();
    const bool haveSecondary =
        detected.nPeaks >= 2 && detected.spikeSites2.size() == nSpikes;

    // Group spikes by cluster once.
    std::vector<std::vector<std::int32_t>> byCluster(clusters.nClusters);
    for (std::size_t i = 0; i < nSpikes; ++i) {
        const std::int32_t c = clusters.spikeClusters[i];
        if (c > 0 && c <= clusters.nClusters)
            byCluster[c - 1].push_back(static_cast<std::int32_t>(i));
    }

    // Spikes to move into a new cluster, per original cluster.
    std::vector<std::vector<std::int32_t>> moved(clusters.nClusters);

    parallelFor(static_cast<std::size_t>(clusters.nClusters), [&](std::size_t index) {
        const int c = static_cast<int>(index);
        const auto& spikes = byCluster[c];
        if (static_cast<int>(spikes.size()) < 2 * splitCriteria.minHalfSize) return;

        const int site = clusters.clusterSites[c];

        std::vector<const float*> features;
        features.reserve(spikes.size());
        for (std::int32_t spike : spikes) {
            const int peak =
                (haveSecondary && detected.spikeSites2[spike] == site) ? 1 : 0;
            features.push_back(detected.featuresFor(spike, peak));
        }

        double separation = 0.0;
        const auto label =
            twoMeans(features, nFeatures, splitCriteria.maxIterations, separation);
        if (label.empty()) return;

        std::vector<std::int32_t> groupA, groupB;
        for (std::size_t i = 0; i < spikes.size(); ++i)
            (label[i] == 0 ? groupA : groupB).push_back(spikes[i]);

        if (static_cast<int>(groupA.size()) < splitCriteria.minHalfSize ||
            static_cast<int>(groupB.size()) < splitCriteria.minHalfSize)
            return;

        if (separation < splitCriteria.minSeparation) return;

        // Judge each half on its own mean waveform, exactly as a whole cluster
        // would be judged.
        const auto shapeOf = [&](const std::vector<std::int32_t>& group) {
            const MeanWaveforms w =
                subsetMeanWaveform(detected, group, site, probe, settings.bitScaling,
                                   settings.filterNorm, nSites);
            return measureShape(w, 0, site, settings.sampleRate,
                                criteria.spreadThreshold);
        };

        const UnitShape shapeA = shapeOf(groupA);
        const UnitShape shapeB = shapeOf(groupB);

        const auto looksLikeNoise = [&](const UnitShape& s) {
            if (s.troughDepth <= 0) return true;
            return s.spikiness > criteria.glitchSpikiness &&
                   s.halfWidthMs <= criteria.glitchHalfWidthMs;
        };

        bool noiseA = looksLikeNoise(shapeA);
        bool noiseB = looksLikeNoise(shapeB);

        // When neither half fails on shape, fall back to relative amplitude:
        // a mode many times larger than its sibling is the artifact population.
        if (!noiseA && !noiseB) {
            const double ratio = shapeA.vpp > 0 && shapeB.vpp > 0
                                     ? std::max(shapeA.vpp, shapeB.vpp) /
                                           std::min(shapeA.vpp, shapeB.vpp)
                                     : 1.0;
            if (ratio >= criteria.vppOutlierFactor) {
                if (shapeA.vpp > shapeB.vpp) noiseA = true;
                else noiseB = true;
            }
        }

        if (splitCriteria.requireOneHalfNoise && !(noiseA ^ noiseB)) return;

        // Move the noise half out into its own cluster.
        moved[c] = noiseA ? groupA : groupB;
    });

    // Apply: each moved group becomes a new cluster appended at the end.
    std::int32_t next = clusters.nClusters;
    for (int c = 0; c < clusters.nClusters; ++c) {
        if (moved[c].empty()) continue;
        ++report.nSplit;
        report.splitParents.push_back(c + 1);
        ++next;
        for (std::int32_t spike : moved[c]) {
            clusters.spikeClusters[spike] = next;
            ++report.nSpikesMovedOut;
        }
    }
    report.nClustersTested = clusters.nClusters;

    if (report.nSplit == 0) return report;

    // Rebuild the per-cluster tables for the enlarged set.
    clusters.nClusters = next;
    clusters.unitCount.assign(next, 0);
    clusters.clusterSites.assign(next, 0);
    clusters.clusterCenters.resize(next, 0);

    std::vector<std::vector<std::int32_t>> votes(next);
    for (std::size_t i = 0; i < nSpikes; ++i) {
        const std::int32_t c = clusters.spikeClusters[i];
        if (c <= 0) continue;
        ++clusters.unitCount[c - 1];
        votes[c - 1].push_back(detected.spikeSites[i]);
    }
    for (std::int32_t c = 0; c < next; ++c) {
        if (votes[c].empty()) continue;
        std::sort(votes[c].begin(), votes[c].end());
        std::int32_t best = votes[c][0], bestRun = 0, current = votes[c][0], run = 0;
        for (std::int32_t s : votes[c]) {
            run = (s == current) ? run + 1 : 1;
            current = s;
            if (run > bestRun) { bestRun = run; best = current; }
        }
        clusters.clusterSites[c] = best;
    }

    return report;
}

std::vector<std::string> annotationStrings(const QualityClassification& classification,
                                           const QualityCriteria& criteria)
{
    std::vector<std::string> notes(classification.verdicts.size());

    for (std::size_t c = 0; c < classification.verdicts.size(); ++c) {
        switch (classification.verdicts[c].quality) {
            case UnitQuality::Noise:
                notes[c] = criteria.noiseLabel;
                break;
            case UnitQuality::Single:
                if (criteria.labelGoodUnits) notes[c] = "single";
                break;
            case UnitQuality::Multi:
                if (criteria.labelGoodUnits) notes[c] = "multi";
                break;
        }
    }

    return notes;
}

}  // namespace sort
}  // namespace jrc
