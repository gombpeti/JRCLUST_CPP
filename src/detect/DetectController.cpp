#include "detect/DetectController.h"

#include "core/Parallel.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <future>
#include <stdexcept>

namespace jrc {
namespace detect {

DetectSettings DetectSettings::fromConfig(const Config& config)
{
    DetectSettings settings;

    settings.filterType = filterTypeFromString(config.getString("filterType", "ndiff"));
    settings.nDiffOrder = config.getInt("nDiffOrder", 2);
    settings.carMode = carModeFromString(config.getString("CARMode", "mean"));

    settings.peaks.qqFactor = config.getNumber("qqFactor", 4.5);
    settings.peaks.minNeighborsDetect = config.getInt("minNeighborsDetect", 0);
    settings.peaks.detectBipolar = config.getInt("detectBipolar", 0) != 0;
    settings.peaks.evtDetectRad = config.getNumber("evtDetectRad", 50.0);

    // refracInt is in milliseconds; detection works in samples.
    const double sampleRate = config.sampleRate();
    const double refracMs = config.getNumber("refracInt", 0.25);
    settings.peaks.refracIntSamples =
        static_cast<int>(std::round(std::fabs(refracMs) * sampleRate / 1000.0));

    settings.peaks.removeCoincident =
        config.getInt("removeCoincident", 1) != 0;
    settings.peaks.coincidenceSamples = static_cast<int>(std::round(
        std::fabs(config.getNumber("coincidenceMs", 0.3)) * sampleRate / 1000.0));
    settings.peaks.coincidenceAmpRatio =
        config.getNumber("coincidenceAmpRatio", 0.7);

    // spikeThreshMax is in µV; compare against raw ADC counts.
    const auto threshMax = config.getNumbers("spikeThreshMax");
    if (!threshMax.empty() && config.bitScaling() != 0.0)
        settings.peaks.spikeThreshMaxRaw =
            std::fabs(threshMax.front()) / config.bitScaling();

    settings.nSamplesPad = config.getInt("nSamplesPad", 100);
    settings.featureSettings = features::FeatureSettings::fromConfig(config);

    // evtWindowRaw, in milliseconds, gives the wider raw waveform window.
    const auto rawWindow = config.getNumbers("evtWindowRaw");
    if (rawWindow.size() >= 2) {
        settings.rawPreSamples = static_cast<int>(
            std::round(std::fabs(rawWindow[0]) * sampleRate / 1000.0));
        settings.rawPostSamples =
            static_cast<int>(std::round(rawWindow[1] * sampleRate / 1000.0));
    }

    // Windows must fit inside the chunk overlap, or a spike near a boundary
    // would have its waveform clamped instead of read from the neighbouring
    // samples. Widen the pad rather than silently truncating waveforms.
    // The raw window is the wider of the two, so it sets the requirement.
    const int needed = std::max({settings.featureSettings.preSamples,
                                 settings.featureSettings.postSamples,
                                 settings.rawPreSamples, settings.rawPostSamples}) + 1;
    settings.nSamplesPad = std::max(settings.nSamplesPad, needed);

    return settings;
}

DetectController::DetectController(const Config& config, const Probe& probe,
                                   const DetectSettings& settings)
    : config_(config), probe_(probe), settings_(settings)
{
}

void DetectController::planChunks(std::int64_t nSamples, int& nChunks,
                                  std::int64_t& chunkSamples) const
{
    chunkSamples = static_cast<std::int64_t>(settings_.chunkSeconds *
                                             config_.sampleRate());
    chunkSamples = std::max<std::int64_t>(chunkSamples, 4 * settings_.nSamplesPad + 16);
    chunkSamples = std::min<std::int64_t>(chunkSamples, nSamples);

    nChunks = static_cast<int>((nSamples + chunkSamples - 1) / chunkSamples);
    if (nChunks < 1) nChunks = 1;
}

void DetectController::loadAndFilter(Recording& recording, std::int64_t begin,
                                     std::int64_t end, std::vector<float>& out,
                                     int& nSamplesOut)
{
    loadChunk(recording, begin, end, out, nSamplesOut);
    filterChunk(out, nSamplesOut);
}

void DetectController::loadChunk(Recording& recording, std::int64_t begin,
                                 std::int64_t end, std::vector<float>& out,
                                 int& nSamplesOut)
{
    nSamplesOut = static_cast<int>(end - begin);

    // The file is channel-ordered and sample-major; detection wants
    // site-ordered and site-major. siteMap gives the channel each site is
    // recorded on (readRawROI(siteMap, ...) in MATLAB). Read, convert and
    // transpose happen together so the samples are touched exactly once.
    const auto tRead = std::chrono::steady_clock::now();
    recording.readSitesTransposed(begin, end, probe_.siteMap(), out);
    timings_.read +=
        std::chrono::duration<double>(std::chrono::steady_clock::now() - tRead).count();
}

void DetectController::filterChunk(std::vector<float>& out, int nSamplesOut)
{
    const int nSites = probe_.nSites();

    const auto tFilter = std::chrono::steady_clock::now();

    switch (settings_.filterType) {
        case FilterType::NDiff:
            ndiffFilter(out, nSamplesOut, nSites, settings_.nDiffOrder);
            break;
        case FilterType::None:
            break;
        default:
            throw std::runtime_error(
                std::string("filterType '") + filterTypeName(settings_.filterType) +
                "' is not implemented yet (only 'ndiff' and 'none' so far)");
    }

    const auto tCar = std::chrono::steady_clock::now();
    timings_.filter += std::chrono::duration<double>(tCar - tFilter).count();

    std::vector<float> channelMeans;
    applyCar(out, nSamplesOut, nSites, settings_.carMode, probe_.ignoredMask(),
             channelMeans);

    timings_.car +=
        std::chrono::duration<double>(std::chrono::steady_clock::now() - tCar).count();
}

std::vector<std::int64_t> DetectController::samplingWindows(
    std::int64_t begin, std::int64_t end, std::int64_t& windowSamples) const
{
    const std::int64_t nSamples = end - begin;

    // Windows are defined in absolute time, not in chunks, so the thresholds
    // (and therefore the spikes) do not move when chunkSeconds changes.
    windowSamples = static_cast<std::int64_t>(settings_.thresholdWindowSeconds *
                                              config_.sampleRate());
    windowSamples = std::max<std::int64_t>(windowSamples, 16);
    windowSamples = std::min(windowSamples, nSamples);

    const int maxWindows =
        static_cast<int>(std::max<std::int64_t>(1, nSamples / windowSamples));
    const int nWindows = (settings_.thresholdWindows <= 0)
                             ? maxWindows
                             : std::min(settings_.thresholdWindows, maxWindows);

    std::vector<std::int64_t> starts;
    starts.reserve(nWindows);
    for (int i = 0; i < nWindows; ++i) {
        // Evenly spaced so drift and varying noise across the session are
        // represented, rather than sampling only the opening minutes.
        const std::int64_t span = nSamples - windowSamples;
        const std::int64_t offset =
            (nWindows == 1) ? span / 2
                            : (span * static_cast<std::int64_t>(i)) / (nWindows - 1);
        starts.push_back(begin + offset);
    }

    return starts;
}

features::PrincipalVectors DetectController::estimatePrincipalVectors(
    Recording& recording, std::int64_t begin, std::int64_t end,
    const std::vector<float>& thresholds)
{
    features::SpikeWindows pooled;
    pooled.nSamples = settings_.featureSettings.windowSamples();
    pooled.nSites = 1;  // principal vectors come from the primary site only
    pooled.nSpikes = 0;

    collectPrincipalVectorSamples(recording, begin, end, thresholds, pooled);
    return features::computePrincipalVectors(pooled, settings_.featureSettings);
}

void DetectController::collectPrincipalVectorSamples(
    Recording& recording, std::int64_t begin, std::int64_t end,
    const std::vector<float>& thresholds, features::SpikeWindows& pooled)
{
    const int nSites = probe_.nSites();

    std::int64_t windowSamples = 0;
    const std::vector<std::int64_t> windowStarts =
        samplingWindows(begin, end, windowSamples);

    const int nUsableWindows = std::max<int>(1, static_cast<int>(windowStarts.size()));
    const int perWindowCap =
        std::max(1, settings_.featureSettings.maxPcaSpikes / nUsableWindows);

    std::vector<float> filtered;

    for (std::int64_t windowBegin : windowStarts) {
        const std::int64_t windowEnd = std::min(windowBegin + windowSamples, end);
        if (windowEnd - windowBegin < 16) continue;

        int nSamplesOut = 0;
        loadAndFilter(recording, windowBegin, windowEnd, filtered, nSamplesOut);

        const SpikeTable spikes = detectSpikes(filtered, nSamplesOut, nSites,
                                               thresholds, probe_, settings_.peaks);
        if (spikes.size() == 0) continue;

        // Take an even spread within the window rather than the first N.
        const int stride =
            std::max<int>(1, static_cast<int>(spikes.size()) / perWindowCap);

        std::vector<std::int64_t> times;
        std::vector<std::int32_t> sites;
        for (std::size_t i = 0; i < spikes.size(); i += stride) {
            times.push_back(spikes.times[i]);
            sites.push_back(spikes.sites[i]);
        }
        if (times.empty()) continue;

        features::SpikeWindows windows =
            features::extractWindows(filtered, nSamplesOut, nSites, times, sites,
                                     probe_, settings_.featureSettings);
        features::applyLocalCar(windows, settings_.featureSettings);

        // Keep only site 0 - the primary trace - which is what spikePrVecs.m
        // builds the covariance from.
        const int nWin = windows.nSamples;
        const std::size_t offset = pooled.data.size();
        pooled.data.resize(offset + static_cast<std::size_t>(windows.nSpikes) * nWin);
        std::memcpy(pooled.data.data() + offset, windows.at(0, 0),
                    static_cast<std::size_t>(windows.nSpikes) * nWin * sizeof(float));
        pooled.nSpikes += windows.nSpikes;
    }
}

std::vector<float> DetectController::estimateThresholds(Recording& recording,
                                                        std::int64_t begin,
                                                        std::int64_t end)
{
    const int nSites = probe_.nSites();

    std::int64_t windowSamples = 0;
    const std::vector<std::int64_t> windowStarts =
        samplingWindows(begin, end, windowSamples);

    std::vector<double> accumulated(nSites, 0.0);
    int nUsed = 0;

    std::vector<float> filtered;
    for (std::int64_t chunkBegin : windowStarts) {
        const std::int64_t chunkEnd = std::min(chunkBegin + windowSamples, end);
        if (chunkEnd - chunkBegin < 16) continue;

        int nSamplesOut = 0;
        loadAndFilter(recording, chunkBegin, chunkEnd, filtered, nSamplesOut);

        const auto tThresh = std::chrono::steady_clock::now();
        const auto thresholds =
            computeThresholds(filtered, nSamplesOut, nSites, settings_.peaks.qqFactor,
                              probe_.ignoredMask());
        timings_.threshold += std::chrono::duration<double>(
                                  std::chrono::steady_clock::now() - tThresh)
                                  .count();

        for (int site = 0; site < nSites; ++site) accumulated[site] += thresholds[site];
        ++nUsed;
    }

    std::vector<float> result(nSites, 0.0f);
    if (nUsed > 0) {
        // detectOneRecording.m averages the per-chunk thresholds the same way.
        for (int site = 0; site < nSites; ++site)
            result[site] = static_cast<float>(accumulated[site] / nUsed);
    }

    return result;
}

void appendDetectResult(DetectResult& combined, const DetectResult& part,
                        std::int64_t sampleOffset)
{
    if (part.size() == 0) return;

    const bool first = combined.size() == 0 && combined.filtSamples == 0;

    if (first) {
        combined.nFeatures = part.nFeatures;
        combined.nPeaks = part.nPeaks;
        combined.filtSamples = part.filtSamples;
        combined.rawSamples = part.rawSamples;
        combined.nWaveformSites = part.nWaveformSites;
        combined.principalVectors = part.principalVectors;
        combined.siteThresholds = part.siteThresholds;
    } else {
        // Features from different bases cannot be clustered together.
        if (part.nFeatures != combined.nFeatures || part.nPeaks != combined.nPeaks)
            throw std::runtime_error(
                "cannot concatenate recordings with different feature layouts");
        if (part.filtSamples != combined.filtSamples ||
            part.nWaveformSites != combined.nWaveformSites)
            throw std::runtime_error(
                "cannot concatenate recordings with different waveform shapes");
    }

    for (std::int64_t t : part.spikeTimes)
        combined.spikeTimes.push_back(t + sampleOffset);

    combined.spikeAmps.insert(combined.spikeAmps.end(), part.spikeAmps.begin(),
                              part.spikeAmps.end());
    combined.spikeSites.insert(combined.spikeSites.end(), part.spikeSites.begin(),
                               part.spikeSites.end());
    combined.spikeSites2.insert(combined.spikeSites2.end(), part.spikeSites2.begin(),
                                part.spikeSites2.end());
    combined.spikeFeatures.insert(combined.spikeFeatures.end(),
                                  part.spikeFeatures.begin(),
                                  part.spikeFeatures.end());
    combined.spikePositions.insert(combined.spikePositions.end(),
                                   part.spikePositions.begin(),
                                   part.spikePositions.end());
    combined.filtWaveforms.insert(combined.filtWaveforms.end(),
                                  part.filtWaveforms.begin(),
                                  part.filtWaveforms.end());
    combined.rawWaveforms.insert(combined.rawWaveforms.end(),
                                 part.rawWaveforms.begin(), part.rawWaveforms.end());

    combined.veto.add(part.veto);

    combined.timings.read += part.timings.read;
    combined.timings.filter += part.timings.filter;
    combined.timings.car += part.timings.car;
    combined.timings.threshold += part.timings.threshold;
    combined.timings.peaks += part.timings.peaks;
    combined.timings.features += part.timings.features;
    combined.timings.total += part.timings.total;
}

DetectResult DetectController::run(Recording& recording, const ProgressFn& progress)
{
    return runRange(recording, 0, recording.nSamples(), progress);
}

DetectResult DetectController::runRange(Recording& recording, std::int64_t begin,
                                        std::int64_t end, const ProgressFn& progress)
{
    const int nSites = probe_.nSites();

    begin = std::max<std::int64_t>(0, begin);
    end = std::min(end, recording.nSamples());
    if (end <= begin) throw std::runtime_error("empty sample range for detection");

    const std::int64_t nSamples = end - begin;

    timings_ = DetectTimings{};
    const auto tStart = std::chrono::steady_clock::now();

    DetectResult result;
    result.siteThresholds = fixedThresholds_ ? *fixedThresholds_
                                             : estimateThresholds(recording, begin, end);

    const auto& fs = settings_.featureSettings;
    const int nPeaks = settings_.computeFeatures ? std::max(1, fs.nPeaksFeatures) : 0;

    if (settings_.computeFeatures) {
        const auto tPca = std::chrono::steady_clock::now();
        result.principalVectors =
            fixedPrincipalVectors_
                ? *fixedPrincipalVectors_
                : estimatePrincipalVectors(recording, begin, end, result.siteThresholds);
        timings_.features += std::chrono::duration<double>(
                                 std::chrono::steady_clock::now() - tPca)
                                 .count();

        result.nFeatures = fs.nFeatures();
        result.nPeaks = nPeaks;
    }

    int nChunks = 0;
    std::int64_t chunkSamples = 0;
    planChunks(nSamples, nChunks, chunkSamples);

    const std::int64_t pad = settings_.nSamplesPad;

    // Reading dominates on a spinning or USB-attached drive, so the next chunk
    // is fetched on a background thread while the current one is filtered and
    // searched. Two buffers ping-pong between the reader and the processor.
    std::vector<float> buffers[2];
    int bufferSamples[2] = {0, 0};

    const auto chunkBounds = [&](int chunk, std::int64_t& loadBegin,
                                 std::int64_t& loadEnd, std::int64_t& coreBegin,
                                 std::int64_t& coreEnd) {
        coreBegin = begin + static_cast<std::int64_t>(chunk) * chunkSamples;
        coreEnd = std::min(coreBegin + chunkSamples, end);
        // Overlap neighbouring chunks so a spike sitting on a boundary is seen
        // with full context; detections inside the pad are discarded because
        // the adjacent chunk owns them.
        loadBegin = std::max<std::int64_t>(begin, coreBegin - pad);
        loadEnd = std::min(end, coreEnd + pad);
    };

    const auto fetch = [&](int chunk, int slot) {
        std::int64_t lb = 0, le = 0, cb = 0, ce = 0;
        chunkBounds(chunk, lb, le, cb, ce);
        if (ce <= cb) {
            bufferSamples[slot] = 0;
            return;
        }
        loadChunk(recording, lb, le, buffers[slot], bufferSamples[slot]);
    };

    fetch(0, 0);
    std::future<void> pending;

    for (int chunk = 0; chunk < nChunks; ++chunk) {
        const int slot = chunk % 2;

        // Kick off the next read before touching this chunk's data.
        if (chunk + 1 < nChunks) {
            pending = std::async(std::launch::async, fetch, chunk + 1, (chunk + 1) % 2);
        }

        std::int64_t loadBegin = 0, loadEnd = 0, coreBegin = 0, coreEnd = 0;
        chunkBounds(chunk, loadBegin, loadEnd, coreBegin, coreEnd);

        int nSamplesOut = bufferSamples[slot];
        if (coreEnd <= coreBegin || nSamplesOut <= 0) {
            if (pending.valid()) pending.get();
            continue;
        }

        std::vector<float>& filtered = buffers[slot];

        // The curation GUI shows unfiltered traces, so keep a copy before the
        // filter overwrites the chunk in place. Only paid for when saving.
        std::vector<float> unfiltered;
        if ((writer_ || settings_.retainWaveforms) && settings_.computeFeatures)
            unfiltered = filtered;

        // Zero-filled recording gaps, found on the RAW chunk before filtering
        // turns their edges into non-zero values. Each edge is a step on every
        // channel at once, which detection would otherwise read as an enormous
        // spike on every site.
        std::vector<std::uint8_t> gapMask;
        if (settings_.veto.rejectGaps) {
            const int pad = settings_.veto.gapPadSamples > 0
                                ? settings_.veto.gapPadSamples
                                : std::max({settings_.rawPreSamples,
                                            settings_.rawPostSamples,
                                            settings_.featureSettings.preSamples,
                                            settings_.featureSettings.postSamples}) + 2;

            const std::int64_t blanked = findRecordingGaps(
                filtered, nSamplesOut, nSites, settings_.veto.gapZeroFraction, pad,
                gapMask);
            result.veto.gapSamples += blanked;
            if (blanked == 0) gapMask.clear();
        }

        filterChunk(filtered, nSamplesOut);

        const auto tPeaks = std::chrono::steady_clock::now();
        const SpikeTable spikes =
            detectSpikes(filtered, nSamplesOut, nSites, result.siteThresholds, probe_,
                         settings_.peaks);
        timings_.peaks +=
            std::chrono::duration<double>(std::chrono::steady_clock::now() - tPeaks)
                .count();

        // Keep only spikes whose time falls in this chunk's own span.
        const std::int64_t coreOffset = coreBegin - loadBegin;
        const std::int64_t coreLength = coreEnd - coreBegin;

        // Gather this chunk's spikes but do NOT commit them yet: the per-spike
        // shape veto below needs their waveforms, and anything it rejects must
        // never reach the result at all.
        std::vector<std::int64_t> keptTimes;
        std::vector<std::int32_t> keptSites;
        std::vector<float> keptAmps;

        for (std::size_t i = 0; i < spikes.size(); ++i) {
            const std::int64_t local = spikes.times[i];
            if (local < coreOffset || local >= coreOffset + coreLength) continue;

            // Anything inside a recording gap is a step edge, not a spike.
            if (!gapMask.empty() && local < static_cast<std::int64_t>(gapMask.size()) &&
                gapMask[static_cast<std::size_t>(local)]) {
                ++result.veto.inGap;
                continue;
            }

            keptTimes.push_back(local);  // chunk-local, for window extraction
            keptSites.push_back(spikes.sites[i]);
            keptAmps.push_back(spikes.amplitudes[i]);
        }

        if (settings_.computeFeatures && !keptTimes.empty()) {
            const auto tFeat = std::chrono::steady_clock::now();

            // Peak 1: windows centred on each spike's own site.
            features::SpikeWindows windows =
                features::extractWindows(filtered, nSamplesOut, nSites, keptTimes,
                                         keptSites, probe_, fs);

            // Per-spike shape veto, while the full neighbourhood is still
            // present (the decay test needs the far sites). Rejected events
            // are dropped here so nothing downstream can adopt them.
            //
            // It runs on the RAW waveform, not the filtered one. The shape
            // tests are physiological - how long a membrane takes to move, how
            // a point source falls off with distance - and ndiff is a
            // differentiator that sharpens every waveform, so judging
            // "too narrow" on filtered data rejects real spikes wholesale.
            if (settings_.veto.enabled && !unfiltered.empty()) {
                features::FeatureSettings vetoFs = fs;
                vetoFs.preSamples = settings_.rawPreSamples;
                vetoFs.postSamples = settings_.rawPostSamples;

                const features::SpikeWindows vetoWindows =
                    features::extractWindows(unfiltered, nSamplesOut, nSites,
                                             keptTimes, keptSites, probe_, vetoFs);

                std::vector<std::uint8_t> keep;
                const VetoCounts counts =
                    vetoSpikes(vetoWindows, windows, settings_.veto,
                               config_.sampleRate(), keep);
                result.veto.add(counts);

                if (counts.rejected > 0) {
                    compactWindows(windows, keep);

                    std::size_t out = 0;
                    for (std::size_t i = 0; i < keptTimes.size(); ++i) {
                        if (!keep[i]) continue;
                        keptTimes[out] = keptTimes[i];
                        keptSites[out] = keptSites[i];
                        keptAmps[out] = keptAmps[i];
                        ++out;
                    }
                    keptTimes.resize(out);
                    keptSites.resize(out);
                    keptAmps.resize(out);
                }

                if (keptTimes.empty()) {
                    if (progress) progress(chunk + 1, nChunks, result.size());
                    if (pending.valid()) pending.get();
                    continue;
                }
            }

            // The secondary site is chosen from the raw neighbourhood, before
            // local CAR trims it, exactly as findSecondaryPeaks.m does.
            std::vector<std::int32_t> secondary;
            if (nPeaks >= 2)
                secondary = features::findSecondarySites(windows, keptSites, probe_, fs);

            // JRCLUST stores the untrimmed neighbourhood waveforms, so capture
            // them before local CAR reduces the site count.
            const bool needWaveforms = writer_ || settings_.retainWaveforms;
            if (needWaveforms) {
                if (writer_) writer_->appendFilt(windows);
                if (settings_.retainWaveforms) {
                    io::packWindowsInt16(windows, result.filtWaveforms);
                    result.filtSamples = windows.nSamples;
                    result.nWaveformSites = windows.nSites;
                }

                features::FeatureSettings rawFs = fs;
                rawFs.preSamples = settings_.rawPreSamples;
                rawFs.postSamples = settings_.rawPostSamples;

                const features::SpikeWindows rawWindows =
                    features::extractWindows(unfiltered, nSamplesOut, nSites, keptTimes,
                                             keptSites, probe_, rawFs);

                if (writer_) writer_->appendRaw(rawWindows);
                if (settings_.retainWaveforms) {
                    io::packWindowsInt16(rawWindows, result.rawWaveforms);
                    result.rawSamples = rawWindows.nSamples;
                }
            }

            features::applyLocalCar(windows, fs);
            const features::SpikeFeatures f1 =
                features::projectFeatures(windows, result.principalVectors, fs);

            features::SpikeFeatures f2;
            if (nPeaks >= 2) {
                features::SpikeWindows windows2 =
                    features::extractWindows(filtered, nSamplesOut, nSites, keptTimes,
                                             secondary, probe_, fs);
                features::applyLocalCar(windows2, fs);
                f2 = features::projectFeatures(windows2, result.principalVectors, fs);

                result.spikeSites2.insert(result.spikeSites2.end(), secondary.begin(),
                                          secondary.end());
            }

            // Commit the survivors. Features are appended in the same order
            // just below, so the two stay aligned.
            for (std::size_t i = 0; i < keptTimes.size(); ++i) {
                result.spikeTimes.push_back(loadBegin + keptTimes[i]);
                result.spikeAmps.push_back(keptAmps[i]);
                result.spikeSites.push_back(keptSites[i]);
            }

            const int nFeat = result.nFeatures;
            const int nKept = static_cast<int>(keptTimes.size());
            const std::size_t before = result.spikeFeatures.size();

            for (int i = 0; i < nKept; ++i) {
                result.spikeFeatures.insert(result.spikeFeatures.end(), f1.at(i),
                                            f1.at(i) + nFeat);
                if (nPeaks >= 2)
                    result.spikeFeatures.insert(result.spikeFeatures.end(), f2.at(i),
                                                f2.at(i) + nFeat);
            }

            // [spike][peak][feature] with feature contiguous is byte-identical
            // to MATLAB's column-major [nFeatures, nPeaks, nSpikes], so this
            // streams straight to disk with no rearrangement.
            if (writer_)
                writer_->appendFeatures(result.spikeFeatures.data() + before,
                                        result.spikeFeatures.size() - before);

            timings_.features += std::chrono::duration<double>(
                                     std::chrono::steady_clock::now() - tFeat)
                                     .count();
        } else {
            // No features requested, so no waveforms were cut and the veto had
            // nothing to inspect. Commit the spikes as detected.
            for (std::size_t i = 0; i < keptTimes.size(); ++i) {
                result.spikeTimes.push_back(loadBegin + keptTimes[i]);
                result.spikeAmps.push_back(keptAmps[i]);
                result.spikeSites.push_back(keptSites[i]);
            }
        }

        if (progress) progress(chunk + 1, nChunks, result.size());

        // Wait for the prefetch so its buffer is complete before the next pass.
        if (pending.valid()) pending.get();
    }

    // Feature-weighted spike positions, for drift-aware mean waveforms.
    // Port of spikePos.m: weight each neighbourhood site by its squared PC1
    // amplitude and take the centroid. The first nSitesEvt features are PC1
    // across the sites, which is the layout projectFeatures produces.
    if (settings_.computeFeatures && result.nFeatures > 0 && !result.spikeFeatures.empty()) {
        const int nSitesEvt = std::min(fs.nSitesEvt(), result.nFeatures);
        const std::size_t nSpikes = result.size();
        result.spikePositions.assign(nSpikes * 2, 0.0f);

        parallelFor(nSpikes, [&](std::size_t i) {
            const float* f = result.featuresFor(i, 0);
            const int* neighbors = probe_.siteNeighbors(result.spikeSites[i]);

            double weight = 0.0, x = 0.0, y = 0.0;
            for (int k = 0; k < nSitesEvt; ++k) {
                const double w = static_cast<double>(f[k]) * f[k];
                weight += w;
                x += w * probe_.x(neighbors[k]);
                y += w * probe_.y(neighbors[k]);
            }

            if (weight > 0) {
                result.spikePositions[i * 2] = static_cast<float>(x / weight);
                result.spikePositions[i * 2 + 1] = static_cast<float>(y / weight);
            } else {
                result.spikePositions[i * 2] =
                    static_cast<float>(probe_.x(result.spikeSites[i]));
                result.spikePositions[i * 2 + 1] =
                    static_cast<float>(probe_.y(result.spikeSites[i]));
            }
        });
    }

    timings_.total =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - tStart).count();
    result.timings = timings_;

    return result;
}

}  // namespace detect
}  // namespace jrc
