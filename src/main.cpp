// jrclust - native C++/CUDA port of the JRCLUST spike sorter.
//
// Stage 1 scope: configuration parsing, raw recording access, and the CUDA
// density-peak kernels carried over unmodified from the MATLAB implementation.

#include "core/Config.h"
#include "core/Parallel.h"
#include "core/Probe.h"
#include "core/Recording.h"
#include "cuda/CudaUtils.h"
#include "cuda/DensityPeak.h"
#include "detect/DetectController.h"
#include "sort/Clustering.h"
#include "sort/PostProcess.h"
#include "sort/UnitClassifier.h"

#include <algorithm>
#include <functional>
#include <chrono>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <map>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <vector>

namespace {

void printUsage()
{
    std::cout
        << "jrclust - native C++/CUDA spike sorter\n\n"
        << "usage:\n"
        << "  jrclust gpuinfo               show CUDA devices\n"
        << "  jrclust info <file.prm>       parse config, open recording, report\n"
        << "  jrclust probe <file.prm>      show probe geometry and site neighbours\n"
        << "  jrclust detect <file.prm>     filter, threshold and detect spikes\n"
        << "  jrclust sort <file.prm>       detect, then cluster (rho/delta on GPU)\n"
        << "  jrclust selftest [nSpikes]    verify CUDA kernels against CPU reference\n"
        << "\n"
        << "detect options:\n"
        << "  --seconds <n>        process only the first n seconds\n"
        << "  --chunk <n>          chunk size in seconds (default 8)\n"
        << "  --recording <n>      which entry of rawRecordings (default 0)\n"
        << "\n"
        << "sort options:\n"
        << "  --seconds <n>        process only the first n seconds\n"
        << "  --recording <n>      sort ONLY entry n of rawRecordings.\n"
        << "                       Default is to concatenate all of them into\n"
        << "                       one sort, as JRCLUST does.\n"
        << "  --recording-path <p> override the recording path in the .prm\n"
        << "  --output-dir <d>     override outputDir from the .prm\n"
        << "  --save               write .jrc files for MATLAB curation\n"
        << "  --no-merge           skip auto-merge\n"
        << "  --rm-refrac          drop refractory-period violations\n"
        << "  --rm-outliers        drop spikes far from their own unit in\n"
        << "                       feature space (outlierThresh in the .prm)\n"
        << "  --no-split           do not split clusters holding two\n"
        << "                       populations (splitting is on by default)\n"
        << "  --veto-spikes        reject individual detected events whose\n"
        << "                       waveform cannot be an action potential:\n"
        << "                       too narrow, single-sample, or no spatial\n"
        << "                       decay across sites\n"
        << "  --veto-decay <r>     far/near amplitude ratio above which an\n"
        << "                       event counts as flat (default 0.65)\n"
        << "  --refrac-merge       veto merges whose spike trains fire\n"
        << "                       independently (better than waveform alone)\n"
        << "  --matlab-lratio      reproduce MATLAB's L-ratio (uses spike count\n"
        << "                       as chi-squared dof; not the correct statistic)\n"
        << "  --cpu                run rho/delta on the CPU reference\n"
        << "  --no-classify        skip noise/single/multi labelling\n"
        << "  --mark-for-delete    label noise 'to_delete' so JRCLUST's\n"
        << "                       Edit > Delete annotated clears them\n"
        << "                       (--arm-delete does the same)\n"
        << "  --label-noise-only   annotate only suspect units, leaving good\n"
        << "                       ones blank, so 'Display annotated units\n"
        << "                       only' shows just the suspects. Default is to\n"
        << "                       label every unit noise/single/multi\n"
        << "  --max-vpp <uV>       treat units above this amplitude as suspect\n"
        << "\n";
}

std::string humanBytes(std::size_t bytes)
{
    const char* units[] = {"B", "KiB", "MiB", "GiB", "TiB"};
    double value = static_cast<double>(bytes);
    int unit = 0;
    while (value >= 1024.0 && unit < 4) {
        value /= 1024.0;
        ++unit;
    }
    char buffer[64];
    std::snprintf(buffer, sizeof(buffer), "%.2f %s", value, units[unit]);
    return buffer;
}

// ------------------------------------------------------------------ gpuinfo --

int commandGpuInfo()
{
    const auto devices = jrc::cuda::enumerateDevices();

    if (devices.empty()) {
        std::cout << "No CUDA devices found.\n";
        return 1;
    }

    for (const auto& device : devices) {
        std::cout << "device " << device.index << ": " << device.name << "\n"
                  << "  compute capability  sm_" << device.major << device.minor << "\n"
                  << "  global memory       " << humanBytes(device.totalMemory) << "\n"
                  << "  multiprocessors     " << device.multiProcessorCount << "\n"
                  << "  max threads/block   " << device.maxThreadsPerBlock << "\n"
                  << "  shared mem/block    " << humanBytes(device.sharedMemPerBlock)
                  << "\n";
    }

    return 0;
}

// --------------------------------------------------------------------- info --

int commandInfo(const std::string& prmPath)
{
    jrc::Config config;
    config.load(prmPath);

    std::cout << "parameter file: " << config.path() << "\n"
              << "  parameters parsed   " << config.all().size() << "\n"
              << "  nChans              " << config.nChans() << "\n"
              << "  sampleRate          " << config.sampleRate() << " Hz\n"
              << "  dataType            " << jrc::dataTypeName(config.dataType()) << "\n"
              << "  headerOffset        " << config.headerOffset() << " bytes\n"
              << "  bitScaling          " << config.bitScaling() << " uV/bit\n"
              << "  rawRecordings       " << config.rawRecordings().size() << "\n\n";

    if (config.rawRecordings().empty()) {
        std::cout << "No raw recordings listed in the parameter file.\n";
        return 0;
    }

    for (const auto& path : config.rawRecordings()) {
        std::cout << "recording: " << path << "\n";

        try {
            jrc::Recording recording;
            recording.open(path, config);

            std::cout << "  file size           " << humanBytes(recording.fileBytes())
                      << "\n"
                      << "  channels            " << recording.nChans() << "\n"
                      << "  samples             " << recording.nSamples() << "\n"
                      << "  duration            " << std::fixed << std::setprecision(2)
                      << recording.duration() << " s\n";

            // Prove the mapping works by pulling the first frame back.
            std::vector<float> frame;
            const std::int64_t nPreview = std::min<std::int64_t>(1, recording.nSamples());
            recording.readSamplesMicrovolts(0, nPreview, frame);

            std::cout << "  first sample (uV)   ";
            const int shown = std::min<int>(8, recording.nChans());
            for (int c = 0; c < shown; ++c)
                std::cout << std::setprecision(1) << frame[c] << " ";
            if (shown < recording.nChans()) std::cout << "...";
            std::cout << "\n";
        } catch (const std::exception& error) {
            std::cout << "  ERROR: " << error.what() << "\n";
        }

        std::cout << "\n";
    }

    return 0;
}

// -------------------------------------------------------------------- probe --

int commandProbe(const std::string& prmPath)
{
    jrc::Config config;
    config.load(prmPath);

    jrc::Probe probe;
    probe.load(config);

    const double detectRad = config.getNumber("evtDetectRad", 50.0);

    std::cout << "probe from: " << config.path() << "\n"
              << "  sites               " << probe.nSites() << "\n"
              << "  ignored sites       " << probe.nIgnored() << "\n"
              << "  active sites        " << probe.nSites() - probe.nIgnored() << "\n"
              << "  evtDetectRad        " << detectRad << " um\n";

    // Shank census.
    std::vector<int> shanks;
    for (int s = 0; s < probe.nSites(); ++s) shanks.push_back(probe.shank(s));
    std::sort(shanks.begin(), shanks.end());
    shanks.erase(std::unique(shanks.begin(), shanks.end()), shanks.end());
    std::cout << "  shanks              " << shanks.size() << "\n";

    // Extent.
    double xMin = 1e30, xMax = -1e30, yMin = 1e30, yMax = -1e30;
    for (int s = 0; s < probe.nSites(); ++s) {
        xMin = std::min(xMin, probe.x(s));
        xMax = std::max(xMax, probe.x(s));
        yMin = std::min(yMin, probe.y(s));
        yMax = std::max(yMax, probe.y(s));
    }
    std::cout << "  x extent            " << std::fixed << std::setprecision(1) << xMin
              << " .. " << xMax << " um\n"
              << "  y extent            " << yMin << " .. " << yMax << " um\n";

    // Neighbourhood size drives detection cost, so report the distribution.
    probe.buildNeighborCache(detectRad);
    std::size_t minN = SIZE_MAX, maxN = 0, totalN = 0;
    for (int s = 0; s < probe.nSites(); ++s) {
        const std::size_t n = probe.cachedNeighbors(s).size();
        minN = std::min(minN, n);
        maxN = std::max(maxN, n);
        totalN += n;
    }
    std::cout << "  neighbours/site     " << minN << " min, "
              << std::setprecision(1)
              << static_cast<double>(totalN) / probe.nSites() << " mean, " << maxN
              << " max\n\n";

    std::cout << "first 8 sites (site -> channel, x, y, shank):\n";
    for (int s = 0; s < std::min(8, probe.nSites()); ++s) {
        std::cout << "  site " << std::setw(4) << (s + 1) << " -> chan "
                  << std::setw(4) << (probe.channelForSite(s) + 1) << "   ("
                  << std::setw(9) << probe.x(s) << ", " << std::setw(9) << probe.y(s)
                  << ")   shank " << probe.shank(s)
                  << (probe.isIgnored(s) ? "   [ignored]" : "") << "\n";
    }

    return 0;
}

// ------------------------------------------------------------------- detect --

int commandDetect(const std::string& prmPath, double limitSeconds, double chunkSeconds,
                  int recordingIndex, int siteReportFrom = 0, int siteReportTo = 0)
{
    jrc::Config config;
    config.load(prmPath);

    jrc::Probe probe;
    probe.load(config);

    auto settings = jrc::detect::DetectSettings::fromConfig(config);
    if (chunkSeconds > 0) settings.chunkSeconds = chunkSeconds;

    probe.buildNeighborCache(settings.peaks.evtDetectRad);
    probe.buildSiteNeighbors(settings.featureSettings.nSitesNeighborhood());

    if (config.rawRecordings().empty()) {
        std::cerr << "error: no rawRecordings listed in the parameter file\n";
        return 1;
    }

    std::cout << "detect: " << config.path() << "\n"
              << "  filter              " << filterTypeName(settings.filterType)
              << " (order " << settings.nDiffOrder << ")\n"
              << "  CAR                 " << config.getString("CARMode", "mean") << "\n"
              << "  qqFactor            " << settings.peaks.qqFactor << "\n"
              << "  refractory          " << settings.peaks.refracIntSamples
              << " samples\n"
              << "  evtDetectRad        " << settings.peaks.evtDetectRad << " um\n"
              << "  sites               " << probe.nSites() << " ("
              << probe.nIgnored() << " ignored)\n";

    const auto& fs = settings.featureSettings;
    std::cout << "  spike window        " << fs.windowSamples() << " samples ("
              << -fs.preSamples << " .. +" << fs.postSamples << ")\n"
              << "  sites per event     " << fs.nSitesEvt() << " of "
              << fs.nSitesNeighborhood() << "\n"
              << "  features            " << fs.nPCsPerSite << " PCs x "
              << fs.nSitesEvt() << " sites = " << fs.nFeatures() << " per peak, "
              << fs.nPeaksFeatures << " peak(s)\n\n";

    if (recordingIndex < 0 ||
        recordingIndex >= static_cast<int>(config.rawRecordings().size())) {
        std::cerr << "error: --recording " << recordingIndex << " is out of range (0.."
                  << config.rawRecordings().size() - 1 << ")\n";
        return 1;
    }
    const std::string& path = config.rawRecordings()[recordingIndex];

    jrc::Recording recording;
    recording.open(path, config);

    std::cout << "recording: " << path << "\n"
              << "  duration            " << std::fixed << std::setprecision(2)
              << recording.duration() << " s\n";

    if (limitSeconds > 0) {
        std::cout << "  limited to          " << limitSeconds << " s\n";
    }
    std::cout << "\n";

    jrc::detect::DetectController controller(config, probe, settings);

    const auto start = std::chrono::steady_clock::now();

    const auto progress = [](int c, int n, std::size_t spikes) {
        std::cout << "\r  chunk " << c << "/" << n << "   spikes " << spikes
                  << "        " << std::flush;
    };

    // --seconds truncates the range rather than the file.
    const std::int64_t lastSample =
        limitSeconds > 0
            ? std::min<std::int64_t>(
                  recording.nSamples(),
                  static_cast<std::int64_t>(limitSeconds * recording.sampleRate()))
            : recording.nSamples();

    const jrc::detect::DetectResult result =
        controller.runRange(recording, 0, lastSample, progress);

    const auto elapsed = std::chrono::duration<double>(
                             std::chrono::steady_clock::now() - start)
                             .count();

    std::cout << "\r" << std::string(60, ' ') << "\r";

    // Threshold summary over active sites only.
    double threshSum = 0.0;
    float threshMin = 1e30f, threshMax = 0.0f;
    int nActive = 0;
    for (int s = 0; s < probe.nSites(); ++s) {
        if (probe.isIgnored(s)) continue;
        const float t = result.siteThresholds[s];
        threshSum += t;
        threshMin = std::min(threshMin, t);
        threshMax = std::max(threshMax, t);
        ++nActive;
    }

    const double processedSeconds =
        limitSeconds > 0 ? std::min(limitSeconds, recording.duration())
                         : recording.duration();

    std::cout << "thresholds (raw ADC counts, active sites):\n"
              << "  min " << std::setprecision(1) << threshMin << "   mean "
              << (nActive ? threshSum / nActive : 0.0) << "   max " << threshMax
              << "\n\n"
              << "spikes detected       " << result.size() << "\n"
              << "  rate                " << std::setprecision(1)
              << (processedSeconds > 0 ? result.size() / processedSeconds : 0.0)
              << " spikes/s\n"
              << "elapsed               " << std::setprecision(2) << elapsed << " s\n"
              << "  throughput          " << std::setprecision(1)
              << (elapsed > 0 ? processedSeconds / elapsed : 0.0)
              << "x realtime\n";

    // Per-site diagnostic: why a site with visible spikes may produce none.
    // Counts both the primary site a spike was assigned to and the secondary,
    // so a site whose events were all attributed to a neighbour is visible as
    // "primary 0, secondary many" rather than looking simply dead.
    if (siteReportTo >= siteReportFrom && siteReportTo > 0) {
        const int nSites = probe.nSites();
        std::vector<std::int64_t> primary(nSites, 0), secondary(nSites, 0);
        for (std::size_t i = 0; i < result.spikeSites.size(); ++i) {
            const int s = result.spikeSites[i];
            if (s >= 0 && s < nSites) ++primary[s];
        }
        for (std::size_t i = 0; i < result.spikeSites2.size(); ++i) {
            const int s = result.spikeSites2[i];
            if (s >= 0 && s < nSites) ++secondary[s];
        }

        std::cout << "\nper-site report (1-based site ids)\n"
                  << "  site  state     thresh   neighbors   primary   secondary\n";
        for (int s = siteReportFrom - 1; s < siteReportTo && s < nSites; ++s) {
            if (s < 0) continue;
            const int* neigh = probe.siteNeighbors(s);
            int nUsable = 0;
            for (int k = 0; k < probe.nSiteNeighbors(); ++k) {
                const int n = neigh[k];
                if (n >= 0 && n < nSites && !probe.isIgnored(n)) ++nUsable;
            }
            std::cout << std::setw(6) << (s + 1) << "  "
                      << std::setw(8) << (probe.isIgnored(s) ? "IGNORED" : "active")
                      << std::setw(9) << std::setprecision(1) << std::fixed
                      << result.siteThresholds[s] << std::setw(12) << nUsable
                      << std::setw(10) << primary[s] << std::setw(12) << secondary[s]
                      << "\n";
        }
        std::cout.unsetf(std::ios::floatfield);
    }

    // Reading runs on a prefetch thread while the previous chunk is filtered
    // and searched, so these stages overlap and do NOT sum to the wall clock.
    // Percentages are therefore given against the busiest single stage.
    const auto& t = result.timings;
    const double compute = t.filter + t.car + t.threshold + t.peaks + t.features;
    const double bound = std::max(t.read, compute);

    const auto stage = [&](const char* name, double seconds) {
        std::cout << "  " << std::left << std::setw(18) << name << std::right
                  << std::setw(8) << std::setprecision(2) << std::fixed << seconds
                  << " s  " << std::setw(5) << std::setprecision(1)
                  << (bound > 0 ? 100.0 * seconds / bound : 0.0) << " %\n";
    };

    std::cout << "\nstage breakdown (read overlaps compute; they do not sum):\n";
    stage("read (I/O)", t.read);
    stage("filter", t.filter);
    stage("CAR", t.car);
    stage("threshold", t.threshold);
    stage("peaks", t.peaks);
    stage("features", t.features);
    std::cout << "  " << std::left << std::setw(18) << "compute total" << std::right
              << std::setw(8) << std::setprecision(2) << compute << " s\n"
              << "  " << std::left << std::setw(18) << "wall clock" << std::right
              << std::setw(8) << t.total << " s\n"
              << "  " << std::left << std::setw(18) << "limited by" << std::right
              << std::setw(8) << (t.read > compute ? "I/O" : "CPU") << "\n"
              << "  " << std::left << std::setw(18) << "threads" << std::right
              << std::setw(8) << jrc::threadCount() << "\n";

    if (recording.fileBytes() > 0 && t.read > 0) {
        std::cout << "  " << std::left << std::setw(18) << "read rate" << std::right
                  << std::setw(8) << std::setprecision(0)
                  << (static_cast<double>(recording.fileBytes()) / (1024.0 * 1024.0) *
                      (processedSeconds / recording.duration()) / t.read)
                  << " MB/s\n";
    }

    if (result.nFeatures > 0 && !result.spikeFeatures.empty()) {
        const auto& pv = result.principalVectors;

        std::cout << "\nfeatures:\n"
                  << "  per spike           " << result.nFeatures << " x "
                  << result.nPeaks << " peak(s) = "
                  << result.nFeatures * result.nPeaks << " floats\n"
                  << "  memory              " << std::setprecision(2)
                  << result.spikeFeatures.size() * sizeof(float) / (1024.0 * 1024.0)
                  << " MiB\n";

        if (!pv.values.empty()) {
            // Variance explained is the standard sanity check on a PCA: the
            // eigenvalues must be descending and the first few dominant.
            double totalVar = 0.0;
            for (double v : pv.values) totalVar += std::max(0.0, v);

            std::cout << "  variance explained  ";
            for (int k = 0; k < std::min<int>(3, static_cast<int>(pv.values.size())); ++k)
                std::cout << "PC" << (k + 1) << " " << std::setprecision(1)
                          << (totalVar > 0 ? 100.0 * pv.values[k] / totalVar : 0.0)
                          << "%  ";
            std::cout << "\n";

            // Each principal vector is z-scored, so its mean must be ~0 and its
            // sample standard deviation ~1.
            for (int k = 0; k < 1; ++k) {
                const float* v = pv.component(k);
                double sum = 0.0, sumsq = 0.0;
                for (int i = 0; i < pv.nSamples; ++i) {
                    sum += v[i];
                    sumsq += double(v[i]) * v[i];
                }
                const double mean = sum / pv.nSamples;
                const double sd = std::sqrt((sumsq - pv.nSamples * mean * mean) /
                                            std::max(1, pv.nSamples - 1));
                std::cout << "  PC1 z-score check   mean " << std::setprecision(3)
                          << mean << ", sd " << sd << " (expect 0, 1)\n";
            }
        }

        // Spread of the first feature, a quick check that projections are live.
        double mean = 0.0;
        for (std::size_t i = 0; i < result.size(); ++i) mean += result.featuresFor(i, 0)[0];
        mean /= std::max<std::size_t>(1, result.size());

        double variance = 0.0;
        for (std::size_t i = 0; i < result.size(); ++i) {
            const double d = result.featuresFor(i, 0)[0] - mean;
            variance += d * d;
        }
        variance /= std::max<std::size_t>(1, result.size());

        std::cout << "  feature[0] spread   mean " << std::setprecision(2) << mean
                  << ", sd " << std::sqrt(variance) << "\n";
    }

    if (!result.spikeTimes.empty()) {
        std::cout << "\nfirst 10 spikes (time_s, site, amplitude_uV):\n";
        for (std::size_t i = 0; i < std::min<std::size_t>(10, result.size()); ++i) {
            std::cout << "  " << std::setprecision(4)
                      << result.spikeTimes[i] / recording.sampleRate() << "   site "
                      << std::setw(4) << (result.spikeSites[i] + 1) << "   "
                      << std::setprecision(1)
                      << result.spikeAmps[i] * recording.bitScaling() << "\n";
        }
    }

    return 0;
}

// --------------------------------------------------------------------- sort --

/// Output prefix JRCLUST would use: outputDir (or the .prm's directory) plus
/// the session name, which is the .prm's base name.
/// Options shared by the sort command, so a stale .prm can be pointed at the
/// data it actually refers to without editing the file.
struct SortOptions {
    double limitSeconds = 0.0;
    int recordingIndex = 0;
    bool useGpu = true;
    bool save = false;
    bool autoMerge = true;
    bool removeRefractory = false;
    bool removeOutliers = false;
    bool refractoryMergeCheck = false;
    bool matlabLRatio = false;
    std::string recordingOverride;
    std::string outputOverride;

    /// -1 concatenates every entry of rawRecordings into one sort, which is
    /// what JRCLUST does. A non-negative value sorts just that one file.
    int singleRecording = -1;

    bool classify = true;          ///< label units noise / single / multi
    bool armDelete = false;        ///< write 'to_delete' instead of 'noise'
    bool labelNoiseOnly = false;   ///< annotate only suspect units, not all
    bool splitBimodal = true;      ///< separate clusters holding two populations
    bool vetoSpikes = false;       ///< reject non-spike-shaped events at detection
    double vetoDecayRatio = 0.65;  ///< far/near amplitude above which it is flat
    double maxVpp = 0.0;           ///< uV ceiling for the amplitude test
};

std::string outputPrefix(const jrc::Config& config, const std::string& overrideDir)
{
    std::string directory = overrideDir.empty() ? config.getString("outputDir", "")
                                                : overrideDir;
    if (directory.empty()) {
        const auto slash = config.path().find_last_of("/\\");
        directory = (slash == std::string::npos) ? "." : config.path().substr(0, slash);
    }

    std::string session = config.path();
    const auto slash = session.find_last_of("/\\");
    if (slash != std::string::npos) session = session.substr(slash + 1);
    const auto dot = session.find_last_of('.');
    if (dot != std::string::npos) session = session.substr(0, dot);

    return directory + "/" + session;
}

int commandSort(const std::string& prmPath, const SortOptions& opts)
{
    jrc::Config config;
    config.load(prmPath);

    if (!opts.recordingOverride.empty())
        config.overrideRecordings({opts.recordingOverride});

    jrc::Probe probe;
    probe.load(config);

    auto settings = jrc::detect::DetectSettings::fromConfig(config);
    settings.retainWaveforms = true;  // post-processing needs them
    settings.veto.enabled = opts.vetoSpikes;
    settings.veto.maxDecayRatio = opts.vetoDecayRatio;
    probe.buildNeighborCache(settings.peaks.evtDetectRad);
    probe.buildSiteNeighbors(settings.featureSettings.nSitesNeighborhood());

    if (config.rawRecordings().empty()) {
        std::cerr << "error: no rawRecordings listed in the parameter file\n";
        return 1;
    }
    // Which files go into this sort. JRCLUST treats every entry of
    // rawRecordings as one continuous recording; --recording picks just one.
    std::vector<std::string> paths;
    if (opts.singleRecording >= 0) {
        if (opts.singleRecording >= static_cast<int>(config.rawRecordings().size())) {
            std::cerr << "error: --recording out of range (0.."
                      << config.rawRecordings().size() - 1 << ")\n";
            return 1;
        }
        paths.push_back(config.rawRecordings()[opts.singleRecording]);
    } else {
        paths = config.rawRecordings();
    }

    const double limitSeconds = opts.limitSeconds;
    const bool save = opts.save;
    const bool useGpu = opts.useGpu;

    // Open everything up front so a bad path fails before any work is done.
    std::vector<jrc::Recording> recordings(paths.size());
    std::vector<std::int64_t> lastSamples(paths.size());
    std::int64_t totalSamples = 0;

    std::cout << "sort: " << config.path() << "\n"
              << (paths.size() > 1 ? "  mode                concatenated ("
                                   : "  mode                single (")
              << paths.size() << " recording" << (paths.size() == 1 ? "" : "s")
              << ")\n";

    for (std::size_t i = 0; i < paths.size(); ++i) {
        recordings[i].open(paths[i], config);

        lastSamples[i] =
            limitSeconds > 0
                ? std::min<std::int64_t>(recordings[i].nSamples(),
                                         static_cast<std::int64_t>(
                                             limitSeconds * recordings[i].sampleRate()))
                : recordings[i].nSamples();

        std::cout << "  [" << i << "] " << paths[i] << "\n"
                  << "      " << std::fixed << std::setprecision(2)
                  << lastSamples[i] / recordings[i].sampleRate() << " s\n";
        totalSamples += lastSamples[i];
    }

    const double sampleRate = recordings.front().sampleRate();
    std::cout << "  total duration      " << totalSamples / sampleRate << " s\n\n";

    // -- detect -------------------------------------------------------------
    jrc::detect::DetectController controller(config, probe, settings);

    const std::string prefix = outputPrefix(config, opts.outputOverride);
    jrc::io::ResultWriter writer;
    if (save) {
        writer.open(prefix, true);
        controller.setWriter(&writer);
        std::cout << "saving to             " << prefix << "_{filt,raw,features}.jrc\n\n";
    }

    const auto tDetect = std::chrono::steady_clock::now();

    // Per-file thresholds, as detect.m computes them: each session has its own
    // noise level, and a shared threshold would bias detection toward whichever
    // file is quieter.
    std::vector<std::vector<float>> thresholds(paths.size());
    for (std::size_t i = 0; i < paths.size(); ++i) {
        std::cout << "\r  thresholds  file " << i + 1 << "/" << paths.size()
                  << "            " << std::flush;
        thresholds[i] = controller.estimateThresholds(recordings[i], 0, lastSamples[i]);
    }

    // ONE principal-vector basis pooled across every file. Features projected
    // onto different bases are not comparable, so per-file bases would make
    // cross-file clustering meaningless - this is the part that actually makes
    // a concatenated sort valid.
    jrc::features::PrincipalVectors basis;
    if (settings.computeFeatures) {
        jrc::features::SpikeWindows pooled;
        pooled.nSamples = settings.featureSettings.windowSamples();
        pooled.nSites = 1;
        pooled.nSpikes = 0;

        for (std::size_t i = 0; i < paths.size(); ++i) {
            std::cout << "\r  PCA basis   file " << i + 1 << "/" << paths.size()
                      << "            " << std::flush;
            controller.collectPrincipalVectorSamples(recordings[i], 0, lastSamples[i],
                                                     thresholds[i], pooled);
        }
        basis = jrc::features::computePrincipalVectors(pooled,
                                                       settings.featureSettings);
    }

    jrc::detect::DetectResult detected;
    std::int64_t offset = 0;
    std::vector<std::int64_t> fileOffsets(paths.size(), 0);

    for (std::size_t i = 0; i < paths.size(); ++i) {
        controller.setPrecomputed(&thresholds[i],
                                  settings.computeFeatures ? &basis : nullptr);

        const std::size_t before = detected.size();
        const jrc::detect::DetectResult part = controller.runRange(
            recordings[i], 0, lastSamples[i],
            [&](int c, int n, std::size_t spikes) {
                std::cout << "\r  detect  file " << i + 1 << "/" << paths.size()
                          << "  chunk " << c << "/" << n << "   spikes "
                          << before + spikes << "        " << std::flush;
            });

        fileOffsets[i] = offset;
        jrc::detect::appendDetectResult(detected, part, offset);
        offset += lastSamples[i];
    }

    controller.setPrecomputed(nullptr, nullptr);

    const double detectSeconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - tDetect).count();

    std::cout << "\r" << std::string(72, ' ') << "\r";
    std::cout << "detected              " << detected.size() << " spikes in "
              << std::setprecision(1) << detectSeconds << " s\n";

    {
        const auto& v = detected.veto;
        if (v.gapSamples > 0 || v.inGap > 0) {
            std::cout << "  recording gaps      " << v.inGap
                      << " events rejected inside zero-filled gaps ("
                      << std::setprecision(1)
                      << v.gapSamples / recordings.front().sampleRate() << " s blanked)\n";
        }
        if (opts.vetoSpikes) {
            std::cout << "  shape veto          " << v.rejected << " of " << v.tested
                      << " events rejected ("
                      << (v.tested ? 100.0 * v.rejected / v.tested : 0.0) << " %)\n"
                      << "    too narrow        " << v.tooNarrow << "\n"
                      << "    single-sample     " << v.tooSpiky << "\n"
                      << "    no spatial decay  " << v.tooFlat << "\n";
        }
    }
    if (paths.size() > 1) {
        std::cout << "  spike times are on one timeline spanning all "
                  << paths.size() << " files\n";
        for (std::size_t i = 0; i < paths.size(); ++i)
            std::cout << "      file " << i << " starts at sample " << fileOffsets[i]
                      << " (" << std::setprecision(1)
                      << fileOffsets[i] / sampleRate << " s)\n";
    }
    std::cout << "\n";

    if (detected.size() < 2) {
        std::cout << "too few spikes to cluster\n";
        return 0;
    }

    // -- cluster ------------------------------------------------------------
    auto clusterSettings = jrc::sort::ClusterSettings::fromConfig(config);
    clusterSettings.useGpu = useGpu;

    std::cout << "clustering (" << (useGpu ? "GPU" : "CPU") << "):\n"
              << "  distCut             " << clusterSettings.distCut << " %\n"
              << "  log10RhoCut         " << clusterSettings.log10RhoCut << "\n"
              << "  log10DeltaCut       " << clusterSettings.log10DeltaCut << "\n"
              << "  minClusterSize      "
              << std::max(clusterSettings.minClusterSize, 2 * detected.nFeatures)
              << "\n";

    jrc::sort::ClusterResult clusters =
        jrc::sort::cluster(detected, probe, clusterSettings,
                           [&](int site, int nSites, int) {
                               if (site % 16 == 0 || site == nSites)
                                   std::cout << "\r  site " << site << "/" << nSites
                                             << "        " << std::flush;
                           });

    std::cout << "\r" << std::string(64, ' ') << "\r";

    const int nBeforeMerge = clusters.nClusters;
    std::cout << "  raw clusters        " << nBeforeMerge << "\n";

    // -- post-process --------------------------------------------------------
    auto postSettings = jrc::sort::PostProcessSettings::fromConfig(config);
    postSettings.refractoryMergeCheck = opts.refractoryMergeCheck;
    postSettings.matlabLRatioDof = opts.matlabLRatio;

    // setprecision is sticky and an earlier line left it at 1 decimal, which
    // rendered maxUnitSim 0.98 as "1.0".
    std::cout << "\npost-processing:\n"
              << std::setprecision(3)
              << "  autoMergeBy         " << postSettings.autoMergeBy << "\n"
              << "  maxUnitSim          " << postSettings.maxUnitSim << "\n"
              << "  nPassesMerge        " << postSettings.nPassesMerge << "\n"
              << "  driftMerge          " << (postSettings.driftMerge ? "on" : "off")
              << "\n"
              << "  auto-merge          " << (opts.autoMerge ? "on" : "off") << "\n"
              << "  refractory cleanup  " << (opts.removeRefractory ? "on" : "off")
              << "\n"
              << "  refractory merge gate  "
              << (opts.refractoryMergeCheck ? "on" : "off") << "\n";

    const jrc::sort::PostProcessResult post = jrc::sort::runPostProcessing(
        detected, clusters, probe, postSettings, opts.autoMerge,
        opts.removeRefractory, opts.removeOutliers,
        [&](const char* stage, int pass, int nPasses) {
            if (nPasses > 0)
                std::cout << "\r  " << stage << " (pass " << pass << "/" << nPasses
                          << ")            " << std::flush;
            else
                std::cout << "\r  " << stage << "                       " << std::flush;
        });

    std::cout << "\r" << std::string(64, ' ') << "\r";

    if (opts.autoMerge) {
        std::cout << "  merged              " << nBeforeMerge << " -> "
                  << clusters.nClusters << " (" << post.nMergedAway
                  << " merged away in " << post.nMergePasses << " pass"
                  << (post.nMergePasses == 1 ? "" : "es") << ")\n";
        if (opts.refractoryMergeCheck)
            std::cout << "  merges vetoed       " << post.nMergesVetoed
                      << " (spike trains fire independently)\n";
    }
    std::cout << "  post-process time   " << std::setprecision(2) << post.seconds
              << " s\n";

    std::size_t assigned = 0;
    for (std::int32_t c : clusters.spikeClusters)
        if (c > 0) ++assigned;

    std::cout << "\nresults:\n"
              << "  clusters            " << clusters.nClusters << "\n"
              << "  assigned spikes     " << assigned << " ("
              << std::setprecision(1)
              << (detected.size() ? 100.0 * assigned / detected.size() : 0.0) << " %)\n"
              << "  unassigned          " << detected.size() - assigned << "\n";

    if (clusters.nClusters > 0) {
        auto counts = clusters.unitCount;
        std::sort(counts.begin(), counts.end(), std::greater<std::int32_t>());
        const std::int32_t median = counts[counts.size() / 2];
        std::cout << "  cluster size        " << counts.back() << " min, " << median
                  << " median, " << counts.front() << " max\n";

        // Quality summary. Medians rather than means: these distributions have
        // long tails and a couple of multi-unit clusters would skew a mean.
        const auto& m = post.metrics;
        const auto medianOf = [](std::vector<double> v) {
            v.erase(std::remove_if(v.begin(), v.end(),
                                   [](double x) { return !std::isfinite(x); }),
                    v.end());
            if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
            std::sort(v.begin(), v.end());
            return v[v.size() / 2];
        };
        const auto medianOfF = [&](const std::vector<float>& v) {
            return medianOf(std::vector<double>(v.begin(), v.end()));
        };

        std::cout << "\nquality (median over units):\n"
                  << "  SNR                 " << std::setprecision(2)
                  << medianOfF(m.unitSNR) << "\n"
                  << "  Vpp                 " << medianOfF(m.unitVpp) << " uV\n"
                  << "  ISI ratio           " << medianOf(m.unitISIRatio) << "\n"
                  << "  isolation distance  " << medianOf(m.unitIsoDist) << "\n"
                  << "  L-ratio             " << medianOf(m.unitLRatio) << "\n";

        int wellIsolated = 0;
        for (int c = 0; c < clusters.nClusters; ++c) {
            const bool goodSnr = c < static_cast<int>(m.unitSNR.size()) &&
                                 m.unitSNR[c] >= 5.0;
            const bool cleanIsi = c < static_cast<int>(m.unitISIRatio.size()) &&
                                  std::isfinite(m.unitISIRatio[c]) &&
                                  m.unitISIRatio[c] <= 0.1;
            if (goodSnr && cleanIsi) ++wellIsolated;
        }
        std::cout << "  units with SNR>=5 and ISI ratio<=0.1: " << wellIsolated
                  << " / " << clusters.nClusters << "\n";

        // The point of a concatenated sort is units followed across sessions,
        // so report how many actually span them rather than assuming they do.
        if (paths.size() > 1) {
            std::vector<std::vector<std::int64_t>> perFile(
                clusters.nClusters, std::vector<std::int64_t>(paths.size(), 0));

            for (std::size_t i = 0; i < detected.size(); ++i) {
                const std::int32_t c = clusters.spikeClusters[i];
                if (c <= 0) continue;

                // Which file this spike's global time falls in.
                std::size_t file = 0;
                for (std::size_t f = 0; f + 1 < paths.size(); ++f)
                    if (detected.spikeTimes[i] >= fileOffsets[f + 1]) file = f + 1;
                ++perFile[c - 1][file];
            }

            int spanAll = 0;
            std::vector<int> presentIn(paths.size() + 1, 0);
            for (int c = 0; c < clusters.nClusters; ++c) {
                int nFilesPresent = 0;
                for (std::size_t f = 0; f < paths.size(); ++f)
                    if (perFile[c][f] > 0) ++nFilesPresent;
                ++presentIn[nFilesPresent];
                if (nFilesPresent == static_cast<int>(paths.size())) ++spanAll;
            }

            std::cout << "\ncross-session:\n"
                      << "  units in all " << paths.size() << " sessions: " << spanAll
                      << " / " << clusters.nClusters << "\n";
            for (std::size_t f = 1; f <= paths.size(); ++f)
                std::cout << "  units present in exactly " << f << " session"
                          << (f == 1 ? "" : "s") << ": " << presentIn[f] << "\n";

            std::cout << "\n  largest 5 units, spikes per session:\n";
            std::vector<int> byCount(clusters.nClusters);
            std::iota(byCount.begin(), byCount.end(), 0);
            std::sort(byCount.begin(), byCount.end(), [&](int a, int b) {
                return clusters.unitCount[a] > clusters.unitCount[b];
            });
            for (int k = 0; k < std::min(5, clusters.nClusters); ++k) {
                const int c = byCount[k];
                std::cout << "    unit " << std::setw(4) << (c + 1) << ":";
                for (std::size_t f = 0; f < paths.size(); ++f)
                    std::cout << std::setw(9) << perFile[c][f];
                std::cout << "\n";
            }
        }
    }

    std::cout << "\ntimings:\n"
              << "  detect              " << std::setprecision(2) << detectSeconds
              << " s\n"
              << "  rho (kernel)        " << clusters.timings.rho << " s\n"
              << "  delta (kernel)      " << clusters.timings.delta << " s\n"
              << "  assign              " << clusters.timings.assign << " s\n"
              << "  cluster total       " << clusters.timings.total << " s\n"
              << "  TOTAL               " << detectSeconds + clusters.timings.total
              << " s\n";

    // -- classify units ------------------------------------------------------
    // After merging, so the verdict is on final units rather than fragments.
    jrc::sort::QualityClassification classification;
    std::vector<std::string> clusterNotes;

    if (opts.classify && clusters.nClusters > 0) {
        jrc::sort::QualityCriteria criteria;
        criteria.maxVpp = opts.maxVpp;
        criteria.labelGoodUnits = !opts.labelNoiseOnly;
        if (opts.armDelete) criteria.noiseLabel = "to_delete";

        // Divide clusters that hold two populations before judging them. A
        // cluster that is mostly neuron has a mostly-neuronal mean waveform
        // however bad its contaminating minority, so it must be separated
        // before any whole-cluster test can see the artifact.
        jrc::sort::MeanWaveforms waveformsForClassify = post.waveforms;
        jrc::sort::UnitMetrics metricsForClassify = post.metrics;

        if (opts.splitBimodal) {
            const jrc::sort::SplitCriteria splitCriteria;
            const auto split = jrc::sort::splitBimodalClusters(
                detected, clusters, post.waveforms, post.metrics, probe, postSettings,
                criteria, splitCriteria);

            std::cout << "\nbimodal split:\n"
                      << "  clusters split      " << split.nSplit << " / "
                      << split.nClustersTested << "\n"
                      << "  spikes separated    " << split.nSpikesMovedOut << "\n";

            if (split.nSplit > 0) {
                // Mean waveforms and metrics describe the old clustering, so
                // recompute them against the new one before labelling.
                waveformsForClassify = jrc::sort::computeMeanWaveforms(
                    detected, clusters, probe, postSettings);
                metricsForClassify = jrc::sort::computeQualityScores(
                    detected, clusters, waveformsForClassify, probe, postSettings);
                std::cout << "  clusters now        " << clusters.nClusters << "\n";
            }
        }

        classification = jrc::sort::classifyUnits(detected, clusters,
                                                  waveformsForClassify,
                                                  metricsForClassify, probe,
                                                  postSettings, criteria);
        clusterNotes = jrc::sort::annotationStrings(classification, criteria);

        std::cout << "\nunit classification (labels go to the curation GUI):\n"
                  << "  0 noise             " << classification.nNoise << "\n"
                  << "  1 single unit       " << classification.nSingle << "\n"
                  << "  2 multi unit        " << classification.nMulti << "\n";

        if (opts.armDelete)
            std::cout << "  noise labelled 'to_delete' -> Edit > Delete annotated\n";
        else
            std::cout << "  noise labelled 'noise' (use --arm-delete for one-click "
                         "removal)\n";

        if (!opts.labelNoiseOnly)
            std::cout << "  every unit carries a label, so the GUI's 'Display\n"
                         "  annotated units only' matches all of them\n"
                         "  (--label-noise-only labels just the suspect units)\n";
        else if (classification.nNoise == 0)
            std::cout << "  no units annotated: the GUI will say 'Subset is empty'\n";
        else
            std::cout << "  only these " << classification.nNoise
                      << " are annotated, so 'Display annotated units only'\n"
                         "  shows exactly them\n";

        // Show why, so the thresholds can be judged rather than trusted.
        std::map<std::string, int> reasonCounts;
        for (const auto& verdict : classification.verdicts)
            if (verdict.quality == jrc::sort::UnitQuality::Noise)
                for (const auto& reason : verdict.reasons) ++reasonCounts[reason];

        if (!reasonCounts.empty()) {
            std::cout << "\n  why units were called noise:\n";
            std::vector<std::pair<std::string, int>> sorted(reasonCounts.begin(),
                                                            reasonCounts.end());
            std::sort(sorted.begin(), sorted.end(),
                      [](const auto& a, const auto& b) { return a.second > b.second; });
            for (const auto& entry : sorted)
                std::cout << "    " << std::setw(4) << entry.second << "  "
                          << entry.first << "\n";
        }

        // Percentiles of each metric, so the thresholds can be judged against
        // the actual distribution instead of taken on faith.
        const auto percentiles = [&](const char* name,
                                     const std::function<double(int)>& get,
                                     int width) {
            std::vector<double> v;
            for (int c = 0; c < clusters.nClusters; ++c) {
                const double x = get(c);
                if (std::isfinite(x)) v.push_back(x);
            }
            if (v.empty()) return;
            std::sort(v.begin(), v.end());
            const auto at = [&](double p) {
                return v[std::min(v.size() - 1,
                                  static_cast<std::size_t>(p * v.size()))];
            };
            std::cout << "    " << std::left << std::setw(16) << name << std::right
                      << std::setprecision(width) << std::setw(10) << at(0.05)
                      << std::setw(10) << at(0.25) << std::setw(10) << at(0.50)
                      << std::setw(10) << at(0.75) << std::setw(10) << at(0.95)
                      << "\n";
        };

        std::cout << "\n  metric distribution across units:\n"
                  << "    " << std::left << std::setw(16) << "" << std::right
                  << std::setw(10) << "p5" << std::setw(10) << "p25" << std::setw(10)
                  << "median" << std::setw(10) << "p75" << std::setw(10) << "p95"
                  << "\n";
        percentiles("Vpp (uV)", [&](int c) { return classification.verdicts[c].shape.vpp; }, 1);
        percentiles("half-width ms", [&](int c) { return classification.verdicts[c].shape.halfWidthMs; }, 3);
        percentiles("spikiness", [&](int c) { return classification.verdicts[c].shape.spikiness; }, 2);
        percentiles("recovery", [&](int c) { return classification.verdicts[c].shape.recoveryRatio; }, 2);
        percentiles("sites", [&](int c) { return double(classification.verdicts[c].shape.spatialSpread); }, 1);
        percentiles("SNR", [&](int c) { return double(post.metrics.unitSNR[c]); }, 2);
        percentiles("ISI ratio", [&](int c) { return post.metrics.unitISIRatio[c]; }, 3);

        // A couple of examples, so the numbers behind a verdict are visible.
        int shown = 0;
        for (int c = 0; c < clusters.nClusters && shown < 3; ++c) {
            const auto& verdict = classification.verdicts[c];
            if (verdict.quality != jrc::sort::UnitQuality::Noise) continue;
            const auto& s = verdict.shape;
            std::cout << (shown == 0 ? "\n  example noise units:\n" : "")
                      << "    unit " << std::setw(4) << (c + 1) << "  Vpp "
                      << std::setw(7) << std::setprecision(1) << s.vpp
                      << " uV  half-width " << std::setprecision(3) << s.halfWidthMs
                      << " ms  spikiness " << std::setprecision(2) << s.spikiness
                      << "  sites " << s.spatialSpread << "\n";
            ++shown;
        }
    }

    if (save) {
        writer.close();

        jrc::io::SidecarData sidecar;
        sidecar.nSpikes = static_cast<std::int64_t>(detected.size());
        sidecar.nSites = probe.nSites();
        sidecar.nClusters = clusters.nClusters;
        sidecar.nFeatures = detected.nFeatures;
        sidecar.nPeaks = detected.nPeaks;
        sidecar.filtSamples = writer.filtSamples();
        sidecar.rawSamples = writer.rawSamples();
        sidecar.nSiteNeighbors = writer.nSiteNeighbors();
        sidecar.sampleRate = sampleRate;
        sidecar.detectSeconds = detectSeconds;
        sidecar.sortSeconds = clusters.timings.total;

        sidecar.spikeTimes = &detected.spikeTimes;
        sidecar.spikeAmps = &detected.spikeAmps;
        sidecar.spikeSites = &detected.spikeSites;
        sidecar.spikeSites2 = &detected.spikeSites2;
        sidecar.spikeRho = &clusters.spikeRho;
        sidecar.spikeDelta = &clusters.spikeDelta;
        sidecar.spikeNeigh = &clusters.spikeNeigh;
        sidecar.spikeClusters = &clusters.spikeClusters;
        sidecar.clusterCenters = &clusters.clusterCenters;
        sidecar.clusterSites = &clusters.clusterSites;
        sidecar.unitCount = &clusters.unitCount;
        sidecar.siteThresh = &detected.siteThresholds;
        sidecar.rhoCutSite = &clusters.siteRhoCut;
        sidecar.spikePositions = &detected.spikePositions;
        if (!clusterNotes.empty()) sidecar.clusterNotes = &clusterNotes;

        jrc::io::writeSidecar(prefix, sidecar);

        std::cout << "\nsaved:\n"
                  << "  " << prefix << "_filt.jrc      [" << writer.filtSamples()
                  << " x " << writer.nSiteNeighbors() << " x " << detected.size()
                  << "] int16\n"
                  << "  " << prefix << "_raw.jrc       [" << writer.rawSamples() << " x "
                  << writer.nSiteNeighbors() << " x " << detected.size()
                  << "] int16\n"
                  << "  " << prefix << "_features.jrc  [" << detected.nFeatures << " x "
                  << detected.nPeaks << " x " << detected.size() << "] single\n"
                  << "  " << prefix << "_cpp.bin       spike table + clustering\n"
                  << "\nnext, in MATLAB:\n"
                  << "  jrc_import_cpp('" << config.path() << "')\n"
                  << "  jrc manual " << config.path() << "\n";
    }

    if (clusters.nClusters > 0) {
        std::cout << "\nlargest 10 clusters (id, spikes, site, firing rate Hz):\n";

        std::vector<int> order(clusters.nClusters);
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](int a, int b) {
            return clusters.unitCount[a] > clusters.unitCount[b];
        });

        // Rates are over the whole concatenation, not one file.
        const double seconds = static_cast<double>(totalSamples) / sampleRate;

        for (int k = 0; k < std::min(10, clusters.nClusters); ++k) {
            const int c = order[k];
            std::cout << "  " << std::setw(4) << (c + 1) << "   " << std::setw(8)
                      << clusters.unitCount[c] << "   site " << std::setw(4)
                      << (clusters.clusterSites[c] + 1) << "   " << std::setprecision(1)
                      << (seconds > 0 ? clusters.unitCount[c] / seconds : 0.0) << "\n";
        }
    }

    return 0;
}

// ----------------------------------------------------------------- selftest --

/// Builds a synthetic site: `nFeatures` x n12 features plus time and rho orders.
struct SyntheticSite {
    std::vector<float> features;
    std::vector<std::int32_t> spikeOrder;
    std::vector<std::int32_t> rhoOrder;
    jrc::cuda::SiteConstants consts;
};

SyntheticSite makeSyntheticSite(int n1, int nFeatures)
{
    SyntheticSite site;

    site.consts.n1 = n1;
    site.consts.n12 = n1 * 2;  // half the spikes come from neighbouring sites
    site.consts.nFeatures = nFeatures;
    site.consts.dnMax = site.consts.n12 / 4;  // as nClusterIntervals = 4 would give
    site.consts.fDcSpk = 0;

    const int n12 = site.consts.n12;

    std::mt19937 rng(20260903);
    std::normal_distribution<float> feature(0.0f, 1.0f);

    site.features.resize(static_cast<std::size_t>(nFeatures) * n12);
    for (auto& value : site.features) value = feature(rng);

    // Orders are permutations of 0..n12-1, matching rankorder() output.
    site.spikeOrder.resize(n12);
    site.rhoOrder.resize(n12);
    std::iota(site.spikeOrder.begin(), site.spikeOrder.end(), 0);
    std::iota(site.rhoOrder.begin(), site.rhoOrder.end(), 0);
    std::shuffle(site.rhoOrder.begin(), site.rhoOrder.end(), rng);

    return site;
}

/// Checks the waveform-shape metrics on waveforms whose answer we know.
///
/// The classifier decides what is and is not a neuron, so the shape measures
/// behind that decision need testing against a known case rather than assumed
/// correct from a recording that happens to contain no artifacts.
bool selfTestShapeMetrics()
{
    constexpr double kSampleRate = 20000.0;
    constexpr int kSamples = 41;
    constexpr int kSites = 4;

    const auto build = [&](const std::vector<float>& trace) {
        jrc::sort::MeanWaveforms w;
        w.nSamples = kSamples;
        w.nSamplesRaw = kSamples;
        w.nSites = kSites;
        w.nLocalSites = kSites;
        w.nClusters = 1;
        w.globalRaw.assign(static_cast<std::size_t>(kSites) * kSamples, 0.0f);
        w.global = w.globalRaw;
        // Put the waveform on site 0 and a half-amplitude copy on site 1, so
        // the spatial-spread count has something realistic to measure.
        for (int s = 0; s < kSamples; ++s) {
            w.globalRaw[s] = trace[s];
            w.globalRaw[static_cast<std::size_t>(kSamples) + s] = 0.5f * trace[s];
        }
        w.global = w.globalRaw;
        return w;
    };

    // A plausible action potential: a trough a few samples wide with a
    // positive rebound after it.
    std::vector<float> spike(kSamples, 0.0f);
    for (int s = 0; s < kSamples; ++s) {
        const double t = (s - 10) / 4.0;           // trough at sample 10
        spike[s] = static_cast<float>(-100.0 * std::exp(-t * t) +
                                      40.0 * std::exp(-std::pow((s - 18) / 5.0, 2)));
    }

    // A one-sample digital glitch: flat, one huge excursion, flat again.
    std::vector<float> glitch(kSamples, 0.0f);
    glitch[10] = -1000.0f;

    const auto spikeShape =
        jrc::sort::measureShape(build(spike), 0, 0, kSampleRate, 0.25);
    const auto glitchShape =
        jrc::sort::measureShape(build(glitch), 0, 0, kSampleRate, 0.25);

    std::cout << "shape metrics:\n"
              << std::fixed << std::setprecision(3)
              << "  real spike    half-width " << spikeShape.halfWidthMs
              << " ms   spikiness " << std::setprecision(2) << spikeShape.spikiness
              << "   recovery " << spikeShape.recoveryRatio << "\n"
              << std::setprecision(3)
              << "  glitch        half-width " << glitchShape.halfWidthMs
              << " ms   spikiness " << std::setprecision(2) << glitchShape.spikiness
              << "   recovery " << glitchShape.recoveryRatio << "\n";

    jrc::sort::QualityCriteria criteria;
    const bool spikePasses = !(spikeShape.spikiness > criteria.glitchSpikiness &&
                               spikeShape.halfWidthMs <= criteria.glitchHalfWidthMs);
    const bool glitchCaught = glitchShape.spikiness > criteria.glitchSpikiness &&
                              glitchShape.halfWidthMs <= criteria.glitchHalfWidthMs;

    std::cout << "  real spike kept    " << (spikePasses ? "PASS" : "FAIL") << "\n"
              << "  glitch flagged     " << (glitchCaught ? "PASS" : "FAIL") << "\n";

    return spikePasses && glitchCaught;
}

int commandSelfTest(int nSpikes)
{
    const int nFeatures = 12;  // within the kernels' MAXDIM of 45
    const float distCut2 = 4.0f;

    std::cout << "self test: " << nSpikes << " primary spikes, " << nFeatures
              << " features\n\n";

    const auto devices = jrc::cuda::enumerateDevices();
    if (devices.empty()) {
        std::cout << "No CUDA device available - cannot run the kernel test.\n";
        return 1;
    }
    std::cout << "using device 0: " << devices.front().name << "\n\n";

    const SyntheticSite site = makeSyntheticSite(nSpikes, nFeatures);
    const int n1 = site.consts.n1;

    // -- rho ----------------------------------------------------------------
    std::vector<float> rhoGpu(n1), rhoCpu(n1);

    jrc::cuda::computeRhoSite(site.features.data(), site.spikeOrder.data(),
                              site.consts, distCut2, rhoGpu.data());
    jrc::cuda::computeRhoSiteCpu(site.features.data(), site.spikeOrder.data(),
                                 site.consts, distCut2, rhoCpu.data());

    double rhoMaxDiff = 0.0;
    for (int i = 0; i < n1; ++i)
        rhoMaxDiff = std::max(rhoMaxDiff, std::fabs(double(rhoGpu[i]) - rhoCpu[i]));

    const double rhoMean =
        std::accumulate(rhoCpu.begin(), rhoCpu.end(), 0.0) / std::max(1, n1);

    std::cout << "rho   mean " << std::fixed << std::setprecision(6) << rhoMean
              << "   max |gpu-cpu| " << std::scientific << rhoMaxDiff << "\n";

    // -- delta --------------------------------------------------------------
    std::vector<float> deltaGpu(n1), deltaCpu(n1);
    std::vector<std::int32_t> neighGpu(n1), neighCpu(n1);

    jrc::cuda::computeDeltaSite(site.features.data(), site.spikeOrder.data(),
                                site.rhoOrder.data(), site.consts, distCut2,
                                deltaGpu.data(), neighGpu.data());
    jrc::cuda::computeDeltaSiteCpu(site.features.data(), site.spikeOrder.data(),
                                   site.rhoOrder.data(), site.consts, distCut2,
                                   deltaCpu.data(), neighCpu.data());

    double deltaMaxDiff = 0.0;
    int neighMismatches = 0;
    for (int i = 0; i < n1; ++i) {
        deltaMaxDiff = std::max(deltaMaxDiff, std::fabs(double(deltaGpu[i]) - deltaCpu[i]));
        if (neighGpu[i] != neighCpu[i]) ++neighMismatches;
    }

    std::cout << "delta max |gpu-cpu| " << std::scientific << deltaMaxDiff
              << "   neighbour mismatches " << neighMismatches << " / " << n1 << "\n\n";

    // Both paths accumulate in a different order, so exact equality is not
    // expected; anything at float rounding level means the kernel is correct.
    const bool rhoOk = rhoMaxDiff < 1e-6;
    const bool deltaOk = deltaMaxDiff < 1e-3;
    const bool neighOk = neighMismatches == 0;

    std::cout << "rho       " << (rhoOk ? "PASS" : "FAIL") << "\n"
              << "delta     " << (deltaOk ? "PASS" : "FAIL") << "\n"
              << "neighbour " << (neighOk ? "PASS" : "FAIL") << "\n";

    std::cout << "\n";
    const bool shapeOk = selfTestShapeMetrics();

    const bool passed = rhoOk && deltaOk && neighOk && shapeOk;
    std::cout << "\n" << (passed ? "SELF TEST PASSED" : "SELF TEST FAILED") << "\n";

    return passed ? 0 : 1;
}

}  // namespace

int main(int argc, char** argv)
{
    if (argc < 2) {
        printUsage();
        return 1;
    }

    const std::string command = argv[1];

    try {
        if (command == "gpuinfo") return commandGpuInfo();

        if (command == "info") {
            if (argc < 3) {
                std::cerr << "error: 'info' needs a .prm file\n";
                return 1;
            }
            return commandInfo(argv[2]);
        }

        if (command == "probe") {
            if (argc < 3) {
                std::cerr << "error: 'probe' needs a .prm file\n";
                return 1;
            }
            return commandProbe(argv[2]);
        }

        if (command == "detect") {
            if (argc < 3) {
                std::cerr << "error: 'detect' needs a .prm file\n";
                return 1;
            }

            double limitSeconds = 0.0;
            double chunkSeconds = 0.0;
            int recordingIndex = 0;
            int siteReportFrom = 0, siteReportTo = 0;
            for (int i = 3; i + 1 < argc; ++i) {
                if (std::strcmp(argv[i], "--seconds") == 0)
                    limitSeconds = std::atof(argv[++i]);
                else if (std::strcmp(argv[i], "--chunk") == 0)
                    chunkSeconds = std::atof(argv[++i]);
                else if (std::strcmp(argv[i], "--recording") == 0)
                    recordingIndex = std::atoi(argv[++i]);
                else if (std::strcmp(argv[i], "--site-report") == 0) {
                    // "from:to", 1-based inclusive
                    const std::string spec = argv[++i];
                    const auto colon = spec.find(':');
                    siteReportFrom = std::atoi(spec.substr(0, colon).c_str());
                    siteReportTo = colon == std::string::npos
                                       ? siteReportFrom
                                       : std::atoi(spec.substr(colon + 1).c_str());
                }
            }

            return commandDetect(argv[2], limitSeconds, chunkSeconds, recordingIndex,
                                 siteReportFrom, siteReportTo);
        }

        if (command == "sort") {
            if (argc < 3) {
                std::cerr << "error: 'sort' needs a .prm file\n";
                return 1;
            }

            SortOptions opts;
            for (int i = 3; i < argc; ++i) {
                if (std::strcmp(argv[i], "--seconds") == 0 && i + 1 < argc)
                    opts.limitSeconds = std::atof(argv[++i]);
                else if (std::strcmp(argv[i], "--recording") == 0 && i + 1 < argc)
                    opts.singleRecording = std::atoi(argv[++i]);
                else if (std::strcmp(argv[i], "--recording-path") == 0 && i + 1 < argc)
                    opts.recordingOverride = argv[++i];
                else if (std::strcmp(argv[i], "--output-dir") == 0 && i + 1 < argc)
                    opts.outputOverride = argv[++i];
                else if (std::strcmp(argv[i], "--cpu") == 0)
                    opts.useGpu = false;
                else if (std::strcmp(argv[i], "--save") == 0)
                    opts.save = true;
                else if (std::strcmp(argv[i], "--no-merge") == 0)
                    opts.autoMerge = false;
                else if (std::strcmp(argv[i], "--rm-refrac") == 0)
                    opts.removeRefractory = true;
                else if (std::strcmp(argv[i], "--rm-outliers") == 0)
                    opts.removeOutliers = true;
                else if (std::strcmp(argv[i], "--no-split") == 0)
                    opts.splitBimodal = false;
                else if (std::strcmp(argv[i], "--veto-spikes") == 0)
                    opts.vetoSpikes = true;
                else if (std::strcmp(argv[i], "--veto-decay") == 0 && i + 1 < argc)
                    opts.vetoDecayRatio = std::atof(argv[++i]);
                else if (std::strcmp(argv[i], "--refrac-merge") == 0)
                    opts.refractoryMergeCheck = true;
                else if (std::strcmp(argv[i], "--matlab-lratio") == 0)
                    opts.matlabLRatio = true;
                else if (std::strcmp(argv[i], "--no-classify") == 0)
                    opts.classify = false;
                else if (std::strcmp(argv[i], "--arm-delete") == 0 ||
                         std::strcmp(argv[i], "--mark-for-delete") == 0)
                    opts.armDelete = true;
                else if (std::strcmp(argv[i], "--label-noise-only") == 0)
                    opts.labelNoiseOnly = true;
                else if (std::strcmp(argv[i], "--max-vpp") == 0 && i + 1 < argc)
                    opts.maxVpp = std::atof(argv[++i]);
            }

            return commandSort(argv[2], opts);
        }

        if (command == "selftest") {
            const int nSpikes = argc > 2 ? std::atoi(argv[2]) : 4096;
            return commandSelfTest(nSpikes);
        }

        if (command == "-h" || command == "--help" || command == "help") {
            printUsage();
            return 0;
        }

        std::cerr << "error: unknown command '" << command << "'\n\n";
        printUsage();
        return 1;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << "\n";
        return 1;
    }
}

