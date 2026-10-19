#include "detect/Filter.h"

#include "core/Parallel.h"

#include <algorithm>
#include <stdexcept>

namespace jrc {
namespace detect {

FilterType filterTypeFromString(const std::string& name)
{
    if (name == "ndiff") return FilterType::NDiff;
    if (name == "sgdiff") return FilterType::SgDiff;
    if (name == "bandpass") return FilterType::Bandpass;
    if (name == "fir1") return FilterType::Fir1;
    if (name == "user") return FilterType::User;
    if (name == "none" || name.empty()) return FilterType::None;
    throw std::runtime_error("unsupported filterType: '" + name + "'");
}

const char* filterTypeName(FilterType type)
{
    switch (type) {
        case FilterType::NDiff: return "ndiff";
        case FilterType::SgDiff: return "sgdiff";
        case FilterType::Bandpass: return "bandpass";
        case FilterType::Fir1: return "fir1";
        case FilterType::User: return "user";
        case FilterType::None: return "none";
    }
    return "?";
}

CarMode carModeFromString(const std::string& name)
{
    if (name == "mean") return CarMode::Mean;
    if (name == "median") return CarMode::Median;
    if (name == "none" || name.empty()) return CarMode::None;
    // 'whiten' and friends are not implemented yet; fail loudly rather than
    // silently doing something different from MATLAB.
    throw std::runtime_error("unsupported CARMode: '" + name + "'");
}

void ndiffFilter(std::vector<float>& samples, int nSamples, int nSites, int order)
{
    // MATLAB returns the input untouched for any order outside 1..3.
    if (order < 1 || order > 3) return;
    if (nSamples <= 0 || nSites <= 0) return;

    // The stencils reach up to 3 samples either side; anything shorter has no
    // interior to write and would underflow the loop bounds below.
    if (nSamples < 8) {
        std::fill(samples.begin(), samples.end(), 0.0f);
        return;
    }

    float* base = samples.data();

    parallelFor(static_cast<std::size_t>(nSites), [&](std::size_t siteIndex) {
        float* y = base + siteIndex * static_cast<std::size_t>(nSamples);

        // MATLAB evaluates the whole right-hand side against the original array
        // before assigning, so filtering in place without a copy would feed
        // already-written values back into later taps.
        std::vector<float> x(y, y + nSamples);

        if (order == 1) {
            // y[i] = x[i+1] - x[i]; last sample zeroed.
            for (int i = 0; i <= nSamples - 2; ++i) y[i] = x[i + 1] - x[i];
            y[nSamples - 1] = 0.0f;
        } else if (order == 2) {
            // y[i] = 2*(x[i+1] - x[i]) + (x[i+2] - x[i-1]), for i in 1..N-3.
            for (int i = 1; i <= nSamples - 3; ++i)
                y[i] = 2.0f * (x[i + 1] - x[i]) + (x[i + 2] - x[i - 1]);
            y[0] = 0.0f;
            y[nSamples - 2] = 0.0f;
            y[nSamples - 1] = 0.0f;
        } else {
            // y[i] = 3*(x[i+1]-x[i]) + 2*(x[i+2]-x[i-1]) + (x[i+3]-x[i-2]),
            // for i in 2..N-4.
            for (int i = 2; i <= nSamples - 4; ++i)
                y[i] = 3.0f * (x[i + 1] - x[i]) + 2.0f * (x[i + 2] - x[i - 1]) +
                       (x[i + 3] - x[i - 2]);
            y[0] = 0.0f;
            y[1] = 0.0f;
            y[nSamples - 3] = 0.0f;
            y[nSamples - 2] = 0.0f;
            y[nSamples - 1] = 0.0f;
        }
    });
}

namespace {

/// Zeroes every ignored site's trace.
void blankIgnoredSites(float* base, int nSamples, int nSites,
                       const std::vector<std::uint8_t>& ignoredMask)
{
    if (ignoredMask.empty()) return;

    for (int site = 0; site < nSites; ++site) {
        if (!ignoredMask[site]) continue;
        float* row = base + static_cast<std::size_t>(site) * nSamples;
        std::fill(row, row + nSamples, 0.0f);
    }
}

}  // namespace

void applyCar(std::vector<float>& samples,
              int nSamples,
              int nSites,
              CarMode mode,
              const std::vector<std::uint8_t>& ignoredMask,
              std::vector<float>& channelMeans)
{
    channelMeans.assign(nSamples, 0.0f);
    float* base = samples.data();

    if (mode == CarMode::None) {
        // MATLAB blanks ignored sites unconditionally, even without CAR.
        blankIgnoredSites(base, nSamples, nSites, ignoredMask);
        return;
    }

    // Sites contributing to the reference: everything not in ignoreSites.
    std::vector<int> goodSites;
    goodSites.reserve(nSites);
    for (int site = 0; site < nSites; ++site)
        if (ignoredMask.empty() || !ignoredMask[site]) goodSites.push_back(site);

    if (goodSites.empty())
        throw std::runtime_error("every site is in ignoreSites; nothing to reference");

    if (mode == CarMode::Mean) {
        // Accumulate per thread over a slice of the sites, then reduce. Each
        // worker walks whole contiguous traces, so this stays sequential in
        // memory despite summing across sites.
        const unsigned workers = static_cast<unsigned>(
            std::min<std::size_t>(threadCount(), goodSites.size()));

        std::vector<std::vector<double>> partials(
            workers, std::vector<double>(nSamples, 0.0));

        parallelFor(workers, [&](std::size_t worker) {
            std::vector<double>& acc = partials[worker];
            for (std::size_t k = worker; k < goodSites.size(); k += workers) {
                const float* row =
                    base + static_cast<std::size_t>(goodSites[k]) * nSamples;
                for (int t = 0; t < nSamples; ++t) acc[t] += row[t];
            }
        });

        const double scale = 1.0 / static_cast<double>(goodSites.size());
        for (int t = 0; t < nSamples; ++t) {
            double sum = 0.0;
            for (unsigned w = 0; w < workers; ++w) sum += partials[w][t];
            channelMeans[t] = static_cast<float>(sum * scale);
        }
    } else {
        // Median needs every site's value at one time point, so it gathers
        // across traces. Parallelised over time blocks to hide the striding.
        parallelFor(static_cast<std::size_t>(nSamples), [&](std::size_t t) {
            std::vector<float> scratch(goodSites.size());
            for (std::size_t i = 0; i < goodSites.size(); ++i)
                scratch[i] =
                    base[static_cast<std::size_t>(goodSites[i]) * nSamples + t];

            // MATLAB's median averages the two middle values for even counts.
            const std::size_t mid = scratch.size() / 2;
            std::nth_element(scratch.begin(), scratch.begin() + mid, scratch.end());
            float median = scratch[mid];
            if (scratch.size() % 2 == 0 && mid > 0) {
                const float lower =
                    *std::max_element(scratch.begin(), scratch.begin() + mid);
                median = 0.5f * (lower + median);
            }
            channelMeans[t] = median;
        });
    }

    // Subtract from every site, then blank the ignored ones - the same order
    // as filtCAR.m, which zeroes after the subtraction.
    const float* means = channelMeans.data();
    parallelFor(static_cast<std::size_t>(nSites), [&](std::size_t site) {
        float* row = base + site * static_cast<std::size_t>(nSamples);
        for (int t = 0; t < nSamples; ++t) row[t] -= means[t];
    });

    blankIgnoredSites(base, nSamples, nSites, ignoredMask);
}

}  // namespace detect
}  // namespace jrc
