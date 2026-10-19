#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace jrc {
namespace detect {

// Sample buffers in the detection pipeline are SITE-MAJOR:
//
//     samples[site * nSamples + t]
//
// so each site's trace is contiguous in time. Filtering, threshold estimation
// and peak finding all walk one site along time, and this layout makes those
// sequential reads instead of strided ones. Only common-average referencing
// needs to cross sites, and it is written to suit.

/// Filter kinds from the `filterType` parameter.
enum class FilterType { NDiff, SgDiff, Bandpass, Fir1, User, None };

FilterType filterTypeFromString(const std::string& name);
const char* filterTypeName(FilterType type);

/// How the common reference is formed, from the `CARMode` parameter.
enum class CarMode { Mean, Median, None };

CarMode carModeFromString(const std::string& name);

/// Applies the differentiator filter in place, one site per worker thread.
///
/// Port of +jrclust/+filters/ndiffFilter.m. Orders outside 1..3 leave the data
/// untouched, as in MATLAB. Edge samples the stencil cannot reach are zeroed,
/// matching MATLAB, so downstream index arithmetic stays aligned.
void ndiffFilter(std::vector<float>& samples, int nSamples, int nSites, int order);

/// Subtracts the common average reference in place and zeroes ignored sites.
///
/// Port of applyCAR() in +jrclust/+filters/filtCAR.m. `channelMeans` receives
/// the per-sample reference (length nSamples), which motion rejection needs.
void applyCar(std::vector<float>& samples,
              int nSamples,
              int nSites,
              CarMode mode,
              const std::vector<std::uint8_t>& ignoredMask,
              std::vector<float>& channelMeans);

}  // namespace detect
}  // namespace jrc
