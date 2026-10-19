#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace jrc {

/// MATLAB's quantile() for a vector.
///
/// MATLAB places the sorted samples at midpoints ((i - 0.5)/n) and linearly
/// interpolates between them, clamping to the extremes outside that range.
/// This differs from the more common (i-1)/(n-1) convention, and the
/// difference matters: it sets the distance cutoff that drives rho, and the
/// position bands that drive drift-aware averaging.
///
/// `values` is sorted in place.
inline float matlabQuantile(std::vector<float>& values, double p)
{
    if (values.empty()) return std::numeric_limits<float>::quiet_NaN();

    std::sort(values.begin(), values.end());
    const std::size_t n = values.size();
    if (n == 1) return values[0];

    const double pos = p * n - 0.5;  // 0-based index
    if (pos <= 0.0) return values.front();
    if (pos >= static_cast<double>(n - 1)) return values.back();

    const std::size_t lo = static_cast<std::size_t>(std::floor(pos));
    const double frac = pos - static_cast<double>(lo);
    return static_cast<float>(values[lo] * (1.0 - frac) + values[lo + 1] * frac);
}

}  // namespace jrc
