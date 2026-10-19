#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace jrc {

/// Eigendecomposition of a small real symmetric matrix, by cyclic Jacobi.
///
/// Replaces MATLAB's eig() in the PCA step. The matrices involved are the
/// covariance of a spike window - 21x21 for a 1 ms window at 20 kHz - so an
/// O(n^3)-per-sweep method with excellent accuracy is the right trade; there is
/// no reason to pull in a full LAPACK dependency for this.
///
/// @param matrix  n x n, row-major, symmetric. Only the lower triangle is read.
/// @param n       matrix dimension.
/// @param values  receives the n eigenvalues, DESCENDING.
/// @param vectors receives the eigenvectors column-major (vectors[i*n + k] is
///                element i of eigenvector k), ordered to match `values`.
inline void symmetricEigen(const std::vector<double>& matrix, int n,
                           std::vector<double>& values,
                           std::vector<double>& vectors)
{
    const std::size_t nn = static_cast<std::size_t>(n) * n;

    std::vector<double> a(matrix.begin(), matrix.begin() + nn);
    std::vector<double> v(nn, 0.0);
    for (int i = 0; i < n; ++i) v[static_cast<std::size_t>(i) * n + i] = 1.0;

    const auto at = [&](int r, int c) -> double& {
        return a[static_cast<std::size_t>(r) * n + c];
    };
    const auto vt = [&](int r, int c) -> double& {
        return v[static_cast<std::size_t>(r) * n + c];
    };

    // Symmetrise defensively: the caller builds this from an outer product, so
    // the halves can differ in the last bit.
    for (int r = 0; r < n; ++r)
        for (int c = r + 1; c < n; ++c) at(r, c) = at(c, r);

    constexpr int kMaxSweeps = 100;

    for (int sweep = 0; sweep < kMaxSweeps; ++sweep) {
        // Convergence: total magnitude of the off-diagonal entries.
        double offDiagonal = 0.0;
        for (int r = 0; r < n; ++r)
            for (int c = r + 1; c < n; ++c) offDiagonal += std::fabs(at(r, c));

        if (offDiagonal < 1e-14) break;

        for (int p = 0; p < n - 1; ++p) {
            for (int q = p + 1; q < n; ++q) {
                const double apq = at(p, q);
                if (std::fabs(apq) < 1e-300) continue;

                // Rotation angle that zeroes a(p,q).
                const double theta = (at(q, q) - at(p, p)) / (2.0 * apq);
                const double t = (theta >= 0.0 ? 1.0 : -1.0) /
                                 (std::fabs(theta) + std::sqrt(theta * theta + 1.0));
                const double c = 1.0 / std::sqrt(t * t + 1.0);
                const double s = t * c;

                for (int k = 0; k < n; ++k) {
                    const double akp = at(k, p);
                    const double akq = at(k, q);
                    at(k, p) = c * akp - s * akq;
                    at(k, q) = s * akp + c * akq;
                }
                for (int k = 0; k < n; ++k) {
                    const double apk = at(p, k);
                    const double aqk = at(q, k);
                    at(p, k) = c * apk - s * aqk;
                    at(q, k) = s * apk + c * aqk;
                }
                for (int k = 0; k < n; ++k) {
                    const double vkp = vt(k, p);
                    const double vkq = vt(k, q);
                    vt(k, p) = c * vkp - s * vkq;
                    vt(k, q) = s * vkp + c * vkq;
                }
            }
        }
    }

    // Sort descending by eigenvalue; MATLAB's eig returns ascending and
    // spikePCA.m flips it.
    std::vector<int> order(n);
    for (int i = 0; i < n; ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&](int x, int y) {
        return a[static_cast<std::size_t>(x) * n + x] >
               a[static_cast<std::size_t>(y) * n + y];
    });

    values.resize(n);
    vectors.assign(nn, 0.0);
    for (int k = 0; k < n; ++k) {
        const int src = order[k];
        values[k] = a[static_cast<std::size_t>(src) * n + src];
        for (int i = 0; i < n; ++i)
            vectors[static_cast<std::size_t>(i) * n + k] =
                v[static_cast<std::size_t>(i) * n + src];
    }
}

}  // namespace jrc
