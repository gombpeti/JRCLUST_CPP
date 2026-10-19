#pragma once

#include "core/Config.h"

#include <cstdint>
#include <vector>

namespace jrc {

/// Probe geometry: which channel each site lives on, where the sites sit, and
/// which shank they belong to.
///
/// Site indices are 0-based throughout. The .prm file stores them 1-based (the
/// MATLAB convention), so `siteMap` is decremented on load and `ignoreSites` is
/// converted at the same time. Anything written back out must add 1 again.
class Probe {
public:
    /// Builds the probe from siteMap / siteLoc / shankMap / ignoreSites.
    /// Throws std::runtime_error if the arrays disagree on site count.
    void load(const Config& config);

    int nSites() const { return nSites_; }

    /// Recording channel index (0-based) that site `site` is stored on.
    int channelForSite(int site) const { return siteMap_[site]; }
    const std::vector<int>& siteMap() const { return siteMap_; }

    double x(int site) const { return siteX_[site]; }
    double y(int site) const { return siteY_[site]; }

    int shank(int site) const { return shankMap_.empty() ? 1 : shankMap_[site]; }

    /// True if the site was listed in ignoreSites and must not produce spikes.
    bool isIgnored(int site) const { return ignored_[site]; }
    const std::vector<std::uint8_t>& ignoredMask() const { return ignored_; }
    int nIgnored() const;

    /// Sites within `radius` micrometres of `site`, including `site` itself.
    /// Port of jrclust.utils.findNearbySites. Sites on another shank are
    /// excluded, matching the MATLAB behaviour of comparing raw coordinates
    /// (shanks are laid out far apart in x).
    std::vector<int> nearbySites(int site, double radius) const;

    /// Precomputes nearbySites() for every site. Detection needs this for all
    /// sites at once and the O(nSites^2) scan is worth doing only once.
    void buildNeighborCache(double radius);
    const std::vector<int>& cachedNeighbors(int site) const
    {
        return neighborCache_[site];
    }
    bool hasNeighborCache() const { return !neighborCache_.empty(); }

    /// Builds the ordered nearest-neighbour table used for feature extraction.
    ///
    /// Port of findSiteNeighbors.m: for each site, the `count` closest sites in
    /// ascending distance, with the site itself first. Sites on another shank
    /// and sites in ignoreSites are pushed to the end by treating them as
    /// infinitely far, so they only appear if there is nothing better.
    void buildSiteNeighbors(int count);

    int nSiteNeighbors() const { return nSiteNeighbors_; }

    /// Ordered neighbours of `site`; entry 0 is `site` itself (unless ignored).
    const int* siteNeighbors(int site) const
    {
        return siteNeighbors_.data() + static_cast<std::size_t>(site) * nSiteNeighbors_;
    }

private:
    int nSites_ = 0;
    std::vector<int> siteMap_;   ///< site -> channel, 0-based
    std::vector<double> siteX_;
    std::vector<double> siteY_;
    std::vector<int> shankMap_;
    std::vector<std::uint8_t> ignored_;
    std::vector<std::vector<int>> neighborCache_;

    int nSiteNeighbors_ = 0;
    std::vector<int> siteNeighbors_;  ///< nSiteNeighbors_ x nSites, site-major
};

}  // namespace jrc
