#include "core/Probe.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>

namespace jrc {

void Probe::load(const Config& config)
{
    const auto siteMap = config.getNumbers("siteMap");
    const auto siteLoc = config.getMatrix("siteLoc", 2);
    const auto shankMap = config.getNumbers("shankMap");

    if (siteMap.empty())
        throw std::runtime_error("siteMap is missing or empty in the .prm");
    if (siteLoc.empty())
        throw std::runtime_error("siteLoc is missing or empty in the .prm");

    nSites_ = static_cast<int>(siteMap.size());

    if (siteLoc.size() / 2 != siteMap.size()) {
        throw std::runtime_error("siteLoc has " + std::to_string(siteLoc.size() / 2) +
                                 " rows but siteMap has " +
                                 std::to_string(siteMap.size()) + " entries");
    }
    if (!shankMap.empty() && shankMap.size() != siteMap.size()) {
        throw std::runtime_error("shankMap has " + std::to_string(shankMap.size()) +
                                 " entries but siteMap has " +
                                 std::to_string(siteMap.size()));
    }

    siteMap_.resize(nSites_);
    siteX_.resize(nSites_);
    siteY_.resize(nSites_);

    const int nChans = config.nChans();

    for (int i = 0; i < nSites_; ++i) {
        // .prm channel numbers are 1-based; store 0-based.
        const int channel = static_cast<int>(siteMap[i]) - 1;
        if (channel < 0 || (nChans > 0 && channel >= nChans)) {
            throw std::runtime_error("siteMap entry " + std::to_string(i + 1) +
                                     " maps to channel " +
                                     std::to_string(channel + 1) +
                                     ", outside 1.." + std::to_string(nChans));
        }
        siteMap_[i] = channel;
        siteX_[i] = siteLoc[2 * static_cast<std::size_t>(i)];
        siteY_[i] = siteLoc[2 * static_cast<std::size_t>(i) + 1];
    }

    shankMap_.clear();
    if (!shankMap.empty()) {
        shankMap_.resize(nSites_);
        for (int i = 0; i < nSites_; ++i)
            shankMap_[i] = static_cast<int>(shankMap[i]);
    }

    ignored_.assign(nSites_, 0);
    for (double value : config.getNumbers("ignoreSites")) {
        const int site = static_cast<int>(value) - 1;  // 1-based in the .prm
        if (site >= 0 && site < nSites_) ignored_[site] = 1;
    }

    neighborCache_.clear();
}

int Probe::nIgnored() const
{
    return static_cast<int>(std::count(ignored_.begin(), ignored_.end(), 1));
}

std::vector<int> Probe::nearbySites(int site, double radius) const
{
    std::vector<int> result;
    const double r2 = radius * radius;

    const double xi = siteX_[site];
    const double yi = siteY_[site];

    for (int j = 0; j < nSites_; ++j) {
        const double dx = siteX_[j] - xi;
        const double dy = siteY_[j] - yi;
        if (dx * dx + dy * dy <= r2) result.push_back(j);
    }

    return result;
}

void Probe::buildNeighborCache(double radius)
{
    neighborCache_.assign(nSites_, {});
    for (int site = 0; site < nSites_; ++site)
        neighborCache_[site] = nearbySites(site, radius);
}

void Probe::buildSiteNeighbors(int count)
{
    nSiteNeighbors_ = std::min(count, nSites_);
    siteNeighbors_.assign(static_cast<std::size_t>(nSites_) * nSiteNeighbors_, 0);

    const bool multiShank = !shankMap_.empty();
    constexpr double kFar = std::numeric_limits<double>::infinity();

    std::vector<std::pair<double, int>> ranked(nSites_);

    for (int site = 0; site < nSites_; ++site) {
        const double xi = siteX_[site];
        const double yi = siteY_[site];

        for (int j = 0; j < nSites_; ++j) {
            double distance;

            // findSiteNeighbors.m sets ignored sites and other shanks to inf so
            // they sort last rather than being removed outright - the table is
            // a fixed size and must always be filled.
            if (ignored_[j] || (multiShank && shankMap_[j] != shankMap_[site])) {
                distance = kFar;
            } else {
                const double dx = siteX_[j] - xi;
                const double dy = siteY_[j] - yi;
                distance = std::sqrt(dx * dx + dy * dy);
            }
            ranked[j] = {distance, j};
        }

        // Ties broken by site index, so the table is deterministic.
        std::partial_sort(ranked.begin(), ranked.begin() + nSiteNeighbors_,
                          ranked.end());

        int* row = siteNeighbors_.data() + static_cast<std::size_t>(site) * nSiteNeighbors_;
        for (int k = 0; k < nSiteNeighbors_; ++k) row[k] = ranked[k].second;
    }
}

}  // namespace jrc
