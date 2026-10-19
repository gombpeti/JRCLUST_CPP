#pragma once

#include <algorithm>
#include <atomic>
#include <cstdlib>
#include <exception>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>

namespace jrc {

/// Worker threads to use. Honours JRC_THREADS when set, else hardware_concurrency.
inline unsigned threadCount()
{
    static const unsigned count = [] {
        if (const char* env = std::getenv("JRC_THREADS")) {
            const int requested = std::atoi(env);
            if (requested > 0) return static_cast<unsigned>(requested);
        }
        const unsigned hw = std::thread::hardware_concurrency();
        return hw > 0 ? hw : 1u;
    }();
    return count;
}

/// Runs `body(i)` for i in [0, count), spread across worker threads.
///
/// Work is claimed atomically one index at a time, which keeps threads busy
/// when per-index cost varies - as it does across sites, where a noisy site
/// yields far more candidate peaks than a quiet one.
///
/// The first exception thrown by any worker is rethrown on the calling thread.
template <typename Body>
void parallelFor(std::size_t count, const Body& body)
{
    if (count == 0) return;

    const unsigned workers = static_cast<unsigned>(
        std::min<std::size_t>(threadCount(), count));

    if (workers <= 1) {
        for (std::size_t i = 0; i < count; ++i) body(i);
        return;
    }

    std::atomic<std::size_t> next{0};
    std::exception_ptr firstError;
    std::mutex errorMutex;

    std::vector<std::thread> threads;
    threads.reserve(workers);

    for (unsigned w = 0; w < workers; ++w) {
        threads.emplace_back([&] {
            try {
                for (;;) {
                    const std::size_t i = next.fetch_add(1, std::memory_order_relaxed);
                    if (i >= count) break;
                    body(i);
                }
            } catch (...) {
                std::lock_guard<std::mutex> lock(errorMutex);
                if (!firstError) firstError = std::current_exception();
            }
        });
    }

    for (auto& thread : threads) thread.join();

    if (firstError) std::rethrow_exception(firstError);
}

}  // namespace jrc
