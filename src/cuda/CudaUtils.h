#pragma once

#include <cuda_runtime.h>

#include <stdexcept>
#include <string>
#include <vector>

namespace jrc {
namespace cuda {

/// Throws std::runtime_error carrying the CUDA error name and the call site.
inline void check(cudaError_t err, const char* what, const char* file, int line)
{
    if (err != cudaSuccess) {
        throw std::runtime_error(std::string("CUDA error in ") + what + " at " +
                                 file + ":" + std::to_string(line) + " - " +
                                 cudaGetErrorString(err));
    }
}

#define JRC_CUDA_CHECK(call) ::jrc::cuda::check((call), #call, __FILE__, __LINE__)

/// RAII device buffer. Frees on scope exit so a throwing launch cannot leak.
template <typename T>
class DeviceBuffer {
public:
    DeviceBuffer() = default;

    explicit DeviceBuffer(std::size_t count) { allocate(count); }

    ~DeviceBuffer() { release(); }

    DeviceBuffer(const DeviceBuffer&) = delete;
    DeviceBuffer& operator=(const DeviceBuffer&) = delete;

    DeviceBuffer(DeviceBuffer&& o) noexcept : ptr_(o.ptr_), count_(o.count_)
    {
        o.ptr_ = nullptr;
        o.count_ = 0;
    }

    DeviceBuffer& operator=(DeviceBuffer&& o) noexcept
    {
        if (this != &o) {
            release();
            ptr_ = o.ptr_;
            count_ = o.count_;
            o.ptr_ = nullptr;
            o.count_ = 0;
        }
        return *this;
    }

    void allocate(std::size_t count)
    {
        release();
        if (count == 0) return;
        JRC_CUDA_CHECK(cudaMalloc(&ptr_, count * sizeof(T)));
        count_ = count;
    }

    void release()
    {
        if (ptr_) {
            cudaFree(ptr_);  // destructor path: cannot throw
            ptr_ = nullptr;
        }
        count_ = 0;
    }

    void upload(const T* host, std::size_t count)
    {
        JRC_CUDA_CHECK(cudaMemcpy(ptr_, host, count * sizeof(T), cudaMemcpyHostToDevice));
    }

    void download(T* host, std::size_t count) const
    {
        JRC_CUDA_CHECK(cudaMemcpy(host, ptr_, count * sizeof(T), cudaMemcpyDeviceToHost));
    }

    void zero()
    {
        if (ptr_) JRC_CUDA_CHECK(cudaMemset(ptr_, 0, count_ * sizeof(T)));
    }

    T* get() { return ptr_; }
    const T* get() const { return ptr_; }
    std::size_t count() const { return count_; }

private:
    T* ptr_ = nullptr;
    std::size_t count_ = 0;
};

/// Summary of a CUDA device, for `jrclust gpuinfo`.
struct DeviceInfo {
    int index = 0;
    std::string name;
    int major = 0;
    int minor = 0;
    std::size_t totalMemory = 0;
    int multiProcessorCount = 0;
    int maxThreadsPerBlock = 0;
    std::size_t sharedMemPerBlock = 0;
};

/// Returns every visible CUDA device. Empty if none (no GPU / no driver).
std::vector<DeviceInfo> enumerateDevices();

}  // namespace cuda
}  // namespace jrc
