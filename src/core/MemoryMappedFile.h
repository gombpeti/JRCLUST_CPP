#pragma once

#include <cstddef>
#include <cstdint>
#include <string>

namespace jrc {

/// Read-only memory map over a file.
///
/// Raw recordings run to tens or hundreds of GB, so they are mapped rather than
/// read — the same reason the MATLAB code used memmapfile() in openRaw.m.
class MemoryMappedFile {
public:
    MemoryMappedFile() = default;
    ~MemoryMappedFile();

    MemoryMappedFile(const MemoryMappedFile&) = delete;
    MemoryMappedFile& operator=(const MemoryMappedFile&) = delete;

    MemoryMappedFile(MemoryMappedFile&& other) noexcept;
    MemoryMappedFile& operator=(MemoryMappedFile&& other) noexcept;

    /// Maps `path` read-only. Throws std::runtime_error on failure.
    void open(const std::string& path);
    void close();

    bool isOpen() const { return data_ != nullptr; }
    const std::uint8_t* data() const { return data_; }
    std::size_t size() const { return size_; }
    const std::string& path() const { return path_; }

private:
    std::uint8_t* data_ = nullptr;
    std::size_t size_ = 0;
    std::string path_;

#ifdef _WIN32
    void* fileHandle_ = nullptr;  // HANDLE
    void* mapHandle_ = nullptr;   // HANDLE
#else
    int fd_ = -1;
#endif
};

}  // namespace jrc
