#include "core/MemoryMappedFile.h"

#include <stdexcept>
#include <utility>

#ifdef _WIN32
#include <windows.h>
#else
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

namespace jrc {

MemoryMappedFile::~MemoryMappedFile() { close(); }

MemoryMappedFile::MemoryMappedFile(MemoryMappedFile&& other) noexcept
    : data_(other.data_), size_(other.size_), path_(std::move(other.path_))
#ifdef _WIN32
      ,
      fileHandle_(other.fileHandle_), mapHandle_(other.mapHandle_)
#else
      ,
      fd_(other.fd_)
#endif
{
    other.data_ = nullptr;
    other.size_ = 0;
#ifdef _WIN32
    other.fileHandle_ = nullptr;
    other.mapHandle_ = nullptr;
#else
    other.fd_ = -1;
#endif
}

MemoryMappedFile& MemoryMappedFile::operator=(MemoryMappedFile&& other) noexcept
{
    if (this != &other) {
        close();
        data_ = other.data_;
        size_ = other.size_;
        path_ = std::move(other.path_);
#ifdef _WIN32
        fileHandle_ = other.fileHandle_;
        mapHandle_ = other.mapHandle_;
        other.fileHandle_ = nullptr;
        other.mapHandle_ = nullptr;
#else
        fd_ = other.fd_;
        other.fd_ = -1;
#endif
        other.data_ = nullptr;
        other.size_ = 0;
    }
    return *this;
}

#ifdef _WIN32

void MemoryMappedFile::open(const std::string& path)
{
    close();

    HANDLE file = ::CreateFileA(path.c_str(), GENERIC_READ, FILE_SHARE_READ, nullptr,
                                OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, nullptr);
    if (file == INVALID_HANDLE_VALUE)
        throw std::runtime_error("cannot open file: " + path);

    LARGE_INTEGER fileSize{};
    if (!::GetFileSizeEx(file, &fileSize)) {
        ::CloseHandle(file);
        throw std::runtime_error("cannot stat file: " + path);
    }
    if (fileSize.QuadPart == 0) {
        ::CloseHandle(file);
        throw std::runtime_error("file is empty: " + path);
    }

    HANDLE mapping = ::CreateFileMappingA(file, nullptr, PAGE_READONLY, 0, 0, nullptr);
    if (!mapping) {
        ::CloseHandle(file);
        throw std::runtime_error("cannot create file mapping: " + path);
    }

    void* view = ::MapViewOfFile(mapping, FILE_MAP_READ, 0, 0, 0);
    if (!view) {
        ::CloseHandle(mapping);
        ::CloseHandle(file);
        throw std::runtime_error("cannot map view of file: " + path);
    }

    fileHandle_ = file;
    mapHandle_ = mapping;
    data_ = static_cast<std::uint8_t*>(view);
    size_ = static_cast<std::size_t>(fileSize.QuadPart);
    path_ = path;
}

void MemoryMappedFile::close()
{
    if (data_) ::UnmapViewOfFile(data_);
    if (mapHandle_) ::CloseHandle(mapHandle_);
    if (fileHandle_ && fileHandle_ != INVALID_HANDLE_VALUE) ::CloseHandle(fileHandle_);

    data_ = nullptr;
    size_ = 0;
    mapHandle_ = nullptr;
    fileHandle_ = nullptr;
    path_.clear();
}

#else  // POSIX

void MemoryMappedFile::open(const std::string& path)
{
    close();

    const int fd = ::open(path.c_str(), O_RDONLY);
    if (fd < 0) throw std::runtime_error("cannot open file: " + path);

    struct stat st {};
    if (::fstat(fd, &st) != 0) {
        ::close(fd);
        throw std::runtime_error("cannot stat file: " + path);
    }
    if (st.st_size == 0) {
        ::close(fd);
        throw std::runtime_error("file is empty: " + path);
    }

    void* view = ::mmap(nullptr, static_cast<std::size_t>(st.st_size), PROT_READ,
                        MAP_PRIVATE, fd, 0);
    if (view == MAP_FAILED) {
        ::close(fd);
        throw std::runtime_error("cannot mmap file: " + path);
    }

    fd_ = fd;
    data_ = static_cast<std::uint8_t*>(view);
    size_ = static_cast<std::size_t>(st.st_size);
    path_ = path;
}

void MemoryMappedFile::close()
{
    if (data_) ::munmap(data_, size_);
    if (fd_ >= 0) ::close(fd_);

    data_ = nullptr;
    size_ = 0;
    fd_ = -1;
    path_.clear();
}

#endif

}  // namespace jrc
