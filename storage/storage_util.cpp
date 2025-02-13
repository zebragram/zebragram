#include "storage_util.hpp"

#include <arpa/inet.h>
#include <fcntl.h>
#include <netinet/in.h>
#include <sys/mman.h>
#include <sys/socket.h>
#include <unistd.h>

#include <iostream>
#include <sstream>

#include "global.hpp"
#ifdef USE_EMP_CHANNEL
#include "channel_util.hpp"
#else
#include <cstring>
#endif

namespace PicoGRAM {

std::vector<std::vector<uint8_t>> shared_memory;
std::unordered_map<std::string, int> filename_to_fid;

int open_file(const char* filename, size_t size) {
  int fd;

  switch (STORAGE_TYPE) {
    case SHARED_MEMORY:
      // check if the file already exists
      if (filename_to_fid.find(filename) != filename_to_fid.end()) {
        return filename_to_fid[filename];
      }
      fd = shared_memory.size();
      // Create a file descriptor for a shared memory file
      shared_memory.push_back(std::vector<uint8_t>());
      filename_to_fid[filename] = fd;
      break;
    case MEMORY:
      // Create a file descriptor for an in-memory file
      fd = memfd_create(filename, MFD_CLOEXEC);
      break;
    case DISK:
      fd = open(filename, O_RDWR | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR);
      break;
    default:
      std::cerr << "Invalid storage type\n";
      return -1;
  }

  if (fd == -1) {
    std::cerr << "Error creating file\n";
    return fd;
  }
  if (size > 0) {
#if STORAGE_TYPE == SHARED_MEMORY
    shared_memory[fd].reserve(size);
#else
    // Resize the in-memory file
    if (ftruncate64(fd, size) == -1) {
      std::cerr << "Error resizing in-memory file\n";
      close(fd);
      return -1;
    }
#endif
  }
  return fd;
}

ssize_t write(int fid, const void* data, ssize_t size) {
#ifdef USE_EMP_CHANNEL
  if (Channel::is_channel_id(fid)) {
    return Channel::write(fid, data, size);
  }
#endif
#if STORAGE_TYPE == SHARED_MEMORY
  if ((size_t)fid >= shared_memory.size()) {
    std::cerr << "Invalid file descriptor\n";
    return -1;
  }
  shared_memory[fid].insert(shared_memory[fid].end(), (const uint8_t*)data,
                            (const uint8_t*)data + size);
  return size;
#else
  return ::write(fid, data, size);
#endif
}

// shadow pread
ssize_t pread(int fid, void* data, ssize_t size, uint64_t offset) {
#if STORAGE_TYPE == SHARED_MEMORY
  if ((size_t)fid >= shared_memory.size()) {
    std::cerr << "Invalid file descriptor\n";
    return -1;
  }
  if (offset + size > shared_memory[fid].size()) {
    std::cerr << "Invalid offset\n";
    return -1;
  }
  memcpy(data, shared_memory[fid].data() + offset, size);
  return size;
#else
  return ::pread(fid, data, size, offset);
#endif
}

int close_file(int fid) {
#ifdef USE_EMP_CHANNEL
  if (Channel::is_channel_id(fid)) {
    return Channel::close(fid);
  }
#endif
#if STORAGE_TYPE == SHARED_MEMORY
  if ((size_t)fid >= shared_memory.size()) {
    std::cerr << "Invalid file descriptor\n";
    return -1;
  }
  std::vector<uint8_t> empty;
  shared_memory[fid].swap(empty);
  return 0;
#else
  return ::close(fid);
#endif
}
}  // namespace PicoGRAM