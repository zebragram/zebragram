#pragma once

#include <sys/types.h>

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace ZebraGRAM {

extern std::vector<std::vector<uint8_t>> shared_memory;
extern std::unordered_map<std::string, int> filename_to_fid;

int open_file(const char* filename = "default.bin", size_t size = 0);

ssize_t write(int fid, const void* data, ssize_t size);

ssize_t pread(int fid, void* data, ssize_t size, uint64_t offset);

ssize_t seek(int fid, uint64_t offset);

ssize_t rewind_write(int fid, uint64_t offset);
int close_file(int fid);
}  // namespace ZebraGRAM