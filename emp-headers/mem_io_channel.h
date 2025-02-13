#ifndef EMP_MEM_IO_CHANNEL
#define EMP_MEM_IO_CHANNEL

#include <cstring>
#include <iostream>

#include "io_channel.h"

namespace emp {
class MemIO : public IOChannel<MemIO> {
 public:
  char* buffer = nullptr;
  int64_t size = 0;
  int64_t read_pos = 0;
  int64_t cap;

  MemIO(int64_t _cap = 1024 * 1024) {
    this->cap = _cap;
    buffer = new char[cap];
    size = 0;
  }
  ~MemIO() {
    if (buffer != nullptr) delete[] buffer;
  }

  void clear() { size = 0; }
  void send_data_internal(const void* data, int64_t len) {
    if (size + len >= cap) {
      char* new_buffer = new char[2 * (cap + len)];
      memcpy(new_buffer, buffer, size);
      delete[] buffer;
      buffer = new_buffer;
      cap = 2 * (cap + len);
    }
    memcpy(buffer + size, data, len);
    size += len;
  }

  void recv_data_internal(void* data, int64_t len) {
    if (read_pos + len <= size) {
      memcpy(data, buffer + read_pos, len);
      read_pos += len;
    } else {
      std::cerr << "Reading out of bound" << std::endl;
      exit(1);
    }
  }
};
}  // namespace emp
#endif  // MEM_IO_CHANNEL_H__
