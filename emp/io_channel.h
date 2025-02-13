#ifndef EMP_IO_CHANNEL_H__
#define EMP_IO_CHANNEL_H__

#include <memory>
const static int NETWORK_BUFFER_SIZE2 = 1024 * 32;
const static int NETWORK_BUFFER_SIZE = 1024 * 1024;
const static int FILE_BUFFER_SIZE = 1024 * 16;
const static int CHECK_BUFFER_SIZE = 1024 * 8;

namespace emp {
template <typename T>
class IOChannel {
 public:
  uint64_t counter = 0;
  void send_data(const void *data, size_t nbyte) {
    counter += nbyte;
    derived().send_data_internal(data, nbyte);
  }

  void recv_data(void *data, size_t nbyte) {
    derived().recv_data_internal(data, nbyte);
  }

 private:
  T &derived() { return *static_cast<T *>(this); }
};
}  // namespace emp
#endif
