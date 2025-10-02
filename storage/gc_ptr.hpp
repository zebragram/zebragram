#pragma once
#include <cstring>

#include "arith_label.hpp"
#include "ec.hpp"
#include "label.hpp"
#include "storage_util.hpp"
#include "util.hpp"

namespace ZebraGRAM {
/**
 * @brief A pointer to the global circuit data
 *
 */
struct GCPtr {
 private:
  uint64_t offset;
  int fid = -1;

 public:
  explicit GCPtr(int fid) : offset(0), fid(fid) {}

  explicit GCPtr(int fid, uint64_t offset) : offset(offset), fid(fid) {
    seek(fid, offset);
  }

  bool is_valid() const { return fid >= 0; }

  void set_offset(uint64_t offset) { this->offset = offset; }

  int get_fid() const { return fid; }

  uint64_t get_offset() const { return offset; }

  // overload operator +
  GCPtr operator+(uint64_t size) const { return GCPtr(fid, offset + size); }

  GCPtr &operator+=(uint64_t size) {
    offset += size;
    return *this;
  }

  GCPtr operator-(uint64_t size) const { return GCPtr(fid, offset - size); }

  GCPtr &operator-=(uint64_t size) {
    offset -= size;
    return *this;
  }

  void rewind_write(uint64_t rewind_size) {
    Assert(rewind_size <= offset);
    offset -= rewind_size;
    ZebraGRAM::rewind_write(fid, rewind_size);
  }

  void write_data(const void *data, uint64_t size) {
    // the garbler should always write data sequentially
    // Assert_eq((int64_t)offset, lseek64(fid, 0, SEEK_CUR));
    // std::cout << "Write " << size << " bytes of data ";
    // for (uint64_t i = 0; i < size; ++i) {
    //   std::cout << (uint16_t)data[i] << " ";
    // }
    // std::cout << std::endl;
    ssize_t bytes_left = size;
    while (bytes_left > 0) {
      ssize_t bytes_written = ZebraGRAM::write(fid, data, bytes_left);
      if (bytes_written < 0) {
        throw std::runtime_error("Failed to write data to file");
      }
      bytes_left -= bytes_written;
      data = (uint8_t *)data + bytes_written;
    }
    offset += size;
  }

  void read_data(void *data, uint64_t size) {
    peak_data(data, 0, size);
    offset += size;
    // std::cout << "Read " << size << " bytes of data ";
    // for (uint64_t i = 0; i < size; ++i) {
    //   std::cout << (uint16_t)data[i] << " ";
    // }
    // std::cout << std::endl;
  }

  void peak_data(void *data, uint64_t peak_offset, uint64_t size) const {
    // increase the offset by peak_offset
    uint64_t read_offset = offset + peak_offset;
    ssize_t bytes_left = size;
    while (bytes_left > 0) {
      ssize_t bytes_read = ZebraGRAM::pread(fid, data, bytes_left, read_offset);
      Assert(bytes_read >= 0);
      if (bytes_read < 0) {
        throw std::runtime_error("Failed to read data from file");
      }
      bytes_left -= bytes_read;
      read_offset += bytes_read;
      data = (uint8_t *)data + bytes_read;
    }
  }

  void skip_data(uint64_t size) {
#ifdef MEASURE_STACK_COST
    if (measure_stack_flag) {
      global_stack_cost += size;
    }
#endif
    offset += size;
  }

  void read_label(Label &label) { read_data(label.get_ptr(), sizeof(Label)); }

  void skip_label() { skip_data(sizeof(Label)); }

  void skip_label(uint64_t count) { skip_data(count * sizeof(Label)); }

  void write_label(const Label &label) {
    write_data(label.get_ptr(), sizeof(Label));
  }

  void read_ec_point(ECPoint &point) {
    read_data(point.bytes, ECPoint::byte_length);
    point.is_temp_point_fresh = false;
  }

  void write_ec_point(const ECPoint &point) {
    write_data(point.bytes, ECPoint::byte_length);
  }

  void skip_ec_point() { skip_data(ECPoint::byte_length); }

  void skip_ec_point(uint64_t count) {
    skip_data(count * ECPoint::byte_length);
  }

  void read_big_int(BigInt &big_int) {
    uint8_t buffer[BigInt::byte_length];
    read_data(buffer, BigInt::byte_length);
    big_int.from_bytes(buffer);
  }

  void peak_big_int(BigInt &big_int, uint64_t index) const {
    uint8_t buffer[BigInt::byte_length];
    peak_data(buffer, index * BigInt::byte_length, BigInt::byte_length);
    big_int.from_bytes(buffer);
  }

  void write_big_int(BigInt &big_int) {
    uint8_t buffer[BigInt::byte_length];
    big_int.to_bytes(buffer);
    write_data(buffer, BigInt::byte_length);
  }

  void skip_big_int() { skip_data(BigInt::byte_length); }

  void skip_big_int(uint64_t count) { skip_data(count * BigInt::byte_length); }

  void read_mac(MAC &mac) { read_data(mac.bits, sizeof(MAC)); }

  void write_mac(const MAC &mac) { write_data(mac.bits, sizeof(MAC)); }

  void skip_mac() { skip_data(sizeof(MAC)); }

  void skip_mac(uint64_t count) { skip_data(count * sizeof(MAC)); }

  void read_bit(uint8_t &bit) {
    read_data(&bit, 1);
    Assert(bit <= 1);
  }

  void write_bit(uint8_t bit) {
    Assert(bit <= 1);
    write_data(&bit, 1);
  }

  void skip_bit() { skip_data(1); }

  template <typename type>
  void write_default(const type &data) {
    static_assert(std::is_trivially_copyable<type>::value);
    write_data((uint8_t *)&data, sizeof(type));
  }

  template <typename type>
  void read_default(type &data) {
    static_assert(std::is_trivially_copyable<type>::value);
    read_data(&data, sizeof(type));
  }

  template <typename type>
  void peak_default(type &data, uint64_t index) const {
    static_assert(std::is_trivially_copyable<type>::value);
    peak_data(&data, index * sizeof(type), sizeof(type));
  }

  template <typename type>
  void peak_default(std::vector<type> &data, uint64_t index) const {
    static_assert(std::is_trivially_copyable<type>::value);
    peak_data(&data[0], index * sizeof(type), data.size() * sizeof(type));
  }

  template <typename type>
  void skip_default() {
    static_assert(std::is_trivially_copyable<type>::value);
    skip_data(sizeof(type));
  }

  template <typename type>
  void skip_default(uint64_t count) {
    static_assert(std::is_trivially_copyable<type>::value);
    skip_data(count * sizeof(type));
  }

  void read_arith_label(ArithLabel &arith_label) {
    uint8_t buffer[ArithLabel::byte_length];
    read_data(buffer, ArithLabel::byte_length);
    arith_label.from_bytes(buffer);
  }

  void write_arith_label(const ArithLabel &arith_label) {
    uint8_t buffer[ArithLabel::byte_length];
    arith_label.to_bytes(buffer);
    // // print out buffer in hex
    // std::cout << "Binary of arith label to write: ";
    // for (uint i = 0; i < ArithLabel::byte_length; i++) {
    //   std::cout << std::hex << (int)buffer[i];
    // }
    // std::cout << std::dec << std::endl;
    write_data(buffer, ArithLabel::byte_length);
  }

  void skip_arith_label() { skip_data(ArithLabel::byte_length); }

  void skip_arith_label(uint64_t count) {
    skip_data(count * ArithLabel::byte_length);
  }

  void peak_arith_label(ArithLabel &arith_label, uint64_t index) const {
    uint8_t buffer[ArithLabel::byte_length];
    peak_data(buffer, index * ArithLabel::byte_length, ArithLabel::byte_length);
    // std::cout << "Binary of arith label to peak at index " << index << ": ";
    // for (uint i = 0; i < ArithLabel::byte_length; i++) {
    //   std::cout << std::hex << (int)buffer[i];
    // }
    // std::cout << std::dec << std::endl;
    arith_label.from_bytes(buffer);
  }

  friend std::ostream &operator<<(std::ostream &os, const GCPtr &ptr) {
    // get file name from file descriptor
    char path[PATH_MAX];
    char proc_path[PATH_MAX];
    snprintf(proc_path, PATH_MAX, "/proc/self/fd/%d", ptr.fid);
    ssize_t len = readlink(proc_path, path, PATH_MAX);
    if (len == -1) {
      throw std::runtime_error("Failed to get file path from file descriptor");
    }
    path[len] = '\0';
    os << path << " offset: " << ptr.get_offset();
    return os;
  }
};
}  // namespace ZebraGRAM