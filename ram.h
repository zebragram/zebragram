#ifndef EMP_RAM_H__
#define EMP_RAM_H__
#include "emp-tool/circuits/bit.h"
#include "emp-tool/circuits/integer.h"
#include "emp-tool/io/file_io_channel.h"
#include "emp-tool/io/highspeed_net_io_channel.h"
#include "emp-tool/io/mem_io_channel.h"
#include "emp-tool/io/net_io_channel.h"
#include "emp/emp_interface.hpp"

namespace emp {
class RAM {
 private:
  PicoGRAM::ORAMType oram;
  uint64_t times_accessed = 0;
  uint64_t max_times_accessed = 0;

  static PicoGRAM::BitType bit_to_bit_type(const Bit& bit) {
    PicoGRAM::BitType bit_type;
    memcpy(bit_type.label, &bit.bit, LAMBDA_BYTES);
    return bit_type;
  }

  static Bit bit_type_to_bit(const PicoGRAM::BitType& bit_type) {
    Bit bit;
    memcpy(&bit.bit, bit_type.label, LAMBDA_BYTES);
    return bit;
  }

  static PicoGRAM::WordType integer_to_word_type(const Integer& integer) {
    PicoGRAM::WordType word_type(integer.size());
    for (size_t i = 0; i < integer.size(); ++i) {
      word_type.bits[i] = bit_to_bit_type(integer.bits[i]);
    }
    return word_type;
  }

  static Integer word_type_to_integer(const PicoGRAM::WordType& word_type) {
    std::vector<Bit> bits(word_type.bits.size());
    for (size_t i = 0; i < word_type.bits.size(); ++i) {
      bits[i] = bit_type_to_bit(word_type.bits[i]);
    }
    return Integer(bits);
  }

 public:
  RAM(const RAM&) = delete;
  RAM& operator=(const RAM&) = delete;

  static void set_global_delta(block delta) {
    PicoGRAM::set_Delta(bit_to_bit_type(Bit(delta)));
  }

  RAM(int addr_width, int word_width, int num_accesses, int party)
      : oram(addr_width, word_width, num_accesses, party == ALICE) {
    Assert(PicoGRAM::Delta_set_flag);
    max_times_accessed = num_accesses;
  }

  template <typename IOChannel>
  void initialize(IOChannel* io) {
    PicoGRAM::IOType io_type;
    // enumerate possible IO types
    if (std::is_same<IOChannel, emp::MemIO>()) {
      io_type = PicoGRAM::IOType::MEM_IO;
    } else if (std::is_same<IOChannel, emp::FileIO>()) {
      io_type = PicoGRAM::IOType::FILE_IO;
    } else if (std::is_same<IOChannel, emp::HighSpeedNetIO>()) {
      io_type = PicoGRAM::IOType::HIGH_SPEED_NET_IO;
    } else if (std::is_same<IOChannel, emp::NetIO>()) {
      io_type = PicoGRAM::IOType::NET_IO;
    } else {
      throw std::runtime_error("Unsupported IO channel type");
    }

    PicoGRAM::Channel::ChannelType channel_type(io, io_type);
    oram.initialize(channel_type);
  }

  Integer access(const Integer& addr, const Bit& is_write,
                 const Integer& new_data) {
    ++times_accessed;
    if (times_accessed > max_times_accessed) {
      throw std::runtime_error("RAM access limit exceeded");
    }
    const PicoGRAM::WordType& addr_ = integer_to_word_type(addr);
    const PicoGRAM::BitType& is_write_ = bit_to_bit_type(is_write);
    const PicoGRAM::WordType& new_data_ = integer_to_word_type(new_data);
    const PicoGRAM::WordType& old_data_ =
        oram.access(addr_, is_write_, new_data_);
    return word_type_to_integer(old_data_);
  }

  ~RAM() {}
};
}  // namespace emp
#endif
