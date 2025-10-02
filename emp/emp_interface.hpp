#pragma once
/**
 * @brief An interface for the ZebraGRAM library compatible with the EMP
 * framework
 *
 */
#include <cstdint>
#include <vector>

#include "../global.hpp"
#include "channel_util.hpp"

namespace ZebraGRAM {
/**
 * @brief A minimal abstraction of a bit garbled with Free-XOR
 *
 */
struct BitType {
  uint8_t label[LAMBDA_BYTES];

  bool operator==(const BitType& other) const {
    for (uint32_t i = 0; i < LAMBDA_BYTES; ++i) {
      if (label[i] != other.label[i]) {
        return false;
      }
    }
    return true;
  }

  bool operator!=(const BitType& other) const { return !(*this == other); }

  friend std::ostream& operator<<(std::ostream& os, const BitType& bit) {
    for (uint32_t i = 0; i < LAMBDA_BYTES; ++i) {
      os << (uint16_t)bit.label[i] << " ";
    }
    return os;
  }
};

/**
 * @brief A minimal abstraction of a word garbled with Free-XOR
 *
 */
struct WordType {
  std::vector<BitType> bits;
  WordType() = default;
  explicit WordType(uint32_t width) : bits(width) {}
  // overload == operator
  bool operator==(const WordType& other) const {
    if (bits.size() != other.bits.size()) {
      return false;
    }
    for (uint32_t i = 0; i < bits.size(); ++i) {
      if (bits[i] != other.bits[i]) {
        return false;
      }
    }
    return true;
  }

  bool operator!=(const WordType& other) const { return !(*this == other); }

  friend std::ostream& operator<<(std::ostream& os, const WordType& word) {
    os << "\n";
    for (const auto& bit : word.bits) {
      os << bit << "\n";
    }
    return os;
  }
};

using ChannelType = Channel::ChannelType;

extern bool Delta_set_flag;

/**
 * @brief Set the global Delta value for the ZebraGRAM library
 *
 * @param Delta the XOR of the zero label and the one label of the same wire,
 * Delta.bit[0] must be odd
 */
void set_Delta(const BitType& Delta);

/**
 * @brief A class representing an Oblivious RAM (ORAM) instance
 *
 */
struct ORAMType {
 private:
  bool is_garbler;
  ChannelType io_channel;
  void* internal_oram;
  uint32_t addr_width;
  uint32_t word_width;
  uint64_t num_accesses;
  uint64_t acccess_counter = 0;
  std::vector<WordType> mocked_input_addrs;
  std::vector<BitType> mocked_input_is_writes;
  std::vector<WordType> mocked_input_new_data;
  std::vector<WordType> mocked_output_old_data;

  std::vector<uint64_t> measure_gadget_gc() const;
  void garble_oram();

 public:
  /**
   * @brief Construct a new ORAMType object
   *
   * @param address_width the width of the address, the memory space is
   * 2^address_width
   * @param word_width the width of the data word
   * @param num_accesses the number of accesses to the ORAM
   * @param is_garbler whether the instance is run by the garbler
   */
  ORAMType(uint32_t address_width, uint32_t word_width, uint64_t num_accesses,
           bool is_garbler);

  // delete copy constructor
  ORAMType(const ORAMType&) = delete;
  ORAMType& operator=(const ORAMType&) = delete;

  /**
   * @brief Initialize an empty ORAM. The garbler will garble the ORAM tree and
   * send the GC to the evaluator. The evaluator will receive the GC and store
   * it locally for later accessees.
   *
   * @param channel the communication channel between the garbler and the
   * evaluator
   *
   */
  void initialize(ChannelType channel);

  /**
   * @brief Perform a read or write access to the ORAM
   *
   * @param addr the address of the data
   * @param is_write whether the operation is a write
   * @param new_data the new data to write, ignored for read
   * @return WordType the old data at addr, for both read and write
   */
  WordType access(const WordType& addr, const BitType& is_write,
                  const WordType& new_data);
  ~ORAMType();
};

}  // namespace ZebraGRAM