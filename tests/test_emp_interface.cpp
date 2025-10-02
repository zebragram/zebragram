#include "emp_interface.hpp"
#include "io_channel_impl.hpp"
#include "test_util.hpp"

static ZebraGRAM::BitType compute_encoding(const ZebraGRAM::BitType& bit,
                                           uint8_t val,
                                           const ZebraGRAM::BitType& Delta) {
  if (val == 0) {
    return bit;
  }
  ZebraGRAM::BitType result;
  for (uint i = 0; i < LAMBDA_BYTES; ++i) {
    result.label[i] = bit.label[i] ^ Delta.label[i];
  }
  return result;
}

static ZebraGRAM::WordType compute_encoding(const ZebraGRAM::WordType& word,
                                            uint64_t val,
                                            const ZebraGRAM::BitType& Delta) {
  ZebraGRAM::WordType result(word.bits.size());
  for (uint i = 0; i < word.bits.size(); ++i) {
    result.bits[i] = compute_encoding(word.bits[i], val & 1, Delta);
    val >>= 1;
  }
  return result;
}

struct ZebraGRAMTester {
  ZebraGRAM::BitType Delta;
  uint addr_width;
  uint word_width;
  uint64_t num_accesses;
  std::vector<uint64_t> addr_vals;
  std::vector<uint64_t> is_write_vals;
  std::vector<uint64_t> new_data_vals;
  std::vector<uint64_t> old_data_vals;
  std::vector<uint64_t> ref_mem_vals;
  std::vector<ZebraGRAM::WordType> addr_labels;
  std::vector<ZebraGRAM::BitType> is_write_labels;
  std::vector<ZebraGRAM::WordType> new_data_labels;
  std::vector<ZebraGRAM::WordType> old_data_labels;

  ZebraGRAMTester(uint addr_width = 4, uint word_width = 8,
                  uint64_t num_accesses = 256)
      : addr_width(addr_width),
        word_width(word_width),
        num_accesses(num_accesses),
        addr_vals(num_accesses),
        is_write_vals(num_accesses),
        new_data_vals(num_accesses),
        old_data_vals(num_accesses),
        ref_mem_vals(1UL << addr_width),
        addr_labels(num_accesses, ZebraGRAM::WordType(addr_width)),
        is_write_labels(num_accesses, ZebraGRAM::BitType()),
        new_data_labels(num_accesses, ZebraGRAM::WordType(word_width)),
        old_data_labels(num_accesses, ZebraGRAM::WordType(word_width)) {
    srand(time(0));
    ZebraGRAM::secure_random(Delta.label, LAMBDA_BYTES);
    Delta.label[0] |= 1;  // make sure it is odd

    ZebraGRAM::set_Delta(Delta);

    // Generate random data for testing
    for (uint64_t i = 0; i < num_accesses; ++i) {
      addr_vals[i] = rand() % (1UL << addr_width);
      is_write_vals[i] = rand() % 2;
      new_data_vals[i] = rand() % (1UL << word_width);
      old_data_vals[i] = ref_mem_vals[addr_vals[i]];
      if (is_write_vals[i]) {
        ref_mem_vals[addr_vals[i]] = new_data_vals[i];
      }
    }

    // Generate labels for testing
    for (uint64_t i = 0; i < num_accesses; ++i) {
      for (uint j = 0; j < addr_width; ++j) {
        ZebraGRAM::secure_random(addr_labels[i].bits[j].label, LAMBDA_BYTES);
      }
      ZebraGRAM::secure_random(is_write_labels[i].label, LAMBDA_BYTES);
      for (uint j = 0; j < word_width; ++j) {
        ZebraGRAM::secure_random(new_data_labels[i].bits[j].label,
                                 LAMBDA_BYTES);
      }
    }
  }

  void garble(ZebraGRAM::ChannelType channel) {
    // Garbler
    ZebraGRAM::ORAMType oram(addr_width, word_width, num_accesses, true);
    oram.initialize(channel);

    for (uint64_t i = 0; i < num_accesses; ++i) {
      ZebraGRAM::WordType addr = addr_labels[i];
      ZebraGRAM::BitType is_write = is_write_labels[i];
      ZebraGRAM::WordType new_data = new_data_labels[i];
      ZebraGRAM::WordType old_data = oram.access(addr, is_write, new_data);
      old_data_labels[i] = old_data;
    }

    std::cout << "Finish garbling" << std::endl;
  }

  void eval(ZebraGRAM::ChannelType channel) {
    ZebraGRAM::ORAMType oram(addr_width, word_width, num_accesses, false);
    oram.initialize(channel);

    for (uint64_t i = 0; i < num_accesses; ++i) {
      ZebraGRAM::WordType addr =
          compute_encoding(addr_labels[i], addr_vals[i], Delta);
      ZebraGRAM::BitType is_write =
          compute_encoding(is_write_labels[i], is_write_vals[i], Delta);
      ZebraGRAM::WordType new_data =
          compute_encoding(new_data_labels[i], new_data_vals[i], Delta);
      ZebraGRAM::WordType old_data = oram.access(addr, is_write, new_data);
      ZebraGRAM::WordType ref_old_data =
          compute_encoding(old_data_labels[i], old_data_vals[i], Delta);
      ASSERT_EQ(old_data, ref_old_data);
    }
  }
};

TEST(EmpORAM, MemIO) {
  emp::MemIO* io_channel = new emp::MemIO();
  ZebraGRAM::ChannelType channel(io_channel, ZebraGRAM::MEM_IO);
  ZebraGRAMTester tester;
  tester.garble(channel);
  tester.eval(channel);
  delete io_channel;
}

TEST(EmpORAM, NetIO) {
  int port = 42345;
  ZebraGRAMTester tester(3, 6, 64);
#pragma omp parallel sections
  {
#pragma omp section
    {
      emp::NetIO* garbler_channel = new emp::NetIO(nullptr, port);
      ZebraGRAM::ChannelType garbler(garbler_channel, ZebraGRAM::NET_IO);
      tester.garble(garbler);
      delete garbler_channel;
    }
#pragma omp section
    {
      emp::NetIO* evaluator_channel = new emp::NetIO("127.0.0.1", port);
      ZebraGRAM::ChannelType evaluator(evaluator_channel, ZebraGRAM::NET_IO);
      tester.eval(evaluator);
      delete evaluator_channel;
    }
  }
}

TEST(EmpORAM, HighSpeedNetIO) {
  int port = 42345;
  ZebraGRAMTester tester(3, 6, 64);
#pragma omp parallel sections
  {
#pragma omp section
    {
      emp::HighSpeedNetIO* garbler_channel =
          new emp::HighSpeedNetIO(nullptr, port, port + 1);
      ZebraGRAM::ChannelType garbler(garbler_channel,
                                     ZebraGRAM::HIGH_SPEED_NET_IO);
      tester.garble(garbler);
      delete garbler_channel;
    }
#pragma omp section
    {
      emp::HighSpeedNetIO* evaluator_channel =
          new emp::HighSpeedNetIO("127.0.0.1", port, port + 1);
      ZebraGRAM::ChannelType evaluator(evaluator_channel,
                                       ZebraGRAM::HIGH_SPEED_NET_IO);
      tester.eval(evaluator);
      delete evaluator_channel;
    }
  }
}
