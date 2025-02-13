#include <algorithm>
#include <random>

#include "circuit_oram.hpp"
#include "test_util.hpp"
using namespace PicoGRAM;
struct CircuitORAMMain : Gadget {
  CircuitORAM oram;
  uint v_addr_width;
  uint position_width;
  uint word_width;
  uint64_t memory_space;
  std::vector<uint64_t> ref_mem;
  std::vector<uint64_t> new_p_addr_list;
  std::vector<uint64_t> pos_map;
  uint64_t new_p_addr_idx;
  std::function<FuncOutput(FuncInput)> main_func = [&](FuncInput) {
    uint64_t v_addr_val = rand() % memory_space;
    Word v_addr = Word::input_dbg(self, v_addr_width, v_addr_val);
    uint64_t p_addr_val = pos_map[v_addr_val];
    Word p_addr = Word::input_dbg(self, position_width, p_addr_val);
    uint8_t is_write_val = rand() % 2;
    Bit is_write = Bit::input_dbg(self, is_write_val);
    uint64_t data_val = rand() % (1UL << word_width);
    Word data = Word::input_dbg(self, word_width, data_val);
    uint64_t new_p_addr_val = new_p_addr_list[new_p_addr_idx++];
    Word new_p_addr = Word::input_g(self, position_width, new_p_addr_val);
    // std::cout << (is_write_val ? "write " + data_val : "read ") << " from "
    //           << p_addr_val << " new addr " << new_p_addr_val << std::endl;
    const Word& output_word =
        oram.read_or_write(p_addr, v_addr, new_p_addr, is_write, data);
    uint64_t output_val = output_word.to_int();
    if (get_mode() == EVAL || get_mode() == DEBUG) {
      EXPECT_EQ(output_val, ref_mem[v_addr_val]);
    }
    if (is_write_val) {
      ref_mem[v_addr_val] = data_val;
    }
    pos_map[v_addr_val] = new_p_addr_val;
    return FuncOutput();
  };
  DEFINE_FUNC(main, {}, main_func);
  CircuitORAMMain(uint64_t T, Mode mode, uint memory_space, uint word_width,
                  const ORAMConstParams& params = ORAMConstParams::dbg_params())
      : Gadget(mode, T),
        oram(this, memory_space, word_width, T, params),
        v_addr_width(log2ceil(memory_space)),
        position_width(log2ceil(memory_space)),
        word_width(word_width),
        memory_space(memory_space) {
    pos_map.resize(memory_space);
    new_p_addr_list.resize(memory_space + T);
    for (uint64_t i = 0; i < memory_space + T; ++i) {
      new_p_addr_list[i] = i;
    }
    int prev_rand_seed = rand();
    srand(12345);  // ensure that the position map is the same for the
    std::random_shuffle(new_p_addr_list.begin(), new_p_addr_list.end());
    srand(prev_rand_seed);  // restore the prng
    std::copy(new_p_addr_list.begin(), new_p_addr_list.begin() + memory_space,
              pos_map.begin());
    oram.init_evict_counter(T % memory_space);
    new_p_addr_idx = memory_space;
    ref_mem.resize(memory_space, 0);
    set_name("CircuitORAMMain");
  }
};

TEST(CircuitORAM, RandTest) {
  test_gadget<CircuitORAMMain>(ANY, {1, 64}, 1024 * 1024 * 16, 16, 8);
}

TEST(CircuitORAM, RandDbgTest) {
  test_gadget_dbg<CircuitORAMMain>(POW2, {8, 128}, 4, 8);
}

TEST(CircuitORAM, RandMeasure) {
  measure_gadget_gc<CircuitORAMMain>(64, 32, 8);
}

TEST(CircuitORAM, MeasurePow24MainMap) {
  auto params = ORAMConstParams::compute_friendly_params(24);
  // params.simd_link_threshold = 32;
  measure_gadget_gc<CircuitORAMMain>(1UL << 24, 1UL << 24, 64, params);
}

TEST(CircuitORAM, Perf) {
  GTEST_SKIP();
  test_gadget<CircuitORAMMain>(ANY, {256, 256}, 1024 * 1024 * 1024, 256, 32);
}

struct CircuitORAMInitMain : Gadget {
  CircuitORAM oram;
  uint v_addr_width;
  uint position_width;
  uint word_width;
  uint64_t memory_space;
  std::vector<uint64_t> ref_mem;
  std::vector<uint64_t> new_p_addr_list;
  std::vector<uint64_t> pos_map;
  uint64_t new_p_addr_idx;
  std::function<FuncOutput(FuncInput)> main_func = [&](FuncInput) {
    uint64_t v_addr_val = rand() % memory_space;
    Word v_addr = Word::input_dbg(self, v_addr_width, v_addr_val);
    uint64_t p_addr_val = pos_map[v_addr_val];
    Word p_addr = Word::input_dbg(self, position_width, p_addr_val);
    uint8_t is_write_val = rand() % 2;
    Bit is_write = Bit::input_dbg(self, is_write_val);
    uint64_t data_val = rand() % (1UL << word_width);
    Word data = Word::input_dbg(self, word_width, data_val);
    uint64_t new_p_addr_val = new_p_addr_list[new_p_addr_idx++];
    Word new_p_addr = Word::input_g(self, position_width, new_p_addr_val);
    const Word& output_word =
        oram.read_or_write(p_addr, v_addr, new_p_addr, is_write, data);
    uint64_t output_val = output_word.to_int();
    if (get_mode() == EVAL || get_mode() == DEBUG) {
      EXPECT_EQ(output_val, ref_mem[v_addr_val]);
    }
    if (is_write_val) {
      ref_mem[v_addr_val] = data_val;
    }
    pos_map[v_addr_val] = new_p_addr_val;
    return FuncOutput();
  };
  DEFINE_FUNC(main, {}, main_func);
  CircuitORAMInitMain(
      uint64_t T, Mode mode, uint memory_space, uint word_width,
      const ORAMConstParams& params = ORAMConstParams::dbg_params())
      : Gadget(mode, T),
        oram(this, memory_space, word_width, T, params),
        v_addr_width(log2ceil(memory_space)),
        position_width(log2ceil(memory_space)),
        word_width(word_width),
        memory_space(memory_space) {
    pos_map.resize(memory_space);
    new_p_addr_list.resize(memory_space + T);
    for (uint64_t i = 0; i < memory_space + T; ++i) {
      new_p_addr_list[i] = i;
    }
    int prev_rand_seed = rand();
    srand(123456);  // ensure that the position map is the same for the garbler
                    // and the evaluator
    std::random_shuffle(new_p_addr_list.begin(), new_p_addr_list.end());

    std::copy(new_p_addr_list.begin(), new_p_addr_list.begin() + memory_space,
              pos_map.begin());
    new_p_addr_idx = memory_space;
    ref_mem.resize(memory_space, 0);
    std::vector<PlaintextBlocksInPos> blocks_in_pos(memory_space);
    for (uint64_t i = 0; i < memory_space; ++i) {
      ref_mem[i] = rand() % (1UL << word_width);
    }
    srand(prev_rand_seed);  // restore the prng
    for (uint64_t i = 0; i < memory_space; ++i) {
      uint64_t p_addr = pos_map[i] % (1UL << position_width);
      uint64_t pos = reverse_bits(p_addr, position_width);
      // std::cout << "assign " << pos << " with " << addr << " " << ref_mem[i]
      //           << std::endl;
      blocks_in_pos[pos].emplace_back(p_addr, i, ref_mem[i], word_width);
    }
    oram.init_ram_g(blocks_in_pos.begin(), blocks_in_pos.end());

    oram.init_evict_counter(T % memory_space);

    set_name("CircuitORAMMain");
  }
};

TEST(CircuitORAM, InitTestTiny) {
  // for valgrind
  test_gadget<CircuitORAMInitMain>(POW2, {8, 8}, 1024 * 1024, 4, 3,
                                   ORAMConstParams{10, 3, 2, 2, 1});
}

TEST(CircuitORAM, InitTest) {
  test_gadget<CircuitORAMInitMain>(POW2, {64, 64}, 1024 * 1024 * 32, 8, 8);
}

TEST(CircuitORAM, InitDbgTest) {
  test_gadget_dbg<CircuitORAMInitMain>(POW2, {512, 512}, 256, 8);
}