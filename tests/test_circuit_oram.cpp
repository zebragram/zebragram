#include <algorithm>
#include <random>

#include "circuit_oram.hpp"
#include "test_util.hpp"

using namespace ZebraGRAM;
struct CircuitORAMMain : Gadget {
  CircuitORAM oram;
  uint v_addr_width;
  uint position_width;
  uint word_width;
  uint payload_width;
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

    uint64_t new_p_addr_val = new_p_addr_list[new_p_addr_idx++];
    Word new_p_addr = Word::input_g(self, position_width, new_p_addr_val);
    // std::cout << (is_write_val ? "write " + data_val : "read ") << " from "
    //           << p_addr_val << " new addr " << new_p_addr_val << std::endl;
    uint64_t output_val = 0;
    if (payload_width) {
      Word data = Word::input_dbg(self, word_width, 0UL);
      std::vector<uint64_t> payload_width_vec(payload_width);
      for (uint i = 0; i < payload_width; ++i) {
        payload_width_vec[i] = data_val * 1000 * (i + 1);
      }
      data.set_payload(ArithWord::input_dbg(self, payload_width_vec));
      const Word& output_word =
          oram.read_or_write(p_addr, v_addr, new_p_addr, is_write, data);
      const ArithWord& payload = output_word.get_payload();
      Assert_eq(payload.width(), payload_width);
      auto revealed = payload.reveal();
      if (get_mode() == EVAL || get_mode() == DEBUG) {
        for (uint i = 0; i < payload_width; ++i) {
          try {
            char* output_c_str = fmpz_get_str(NULL, 10, revealed[i].raw);
            output_val = std::stoull(output_c_str);
            flint_free(output_c_str);
          } catch (const std::out_of_range& e) {
            std::cerr << "Out of range error" << std::endl;
            // std::cerr << "Value: " << output_c_str << std::endl;
            // output_val = UINT64_MAX; // Assign a default value or handle the
            // error as needed
          }
          // EXPECT_EQ(output_val, ref_mem[v_addr_val] * 1000 * (i + 1));
          if (output_val != ref_mem[v_addr_val] * 1000 * (i + 1)) {
            std::cerr << "Mismatch at index " << i << ": expected "
                      << ref_mem[v_addr_val] * 1000 * (i + 1) << ", got "
                      << output_val << std::endl;
          }
        }
      }
    } else {
      Word data = Word::input_dbg(self, word_width, data_val);
      const Word& output_word =
          oram.read_or_write(p_addr, v_addr, new_p_addr, is_write, data);
      output_val = output_word.to_int();
    }

    if (is_write_val) {
      ref_mem[v_addr_val] = data_val;
    }
    pos_map[v_addr_val] = new_p_addr_val;
    return FuncOutput();
  };
  DEFINE_FUNC(main, {}, main_func);
  CircuitORAMMain(uint64_t T, Mode mode, uint memory_space, uint word_width,
                  uint payload_width = 0,
                  const ORAMConstParams& params = ORAMConstParams::dbg_params())
      : Gadget(mode, T),
        oram(this, memory_space, word_width, payload_width, T, params),
        v_addr_width(log2ceil(memory_space)),
        position_width(log2ceil(memory_space)),
        word_width(word_width),
        payload_width(payload_width),
        memory_space(memory_space) {
    Assert(payload_width == 0 || params.simd_link_threshold > 4 * T);
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

const uint64_t payload_width_4k =
    (4096 * 8 + LEN_ARITH_DIGIT - 1) / LEN_ARITH_DIGIT;

const uint64_t payload_width_1k =
    (1024 * 8 + LEN_ARITH_DIGIT - 1) / LEN_ARITH_DIGIT;

TEST(CircuitORAM, RandPayloadTestTiny) {
  test_gadget<CircuitORAMMain>(ANY, {8, 8}, 1024 * 1024 * 16, 8, 1, 1,
                               ORAMConstParams::dbg_no_simd_params());
}

TEST(CircuitORAM, RandPayloadMedium) {
  ORAMConstParams params = ORAMConstParams::recommended_params(64, UINT32_MAX);
  test_gadget<CircuitORAMMain>(ANY, {64, 64}, 1024UL * 1024UL * 2048UL, 64, 1,
                               payload_width_4k, params);
}

TEST(CircuitORAM, RandPayloadLarge) {
  ORAMConstParams params =
      ORAMConstParams::recommended_params(4096, UINT32_MAX);
  test_gadget<CircuitORAMMain>(ANY, {4096, 4096},
                               1024UL * 1024UL * 1048UL * 300UL, 4096, 1,
                               payload_width_4k, params);
}

TEST(CircuitORAM, RandPayloadHuge) {
  paillier_keygen(default_pk, default_sk, 4096);
  ORAMConstParams params =
      ORAMConstParams::recommended_params(16384, UINT32_MAX);
  test_gadget<CircuitORAMMain>(ANY, {16384, 16384},
                               1024UL * 1024UL * 1024UL * 800UL, 16384, 1,
                               payload_width_4k, params);
}

TEST(CircuitORAM, Payload256M4k) {
  paillier_keygen(default_pk, default_sk, 4096);
  ORAMConstParams params =
      ORAMConstParams::recommended_params(65536, UINT32_MAX);
  test_gadget<CircuitORAMMain>(ANY, {65536, 65536},
                               1024UL * 1024UL * 1024L * 200UL, 65536, 1,
                               payload_width_4k, params);
}

TEST(CircuitORAM, RandPayload1kLarge) {
  paillier_keygen(default_pk, default_sk, 8192);
  ORAMConstParams params =
      ORAMConstParams::recommended_params(65536 * 4, UINT32_MAX);
  test_gadget<CircuitORAMMain>(ANY, {16384, 16384},
                               1024UL * 1024UL * 1024L * 200UL, 16384, 1,
                               payload_width_1k, params);
}

TEST(CircuitORAM, Payload256M1k) {
  paillier_keygen(default_pk, default_sk, 8192);
  ORAMConstParams params =
      ORAMConstParams::recommended_params(65536 * 4, UINT32_MAX);
  test_gadget<CircuitORAMMain>(ANY, {65536 * 4, 65536 * 4},
                               1024UL * 1024UL * 1024L * 300UL, 65536 * 4, 1,
                               payload_width_1k, params);
}

TEST(CircuitORAM, RandDbgTest) {
  test_gadget_dbg<CircuitORAMMain>(POW2, {8, 128}, 4, 8, 0);
}

TEST(CircuitORAM, RandMeasure) {
  measure_gadget_gc<CircuitORAMMain>(64, 32, 8, 0);
}

TEST(CircuitORAM, MeasurePow24MainMap) {
  auto params = ORAMConstParams::compute_friendly_params(24);
  // params.simd_link_threshold = 32;
  measure_gadget_gc<CircuitORAMMain>(1UL << 24, 1UL << 24, 64, 0, params);
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
        oram(this, memory_space, word_width, 0, T, params),
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