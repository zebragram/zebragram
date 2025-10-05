#include "rec_oram.hpp"
#include "test_util.hpp"
using namespace ZebraGRAM;
struct RecORAMMain : Gadget {
  RecursiveORAM oram;
  uint addr_width;
  uint word_width;
  uint64_t memory_space;
  std::vector<uint64_t> ref_mem;
  std::function<FuncOutput(FuncInput)> main_func = [&](FuncInput) {
    uint64_t addr_val = rand() % memory_space;
    Word addr = Word::input_dbg(self, addr_width, addr_val);
    uint8_t is_write_val = rand() % 2;
    Bit is_write = Bit::input_dbg(self, is_write_val);
    uint64_t data_val = rand() % (1UL << word_width);
    Word data = Word::input_dbg(self, word_width, data_val);

    // std::cout << (is_write_val ? "write " + data_val : "read ") << " from "
    //           << p_addr_val << " new addr " << new_addr_val << std::endl;
    const Word& output_word = oram.read_or_write(addr, is_write, data);
    uint64_t output_val = output_word.to_int();
    if (get_mode() == EVAL || get_mode() == DEBUG) {
      EXPECT_EQ(output_val, ref_mem[addr_val]);
      //   std::cout << "output val " << output_val << ", ref " <<
      //   ref_mem[addr_val]
      //             << std::endl;
    }
    if (is_write_val) {
      ref_mem[addr_val] = data_val;
    }
    return FuncOutput();
  };
  DEFINE_FUNC(main, {}, main_func);
  RecORAMMain(uint64_t T, Mode mode, uint memory_space, uint word_width,
              uint64_t simd_link_threshold = 8,
              uint64_t linear_oram_threshold = 4)
      : Gadget(mode, T),
        oram(this, memory_space, word_width, 0, T, simd_link_threshold,
             linear_oram_threshold),
        addr_width(log2ceil(memory_space)),
        word_width(word_width),
        memory_space(memory_space) {
    oram.init_empty_ram();

    ref_mem.resize(memory_space, 0);
    set_name("CircuitORAMMain");
  }
};

struct RecORAMInitMain : Gadget {
  RecursiveORAM oram;
  uint addr_width;
  uint word_width;
  uint64_t memory_space;
  std::vector<uint64_t> ref_mem;
  std::function<FuncOutput(FuncInput)> main_func = [&](FuncInput) {
    if (t == 0) {
      std::vector<Word> data;
      data.reserve(memory_space);
      for (uint64_t i = 0; i < memory_space; ++i) {
        if (i % 2 == 0) {
          data.emplace_back(Word::input_g(this, word_width, ref_mem[i]));
        } else {
          data.emplace_back(Word::input_dbg(this, word_width, ref_mem[i]));
        }
      }
      oram.init_ram(data);
    }
    uint64_t addr_val = rand() % memory_space;
    Word addr = Word::input_dbg(self, addr_width, addr_val);
    uint8_t is_write_val = rand() % 2;
    Bit is_write = Bit::input_dbg(self, is_write_val);
    uint64_t data_val = rand() % (1UL << word_width);
    Word data = Word::input_dbg(self, word_width, data_val);

    // std::cout << (is_write_val ? "write " + data_val : "read ") << " from "
    //           << p_addr_val << " new addr " << new_addr_val << std::endl;
    const Word& output_word = oram.read_or_write(addr, is_write, data);
    uint64_t output_val = output_word.to_int();
    if (get_mode() == EVAL || get_mode() == DEBUG) {
      EXPECT_EQ(output_val, ref_mem[addr_val]);
      //   std::cout << "output val " << output_val << ", ref " <<
      //   ref_mem[addr_val]
      //             << std::endl;
    }
    if (is_write_val) {
      ref_mem[addr_val] = data_val;
    }
    return FuncOutput();
  };
  DEFINE_FUNC(main, {}, main_func);
  RecORAMInitMain(uint64_t T, Mode mode, uint memory_space, uint word_width,
                  uint64_t simd_link_threshold = 16,
                  uint64_t linear_oram_threshold = 4)
      : Gadget(mode, T),
        oram(this, memory_space, word_width, 0, T, simd_link_threshold,
             linear_oram_threshold),
        addr_width(log2ceil(memory_space)),
        word_width(word_width),
        memory_space(memory_space) {
    ref_mem.resize(memory_space);
    for (uint64_t i = 0; i < memory_space; ++i) {
      ref_mem[i] = rand() % (1UL << word_width);
    }
    set_name("CircuitORAMInitMain");
  }
};

TEST(RecORAM, RandDbgTest) {
  test_gadget_dbg<RecORAMMain>(POW2, {512, 512}, 64, 8);
}

TEST(RecORAM, MeasureSmall) {
  measure_gadget_gc<RecORAMMain>(1024, 256, 16, 8, 512);
}

TEST(RecORAM, MeasureMedium) {
  GTEST_SKIP();
  measure_gadget_gc<RecORAMMain>(8192, 4096, 32, 8, 512);
}

TEST(RecORAM, MeasureNonSIMDMedium) {
  GTEST_SKIP();
  measure_gadget_gc<RecORAMMain>(8192, 4096, 32, UINT32_MAX, 512);
}

TEST(RecORAM, TestMeasureCorrectness) {
  int seed = time(0);
  // std::cout << "seed: " << seed << std::endl;
  srand(seed);
  uint64_t T = rand() % 256 + 1;
  uint memory_space = 1UL << (rand() % 6 + 1);
  uint word_width = rand() % 8 + 1;
  std::cout << "T: " << T << ", memory space: " << memory_space
            << ", word width: " << word_width << std::endl;

  GCPtr measure_gc(-1);  // dummy storage
  RecORAMMain rec_oram_measure(T, MEASURE, memory_space, word_width);
  uint64_t measure_gc_size = rec_oram_measure.garble(measure_gc).get_offset();
  std::vector<GCPtr> measure_gc_ptrs = rec_oram_measure.get_init_gc_ptrs();

  std::cout << "measure gc size: " << measure_gc_size << " bytes" << std::endl;
  int fid = open_file("rec_oram_test_data.bin", 1024 * 1024 * 128);
  GCPtr real_gc(fid);
  RecORAMMain rec_oram_real(T, GARBLE, memory_space, word_width);
  uint64_t real_gc_size = rec_oram_real.garble(real_gc).get_offset();
  std::vector<GCPtr> real_gc_ptrs = rec_oram_real.get_init_gc_ptrs();

  EXPECT_EQ(measure_gc_size, real_gc_size);
  close_file(fid);
}

TEST(RecORAM, RandTest) {
  test_gadget<RecORAMMain>(ANY, {1, 128}, 1024 * 1024 * 128, 32, 8);
}

TEST(RecORAM, RandInitTest) {
  test_gadget<RecORAMInitMain>(ANY, {1, 128}, 1024 * 1024 * 8, 32, 8);
}

TEST(RecORAM, PerfSIMDInitSmall) {
  test_gadget_with_workers<RecORAMInitMain>(ANY, {512, 512}, 1UL << 30, 4, 512,
                                            32, 32, 512);
}

TEST(RecORAM, PerfSIMDSmall) {
  test_gadget_with_workers<RecORAMMain>(ANY, {512, 512}, 1UL << 30, 4, 512, 32,
                                        32, 512);
}

TEST(RecORAM, PerfSIMDMedium4threads) {
  test_gadget_with_workers<RecORAMMain>(ANY, {4096, 4096}, 8UL << 30, 4, 4096,
                                        32, 32, 512);
}

TEST(RecORAM, PerfSIMDMedium8threads) {
  test_gadget_with_workers<RecORAMMain>(ANY, {4096, 4096}, 8UL << 30, 8, 4096,
                                        32, 32, 512);
}

/* ================== Performance tests in paper ================== */
// Warning: 1TB memory required to run the largest benchmark
// alternatively, change STORAGE_TYPE to DISK in global.hpp
// measure ZebraGRAM communication
TEST(RecORAM, IncrPerfSIMDCommCost) {
  for (uint64_t N = 1024; N <= (1UL << 24); N *= 2) {
    measure_gadget_gc<RecORAMMain>(N, N, 64, 32, 512);
  }
}

// Picogram perf
TEST(RecORAM, IncrPerfSIMD8threads) {
  /**
   */
  for (uint64_t N = 1024; N <= (1UL << 16); N *= 2) {
    test_gadget_with_workers<RecORAMMain>(ANY, {N, N}, N << 22, 7, N, 64, 32,
                                          512);
  }
}

// the ZebraGRAM no SIMD baseline communication
TEST(RecORAM, IncrPerfNonSIMDCommCost) {
  for (uint64_t N = 1024; N <= (1UL << 24); N *= 2) {
    measure_gadget_gc<RecORAMMain>(N, N, 64, UINT32_MAX, 512);
  }
}

// Picogram communication at 2^16 for different parameters
TEST(RecORAM, IncrMeasureSIMDAdaptive65536) {
  uint64_t N = (1UL << 16);
  for (uint64_t threshold = 8; threshold <= N * 2; threshold *= 4) {
    measure_gadget_gc<RecORAMMain>(N, N, 64, threshold, 512);
  }
}

// Picogram perf at 2^16 for different parameters
TEST(RecORAM, IncrPerfSIMD8threadsAdaptive65536) {
  uint64_t N = (1UL << 16);
  for (uint64_t threshold = 8; threshold <= N * 2; threshold *= 4) {
    uint64_t gc_size = 0;
    switch (threshold) {
      case 128:
        gc_size = 265UL << 30;
        break;
      case 512:
        gc_size = 290UL << 30;
        break;
      case 2048:
        gc_size = 330UL << 30;
        break;
      case 8192:
        gc_size = 385UL << 30;
        break;
      case 32768:
        gc_size = 449UL << 30;
        break;
      case 131072:
        gc_size = 610UL << 30;
        break;
      default:
        gc_size = 400UL << 30;
        break;
    }
    test_gadget_with_workers<RecORAMMain>(ANY, {N, N}, gc_size, 7, N, 64,
                                          threshold, 512);
  }
}
