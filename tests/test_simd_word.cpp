#include "gadget.hpp"
#include "simd_word.hpp"
#include "test_util.hpp"
using namespace ZebraGRAM;
struct TestSIMDWordMain : Gadget {
  uint max_word_width, min_word_width;
  std::function<FuncOutput(FuncInput)> main_func = [&](FuncInput) {
    uint word_width =
        min_word_width + get_time() % (max_word_width + 1 - min_word_width);
    uint64_t word_val = rand() % (1UL << word_width);
    Word word = Word::input_dbg(this, word_width, word_val);
    SIMDWord simd_word(word, 0);
    Word convert_back_word = simd_word.to_word();
    Assert_eq(convert_back_word.width(), word_width);
    for (uint i = 0; i < word_width; ++i) {
      Assert(!convert_back_word[i].is_pub_e());
    }
    uint64_t convert_back_val = convert_back_word.to_int();
    if (get_mode() == EVAL || get_mode() == DEBUG) {
      EXPECT_EQ(convert_back_val, word_val);
    }

    uint shift_amount = (get_time()) % word_width;
    simd_word >>= shift_amount;
    Word convert_back_shifted_word = simd_word.to_word();
    uint64_t convert_back_shifted_val = convert_back_shifted_word.to_int();
    if (get_mode() == EVAL || get_mode() == DEBUG) {
      EXPECT_EQ(convert_back_shifted_val, word_val >> shift_amount);
    }

    return FuncOutput();
  };
  DEFINE_FUNC(main, {}, main_func);
  TestSIMDWordMain(uint64_t T, Mode mode, uint max_word_width = 30,
                   uint min_word_width = 2)
      : Gadget(mode, T),
        max_word_width(max_word_width),
        min_word_width(min_word_width) {}
};

TEST(TestSIMDWord, MainTest) { test_gadget<TestSIMDWordMain>(ANY, {100, 100}); }
TEST(TestSIMDWord, MainDbgTest) {
  test_gadget_dbg<TestSIMDWordMain>(ANY, {900, 900});
}

struct TestSIMDWordPerfMain : Gadget {
  uint max_word_width, min_word_width;
  std::function<FuncOutput(FuncInput)> main_func = [&](FuncInput) {
    uint word_width =
        min_word_width + get_time() % (max_word_width + 1 - min_word_width);
    uint64_t word_val = rand() % (1UL << word_width);
    Word word = Word::input_dbg(this, word_width, word_val);
    SIMDWord simd_word(word, 0);
    Word convert_back_word = simd_word.to_word();
    Assert_eq(convert_back_word.width(), word_width);
    for (uint i = 0; i < word_width; ++i) {
      Assert(!convert_back_word[i].is_pub_e());
    }

    return FuncOutput();
  };
  DEFINE_FUNC(main, {}, main_func);
  TestSIMDWordPerfMain(uint64_t T, Mode mode, uint max_word_width = 30,
                       uint min_word_width = 2)
      : Gadget(mode, T),
        max_word_width(max_word_width),
        min_word_width(min_word_width) {}
};

TEST(TestSIMDWord, Perf) {
  //   T: 50000
  // gc size: 190.326 MB
  // garbling time: 89775 ms
  // eval time: 57877 ms

  // without pow
  // garbling time: 16069 ms
  // eval time: 10315 ms

  // pow twice
  // garbling time: 139653 ms
  // eval time: 91470 ms

  // 2 threads
  // garbling time: 60893 ms
  // eval time: 43606 ms

  // 4 threads
  // garbling time: 49612 ms
  // eval time: 37644 ms
  // GTEST_SKIP();
  test_gadget<TestSIMDWordPerfMain>(ANY, {5000, 5000}, 64, 64);
  // SIMDWord::stop_workers_g();
}

TEST(TestSIMDWord, Worker) {
  SIMDWord::start_workers_g(4);
  SIMDWord::stop_workers_g();
}

struct TestSIMDWordAggrMain : Gadget {
  std::function<FuncOutput(FuncInput)> main_func = [&](FuncInput) {
    uint num_words = 1 + t % 10;
    std::vector<SIMDWord> simd_words;
    uint64_t ref_val = 0;
    uint64_t bit_offset = 0;
    uint non_zero_word_idx = t % num_words;
    uint word_width = t % 6 + 1;
    BigInt label = BigInt::random();
    for (uint i = 0; i < num_words; ++i) {
      uint64_t word_val =
          i == non_zero_word_idx ? rand() % (1UL << word_width) : 0;
      ref_val |= word_val;
      Word word = Word::input_dbg(this, word_width, word_val);
      simd_words.emplace_back(word, bit_offset);
      bit_offset += word_width;
    }
    std::vector<SIMDWord> aggr_simd_words =
        SIMDWord::aggr_simd_words(simd_words, 1);
    std::vector<Word> aggr_words = SIMDWord::to_words(aggr_simd_words);
    Assert_eq(aggr_words.size(), 1);
    uint64_t aggr_val = aggr_words[0].to_int();
    if (get_mode() == EVAL || get_mode() == DEBUG) {
      EXPECT_EQ(aggr_val, ref_val);
    }
    return FuncOutput();
  };
  DEFINE_FUNC(main, {}, main_func);
  TestSIMDWordAggrMain(uint64_t T, Mode mode) : Gadget(mode, T) {}
};

TEST(TestSIMDWord, AggrMain) {
  test_gadget<TestSIMDWordAggrMain>(ANY, {10, 100});
}