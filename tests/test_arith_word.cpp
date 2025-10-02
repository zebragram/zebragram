#include "arith_word.hpp"
#include "bit.hpp"
#include "gadget.hpp"
#include "test_util.hpp"
using namespace ZebraGRAM;

// TEST(TestArithWord, Join)
// {
//     test_gadget<TestGadgetJoin<ArithWord>>(ANY, {2, 2});
// }

struct TestArithWordMain : Gadget {
  std::function<FuncOutput(FuncInput)> main_func = [&](FuncInput) {
    std::vector<int> bit_vals;
    std::vector<std::vector<uint64_t>> word_vals;
    std::vector<Bit> bits;
    std::vector<ArithWord> words;
    int batch_size = 1 + get_time() % 6;
    int word_width = 2 + get_time();
    for (int i = 0; i < batch_size; ++i) {
      bit_vals.push_back(get_time() % 2);
      bits.push_back(Bit::input_dbg(self, bit_vals.back()));
      word_vals.push_back(std::vector<uint64_t>());
      std::vector<uint64_t> &vals = word_vals.back();
      for (int j = 0; j < word_width; ++j) {
        vals.push_back((uint64_t)(get_time() * 1000));
      }
      words.push_back(ArithWord::input_dbg(self, vals));
    }

    // std::cout << "batch size: " << batch_size << ", word width: " <<
    // word_width << std::endl;
    std::vector<ArithWord> results = batch_bit_arith_word_mul(bits, words);
    // std::cout << "batch mul done for mode " << get_mode() << std::endl;
    for (int i = 0; i < batch_size; ++i) {
      std::vector<ArithLabel> revealed = results[i].reveal();
      if (get_mode() == EVAL || get_mode() == DEBUG) {
        for (int j = 0; j < word_width; ++j) {
          uint64_t expected = bit_vals[i] * word_vals[i][j];
          uint64_t actual = fmpz_get_ui(revealed[j].raw);
          fmpz_t expected_auth;
          fmpz_init(expected_auth);
          fmpz_mul(expected_auth, revealed[j].raw, default_sk.K);
          EXPECT_EQ(fmpz_cmp(expected_auth, revealed[j].auth), 0);
          fmpz_clear(expected_auth);
          EXPECT_EQ(actual, expected);
        }
      }
    }
    return FuncOutput();
  };
  DEFINE_FUNC(main, {}, main_func);
  TestArithWordMain(uint64_t T, Mode mode) : Gadget(mode, T) {}
};

struct BenchArithWordMain : Gadget {
  std::function<FuncOutput(FuncInput)> main_func = [&](FuncInput) {
    std::vector<int> bit_vals;
    std::vector<std::vector<uint64_t>> word_vals;
    std::vector<Bit> bits;
    std::vector<ArithWord> words;
    int batch_size = 4;
    int word_width = (4096 * 8 + LEN_ARITH_DIGIT - 1) / LEN_ARITH_DIGIT;
    for (int i = 0; i < batch_size; ++i) {
      bit_vals.push_back(get_time() % 2);
      bits.push_back(Bit::input_dbg(self, bit_vals.back()));
      word_vals.push_back(std::vector<uint64_t>());
      std::vector<uint64_t> &vals = word_vals.back();
      for (int j = 0; j < word_width; ++j) {
        vals.push_back((uint64_t)(get_time() * 1000));
      }
      words.push_back(ArithWord::input_dbg(self, vals));
    }

    // std::cout << "batch size: " << batch_size << ", word width: " <<
    // word_width << std::endl;
    std::vector<ArithWord> results = batch_bit_arith_word_mul(bits, words);
    // std::cout << "batch mul done for mode " << get_mode() << std::endl;
    for (int i = 0; i < batch_size; ++i) {
      std::vector<ArithLabel> revealed = results[i].reveal();
      if (get_mode() == EVAL || get_mode() == DEBUG) {
        for (int j = 0; j < word_width; ++j) {
          uint64_t expected = bit_vals[i] * word_vals[i][j];
          uint64_t actual = fmpz_get_ui(revealed[j].raw);
          fmpz_t expected_auth;
          fmpz_init(expected_auth);
          fmpz_mul(expected_auth, revealed[j].raw, default_sk.K);
          EXPECT_EQ(fmpz_cmp(expected_auth, revealed[j].auth), 0);
          fmpz_clear(expected_auth);
          EXPECT_EQ(actual, expected);
        }
      }
    }
    return FuncOutput();
  };
  DEFINE_FUNC(main, {}, main_func);
  BenchArithWordMain(uint64_t T, Mode mode) : Gadget(mode, T) {}
};

TEST(TestArithWord, Main) {
  omp_set_nested(2);
  test_gadget<TestArithWordMain>(ANY, {16, 16});
}

TEST(TestArithWord, bench) {
  paillier_keygen(default_pk, default_sk, 1024);
  omp_set_nested(2);
  test_gadget<BenchArithWordMain>(ANY, {100, 100});
}
