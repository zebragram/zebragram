#include "gadget.hpp"
#include "test_util.hpp"
#include "word.hpp"
using namespace PicoGRAM;
struct TestWordMain : Gadget {
  std::function<FuncOutput(FuncInput)> main_func = [&](FuncInput) {
    uint len_a = 1 + get_time() % 30;
    uint len_b = 1 + get_time() % 29;
    uint val_a = rand() % (1 << len_a);
    uint val_b = rand() % (1 << len_b);
    Word a = Word::input_dbg(self, len_a, val_a);
    Word b = Word::input_dbg(self, len_b, val_b);
    // test unary operations
    Bit not_a = !a;
    Word neg_a = ~a;
    Word count_one_a = a.count_ones();
    Word count_zero_a = a.count_zeros();
    Word trailing_zeros_a = a.trailing_zeros();
    uint not_a_val = not_a.to_int();
    uint neg_a_val = neg_a.to_int();
    uint count_one_a_val = count_one_a.to_int();
    uint count_zero_a_val = count_zero_a.to_int();
    uint trailing_zeros_a_val = trailing_zeros_a.to_int();
    if (get_mode() == EVAL || get_mode() == DEBUG) {
      EXPECT_EQ(not_a_val, !val_a);
      EXPECT_EQ(neg_a_val, (~val_a) & ((1 << len_a) - 1));
      EXPECT_EQ(count_one_a_val, __builtin_popcount(val_a));
      EXPECT_EQ(count_zero_a_val, a.width() - __builtin_popcount(val_a));
      if (val_a == 0) {
        EXPECT_EQ(trailing_zeros_a_val, a.width());
      } else {
        EXPECT_EQ(trailing_zeros_a_val, __builtin_ctz(val_a));
      }
    }

    // test add, and bitwise binary operations
    Word sum = a + b;
    Word xor_ = a ^ b;
    Word and_ = a & b;
    Word or_ = a | b;
    uint sum_val = sum.to_int();
    uint xor_val = xor_.to_int();
    uint and_val = and_.to_int();
    uint or_val = or_.to_int();
    if (get_mode() == EVAL || get_mode() == DEBUG) {
      EXPECT_EQ(sum_val, val_a + val_b);
      EXPECT_EQ(xor_val, val_a ^ val_b);
      EXPECT_EQ(and_val, val_a & val_b);
      EXPECT_EQ(or_val, val_a | val_b);
    }
    // test mult bit
    uint8_t bit_val = rand() % 2;
    Bit mult_bit = Bit::input_dbg(self, bit_val);
    Word mult = a * mult_bit;
    uint mult_val = mult.to_int();
    if (get_mode() == EVAL || get_mode() == DEBUG) {
      EXPECT_EQ(mult_val, val_a * bit_val);
    }

    uint64_t const_val = get_time();
    // test compare
    Bit is_less = a < b;
    Bit is_eq = a == b;
    Bit is_eq_const = a == const_val;
    Bit is_neq = a != b;
    Bit is_neq_const = a != const_val;
    Bit is_greater = a > b;
    Bit is_leq = a <= b;
    Bit is_geq = a >= b;
    uint is_less_val = is_less.to_int();
    uint is_eq_val = is_eq.to_int();
    uint is_eq_const_val = is_eq_const.to_int();
    uint is_neq_val = is_neq.to_int();
    uint is_neq_const_val = is_neq_const.to_int();
    uint is_greater_val = is_greater.to_int();
    uint is_leq_val = is_leq.to_int();
    uint is_geq_val = is_geq.to_int();
    if (get_mode() == EVAL || get_mode() == DEBUG) {
      EXPECT_EQ(is_less_val, val_a < val_b);
      EXPECT_EQ(is_eq_val, val_a == val_b);
      EXPECT_EQ(is_eq_const_val, val_a == const_val);
      EXPECT_EQ(is_neq_val, val_a != val_b);
      EXPECT_EQ(is_neq_const_val, val_a != const_val);
      EXPECT_EQ(is_greater_val, val_a > val_b);
      EXPECT_EQ(is_leq_val, val_a <= val_b);
      EXPECT_EQ(is_geq_val, val_a >= val_b);
    }

    // test shift
    uint shift = get_time() % 10;
    Word left_shift = a << shift;
    Word right_shift = a >> shift;
    uint64_t left_shift_val = left_shift.to_int();
    uint64_t right_shift_val = right_shift.to_int();
    if (get_mode() == EVAL || get_mode() == DEBUG) {
      EXPECT_EQ(left_shift_val, (uint64_t)val_a << shift);
      EXPECT_EQ(right_shift_val, (uint64_t)val_a >> shift);
    }

    // test conditional swap
    uint8_t cond_val = rand() % 2;
    Bit cond = Bit::input_dbg(self, cond_val);
    Word swap_a = a;
    Word swap_b = b;
    Word::cond_swap(cond, swap_a, swap_b);
    uint swap_a_val = swap_a.to_int();
    uint swap_b_val = swap_b.to_int();
    if (get_mode() == EVAL || get_mode() == DEBUG) {
      if (cond_val) {
        EXPECT_EQ(swap_a_val, val_b);
        EXPECT_EQ(swap_b_val, val_a);
      } else {
        EXPECT_EQ(swap_a_val, val_a);
        EXPECT_EQ(swap_b_val, val_b);
      }
    }
    return FuncOutput();
  };
  DEFINE_FUNC(main, {}, main_func);
  TestWordMain(uint64_t T, Mode mode) : Gadget(mode, T) {}
};

TEST(TestWord, Main) { test_gadget<TestWordMain>(ANY, {900, 900}); }
TEST(TestWord, MainDbg) { test_gadget_dbg<TestWordMain>(ANY, {900, 900}); }
TEST(TestWord, Measure) { measure_gadget_gc<TestWordMain>(900); }

struct TestWordSimple : Gadget {
  std::function<FuncOutput(FuncInput)> main_func = [&](FuncInput) {
    Bit a = Bit::input_dbg(self, t % 2);
    Bit b = Bit::input_dbg(self, (t >> 1) % 2);
    Bit c = a & b;
    bool c_val = c.to_int();
    if (get_mode() == EVAL || get_mode() == DEBUG) {
      EXPECT_EQ(c_val, (t % 2) & ((t >> 1) % 2));
    }

    return FuncOutput();
  };
  DEFINE_FUNC(main, {}, main_func);
  TestWordSimple(uint64_t T, Mode mode) : Gadget(mode, T) {}
};

TEST(TestWord, Simple) { test_gadget<TestWordSimple>(ANY, {4, 4}); }