#include <gtest/gtest.h>

#include "util.hpp"
using namespace PicoGRAM;
TEST(TestUtil, Log) {
  srand(time(0));
  for (uint64_t r = 0; r < 10000; ++r) {
    uint max_bit_width = (rand() % 64) + 1;
    uint64_t n = rand64() & ((1UL << max_bit_width) - 1UL);
    uint log_n = PicoGRAM::log2(n);
    if (n) {
      ASSERT_LE(1UL << log_n, n);
    } else {
      ASSERT_EQ(log_n, 0);
    }
    ASSERT_GE((1UL << (log_n + 1)) - 1, n);
  }
}

TEST(TestUtil, ReverseBits) {
  srand(time(0));
  for (uint64_t r = 0; r < 10000; ++r) {
    uint max_bit_width = (rand() % 64) + 1;
    uint64_t n = rand64() & ((1UL << max_bit_width) - 1UL);
    uint width = bit_width(n);
    uint64_t reverse_n = reverse_bits(n, width);
    ASSERT_EQ(n, reverse_bits(reverse_n, width));
    uint rotate = rand() % max_bit_width;
    uint64_t rotate_n = reverse_bits(
        reverse_bits(n, rotate) |
            (reverse_bits(n >> rotate, max_bit_width - rotate) << rotate),
        max_bit_width);

    uint64_t rotate_n_ref = (n >> rotate) | ((n << (max_bit_width - rotate)) &
                                             ((1UL << max_bit_width) - 1));
    ASSERT_EQ(rotate_n, rotate_n_ref);
  }
}