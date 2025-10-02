#include <gtest/gtest.h>

#include <cmath>

#include "util.hpp"
using namespace ZebraGRAM;
TEST(TestUtil, Log) {
  srand(time(0));
  for (uint64_t r = 0; r < 10000; ++r) {
    uint max_bit_width = (rand() % 64) + 1;
    uint64_t n = rand64() & ((1UL << max_bit_width) - 1UL);
    uint log_n = ZebraGRAM::log2(n);
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

TEST(TestUtil, Chernoff) {
  std::cout << chernoff_upper_bound(1UL << 20, 1UL << 10, 128) << std::endl;
  std::cout << chernoff_upper_bound(1UL << 20, 1UL << 19, 128) << std::endl;
  std::cout << chernoff_upper_bound(1UL << 20, 1UL << 19, 1UL << 18)
            << std::endl;
  // consider circuit oram scenario
  uint64_t N = 1UL << 17;  // number of leaves
  uint64_t M = 1UL << 16;  // number of real accesses
  uint64_t way = 4;
  std::vector<uint64_t> bounds;
  std::vector<uint64_t> naive_bounds;
  for (uint64_t K = 1; K <= N; K *= way) {
    uint64_t bound = ceil(chernoff_upper_bound(N, M, K));
    bounds.push_back(bound);
    naive_bounds.push_back(std::min(M, (uint64_t)ceil(N / K)));
    std::cout << "K = " << K << ", bound = " << bound
              << ", naive_bound = " << naive_bounds.back() << std::endl;
  }
  for (uint64_t i = 0; i < bounds.size(); ++i) {
    ASSERT_LE(bounds[i], naive_bounds[i]);
  }
}