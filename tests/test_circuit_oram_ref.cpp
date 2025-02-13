#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>

#include "circuit_oram_ref.hpp"
using namespace PicoGRAM;
struct CircuitORAMRefMain {
  CircuitORAMRef oram;
  uint v_addr_width;
  uint position_width;
  uint word_width;
  uint64_t memory_space;
  std::vector<uint64_t> ref_mem;
  std::vector<uint64_t> new_p_addr_list;
  std::vector<uint64_t> pos_map;
  uint64_t new_p_addr_idx;
  std::function<FuncOutput(FuncInput)> main = [&](FuncInput) {
    uint64_t v_addr_val = rand() % memory_space;
    Word v_addr = v_addr_val;
    uint64_t p_addr_val = pos_map[v_addr_val];
    Word p_addr = p_addr_val;
    uint8_t is_write_val = rand() % 2;
    Bit is_write = is_write_val;
    uint64_t data_val = rand() % (1UL << word_width);
    Word data = data_val;
    uint64_t new_p_addr_val = new_p_addr_list[new_p_addr_idx++];
    Word new_p_addr = new_p_addr_val;
    // std::cout << (is_write_val ? "write " + data_val : "read ") << " from "
    //           << p_addr_val << " new addr " << new_p_addr_val << std::endl;
    const Word& output_word =
        oram.read_or_write(p_addr, v_addr, new_p_addr, is_write, data);
    uint64_t output_val = output_word;

    EXPECT_EQ(output_val, ref_mem[v_addr_val]);

    if (is_write_val) {
      ref_mem[v_addr_val] = data_val;
    }
    pos_map[v_addr_val] = new_p_addr_val;
    return FuncOutput();
  };

  CircuitORAMRefMain(uint64_t T, uint memory_space, uint word_width,
                     uint stash_size, uint bkt_size, uint evict_freq,
                     uint log_way)
      : oram(T, memory_space, word_width, stash_size, bkt_size, evict_freq,
             log_way),
        v_addr_width(log2ceil(memory_space)),
        position_width(log2ceil(memory_space)),
        word_width(word_width),
        memory_space(memory_space) {
    Assert(stash_size <= 64);
    Assert(bkt_size <= 64);
    pos_map.resize(memory_space);
    new_p_addr_list.resize(memory_space + T);
    for (uint64_t i = 0; i < memory_space + T; ++i) {
      new_p_addr_list[i] = i % memory_space;
    }
    std::random_shuffle(new_p_addr_list.begin(), new_p_addr_list.end());

    std::copy(new_p_addr_list.begin(), new_p_addr_list.begin() + memory_space,
              pos_map.begin());
    new_p_addr_idx = memory_space;
    ref_mem.resize(memory_space, 0);
  }

  uint stash_load() const { return oram.stash_load(); }
};

#include <utility>
#include <vector>

// Linear regression function
// X = 1, 2, ..., n
std::pair<double, double> linearRegression(const std::vector<double>& Y) {
  // Ensure there is data
  Assert(!Y.empty());

  size_t n = Y.size();
  double sumX = 0.0;
  double sumY = 0.0;
  double sumXY = 0.0;
  double sumX2 = 0.0;

  // Calculate sums
  for (size_t i = 0; i < n; ++i) {
    double x = static_cast<double>(i + 1);  // X = 1, 2, ..., n
    sumX += x;
    sumY += Y[i];
    sumXY += x * Y[i];
    sumX2 += x * x;
  }

  // Calculate slope (m) and intercept (b)
  double denominator = n * sumX2 - sumX * sumX;
  if (denominator == 0) {
    throw std::runtime_error(
        "Denominator is zero, cannot perform linear regression.");
  }

  double slope = (n * sumXY - sumX * sumY) / denominator;
  double intercept = (sumY - slope * sumX) / n;

  return {slope, intercept};
}

TEST(CircuitORAMRef, RandTest) {
  srand(time(0));
  uint64_t T = 1UL << 22;
  uint stash_size = 63;
  CircuitORAMRefMain main(T, 4096, 1, stash_size, 5, 2, 2);
  std::vector<uint64_t> stash_load_stats(64, 0);
  for (uint64_t t = 0; t < T; ++t) {
    main.main({});
    uint load = main.stash_load();
    ASSERT_LT(load, 64);
    ++stash_load_stats[load];
  }
  for (uint i = 0; i < stash_size; ++i) {
    uint64_t count = stash_load_stats[i];
    double ratio = count / (double)T;
    double log2_ratio = std::log2(ratio);
    std::cout << "stash load " << i << ": " << stash_load_stats[i]
              << ", log2 freq = " << log2_ratio << std::endl;
  }
  std::vector<double> log2_freqs;
  for (uint i = 1; i < stash_size; ++i) {
    if (stash_load_stats[i] <= 20) {
      break;
    }
    log2_freqs.push_back(std::log2(stash_load_stats[i] / (double)T));
  }
  auto slope_and_intercept = linearRegression(log2_freqs);
  double slope = slope_and_intercept.first;
  double intercept = slope_and_intercept.second;
  std::cout << "slope: " << slope << ", intercept: " << intercept << std::endl;
  double target_prob = -60.0;  // target failure probability per access
  double target_overflow_one_prob =
      target_prob + std::log2(1.0 - std::pow(2.0, slope));
  std::cout << "target overflow one prob: " << target_overflow_one_prob
            << std::endl;
  double min_stash_size = (target_overflow_one_prob - intercept) / slope;
  std::cout << "min stash size: " << min_stash_size << std::endl;
}