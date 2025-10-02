#include "util.hpp"

#include <dlfcn.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
namespace ZebraGRAM {
uint bit_width(uint64_t n) { return log2(n) + 1; }

uint64_t reverse_bits(uint64_t word, uint count) {
  uint64_t reverse_word = 0;
  for (uint i = 0; i < count; ++i) {
    reverse_word |= ((word >> i) & 1) << ((count - 1) - i);
  }
  return reverse_word;
}

uint log2(uint64_t n) {
  if (n < 256) {
    uint result = 0;
    while (n > 1) {
      n >>= 1;
      ++result;
    }
    return result;
  }
  static const uint64_t masks[6] = {0x2,    0xC,        0xF0,
                                    0xFF00, 0xFFFF0000, 0xFFFFFFFF00000000};
  uint c = 32;
  uint r = 0;
  for (uint i = 6; i;) {
    const bool cond = n & masks[--i];
    if (cond) {
      n >>= c;
      r |= c;
    }
    c >>= 1;
  }
  return r;
}

uint log2ceil(uint64_t n) { return log2(n - 1) + 1; }

uint64_t rand64() { return (uint64_t(rand()) << 32) | uint64_t(rand()); }

constexpr double LN2 = 0.693147180559945309417232121458176568;  // ln(2)

// epsilon_exp is the base-2 exponent: e.g. -60 means failure prob = 2^-60
double chernoff_upper_bound(double N, double M, double K, int epsilon_exp) {
  if (K == 0 || N == 0 || M == 0) return 0.0;

  // Number of "good" items (< N/K)
  double S = N / K;  // floor
  double mu = M * S / N;

  // log(1/epsilon) = -epsilon_exp * ln(2)
  double log_term = -static_cast<double>(epsilon_exp) * LN2;

  // Chernoff bound: mu + sqrt(3 * mu * log(1/epsilon))
  double bound = ceil(mu + std::sqrt(3.0 * mu * log_term));

  // Can't exceed min(M, S)
  return std::min(bound, static_cast<double>(std::min(M, S)));
}

}  // namespace ZebraGRAM

#include <boost/stacktrace.hpp>
#include <iostream>
void print_stack_trace() { std::cout << boost::stacktrace::stacktrace(); }

Timer ParTracker::par_timer;
bool ParTracker::parallel_enabled = false;
Timer test_timer;

#ifdef MEASURE_STACK_COST
uint64_t global_stack_cost = 0;
bool measure_stack_flag = false;
#endif

#ifdef MEASURE_TSC_STACK
uint64_t global_tsc_stack_cost = 0;
#endif

#ifdef PRINT_MEMORY_USE
#include <signal.h>

void raise_sigusr1() { raise(SIGUSR1); }
#else
void raise_sigusr1() {}
#endif

#include <flint/flint.h>
#include <omp.h>
void flint_cleanup_all_threads() {
#pragma omp parallel num_threads(omp_get_max_threads())
  {
    flint_cleanup();
  }
}

#ifdef PRINT_PROGRESS_GRANULARITY
// record the number of garbled bytes for printing progress
#include <chrono>
#include <ctime>
#include <string>

uint64_t next_gc_ckpt = 0;

void print_curr_gc_offset(uint64_t curr_gc_offset) {
  if (curr_gc_offset >= next_gc_ckpt) {
    next_gc_ckpt = curr_gc_offset / PRINT_PROGRESS_GRANULARITY *
                       PRINT_PROGRESS_GRANULARITY +
                   PRINT_PROGRESS_GRANULARITY;
    flint_cleanup_all_threads();
    auto now = std::chrono::system_clock::now();
    std::time_t now_c = std::chrono::system_clock::to_time_t(now);
    std::string time_str = std::ctime(&now_c);
    time_str.pop_back();
    std::cout << "Garbled " << curr_gc_offset
              << " bytes so far, current time: " << time_str << std::endl;
#ifdef PRINT_MEMORY_USE
    raise_sigusr1();
#endif
  }
}
#endif