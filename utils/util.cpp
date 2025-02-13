#include "util.hpp"

#include <dlfcn.h>

#include <fstream>
#include <iomanip>
namespace PicoGRAM {
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

}  // namespace PicoGRAM

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
