#pragma once
#include <gtest/gtest.h>

#include <chrono>

#include "gadget.hpp"
#include "gc_ptr.hpp"
#include "global.hpp"
#include "simd_word.hpp"
#include "util.hpp"
namespace ZebraGRAM {

// Helper function to print size with 4 decimal places using integer operations
// only
static void print_size_with_unit(uint64_t size_bytes, const char* unit,
                                 uint64_t divisor) {
  uint64_t integer_part = size_bytes / divisor;
  uint64_t remainder = size_bytes % divisor;
  uint64_t decimal_part = (remainder * 10000) / divisor;  // 4 decimal places
  std::cout << integer_part << "." << std::setfill('0') << std::setw(4)
            << decimal_part << " " << unit << std::endl;
}

static void print_gc_size(uint64_t gc_size) {
  if (gc_size > 1UL << 30) {
    print_size_with_unit(gc_size, "GB", 1UL << 30);
  } else if (gc_size > 1UL << 20) {
    print_size_with_unit(gc_size, "MB", 1UL << 20);
  } else {
    print_size_with_unit(gc_size, "KB", 1UL << 10);
  }
}

enum TRangeType { POW2, EVEN, ANY };

static uint64_t get_rand_T(TRangeType T_range_type,
                           std::pair<uint32_t, uint32_t> T_range) {
  uint64_t T;
  if (T_range.first >= T_range.second) {
    T = T_range.first;
  } else {
    switch (T_range_type) {
      case POW2: {
        uint log_T_min = log2ceil(T_range.first);
        uint log_T_max = bit_width(T_range.second);
        uint log_T = rand() % (log_T_max - log_T_min) + log_T_min;
        T = 1UL << log_T;
      } break;
      case EVEN:
        T = (rand() % (T_range.second - T_range.first) + T_range.first) / 2 * 2;
        if (T < T_range.first) {
          T = T_range.first + 1;
        }
        break;
      case ANY:
        T = rand() % (T_range.second - T_range.first) + T_range.first;
        break;
      default:
        T = T_range.first;
    }
  }
  return T;
}

template <typename MainGadget, typename... Args>
void test_gadget_with_workers(TRangeType T_range_type = ANY,
                              std::pair<uint32_t, uint32_t> T_range = {1, 1000},
                              uint64_t buffer_size = 1024 * 1024 * 512,
                              uint num_workers = 0, Args... args) {
  uint64_t rand_seed = time(0);
  srand(rand_seed);
  uint64_t T = get_rand_T(T_range_type, T_range);
  std::cout << "T: " << T << std::endl;
  srand(rand_seed + 1);
  int fid = open_file("test_data.bin", (size_t)buffer_size);
  ASSERT_GE(fid, 0);
  GCPtr gc_ptr(fid);
  // timing garbling time
  std::vector<GCPtr> gc_ptrs;
  SIMDWord::start_workers_g(num_workers);
  _fmpz_cleanup();
  std::cout << "Before garbling: " << std::endl;
  raise_sigusr1();
  {
    MainGadget parent_g(T, GARBLE, args...);
    std::cout << "After creating gadget: " << std::endl;
    raise_sigusr1();
    auto start = std::chrono::high_resolution_clock::now();
    GCPtr output = parent_g.garble(gc_ptr);
    auto end = std::chrono::high_resolution_clock::now();
    SIMDWord::stop_workers_g();
    uint64_t gc_size = output.get_offset();
#ifdef MAX_GADGET_TIME
    std::cout << "GC used for testing: ";
    print_gc_size(gc_size);
    gc_size += overwritten_bytes;
    overwritten_bytes = 0;
#endif
    std::cout << "Actual GC used: ";
    print_gc_size(gc_size);
    double garble_time =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    std::cout << "garbling time: " << garble_time << " ms" << std::endl;
    std::cout << "avg time per time step: " << garble_time / T << " ms"
              << std::endl;

    std::cout << "parallel time: " << ParTracker::get_time() * 1e3 << " ms"
              << std::endl;
    ParTracker::reset();
    gc_ptrs = parent_g.get_init_gc_ptrs();
  }
  flint_cleanup_all_threads();
  std::cout << "After garbling: " << std::endl;
  raise_sigusr1();
  srand(rand_seed + 1);
  MainGadget parent_e(T, EVAL, args...);

  parent_e.set_init_gc_ptrs(gc_ptrs);

  srand(rand_seed + 2);
  // std::cout << "begin eval" << std::endl;
  SIMDWord::start_workers_e(num_workers);
#ifdef MAX_GADGET_TIME
  uint64_t eval_T = std::min(T, (uint64_t)MAX_GADGET_TIME / 2);
#else
  uint64_t eval_T = std::min(T, 100UL);
#endif
  auto start = std::chrono::high_resolution_clock::now();
  for (uint64_t i = 0; i < eval_T; ++i) {
    // std::cout << "time: " << i << " / " << T << std::endl;
    parent_e.main();
  }
  auto end = std::chrono::high_resolution_clock::now();
  SIMDWord::stop_workers_e();
  int64_t eval_time =
      std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
          .count();
  std::cout << "eval time: " << eval_time << " ms" << std::endl;
  std::cout << "avg time per time step: " << eval_time / eval_T << " ms"
            << std::endl;
  std::cout << "parallel time: " << ParTracker::get_time() * 1e3 << " ms"
            << std::endl;
  ParTracker::reset();
  test_timer.print_laps();
  if (test_timer.get_time() > 0) {
    std::cout << "test timer time: " << test_timer.get_time() * 1e3 << " ms"
              << std::endl;
    test_timer.reset();
  }
  close_file(fid);
#ifdef PRINT_PROGRESS_GRANULARITY
  next_gc_ckpt = 0;
#endif
}

template <typename MainGadget, typename... Args>
void test_gadget(TRangeType T_range_type = ANY,
                 std::pair<uint32_t, uint32_t> T_range = {1, 1000},
                 uint64_t buffer_size = 1024 * 1024 * 512, Args... args) {
  test_gadget_with_workers<MainGadget>(T_range_type, T_range, buffer_size, 0,
                                       args...);
}

template <typename MainGadget, typename... Args>
void test_gadget_dbg(TRangeType T_range_type = ANY,
                     std::pair<uint32_t, uint32_t> T_range = {1, 1000},
                     Args... args) {
  uint64_t rand_seed = time(0);
  // std::cout << "rand seed: " << rand_seed << std::endl;
  srand(rand_seed);
  uint64_t T = get_rand_T(T_range_type, T_range);
  srand(rand_seed + 1);
  MainGadget parent(T, DEBUG, args...);
  srand(rand_seed + 2);
  for (uint64_t i = 0; i < T; ++i) {
    // std::cout << "time: " << i << " / " << T << std::endl;
    parent.main();
  }
}

template <typename MainGadget, typename... Args>
void measure_gadget_gc(uint64_t T, Args... args) {
  uint64_t rand_seed = time(0);
  // srand(rand_seed);
  std::cout << "T: " << T << std::endl;
  srand(rand_seed + 1);
#ifdef MEASURE_STACK_COST
  global_stack_cost = 0;
#endif
#ifdef MEASURE_TSC_STACK
  global_tsc_stack_cost = 0;
#endif
  MainGadget parent(T, MEASURE, args...);
  GCPtr gc_ptr(-1);
  uint64_t gc_size = parent.garble(gc_ptr).get_offset();
  if (gc_size > 1UL << 30) {
    print_size_with_unit(gc_size, "GB", 1UL << 30);
#ifdef MEASURE_STACK_COST
    print_size_with_unit(global_stack_cost, "GB (stack cost)", 1UL << 30);
#endif
#ifdef MEASURE_TSC_STACK
    print_size_with_unit(global_tsc_stack_cost, "GB (tsc stack cost)",
                         1UL << 30);
#endif
  } else if (gc_size > 1UL << 20) {
    print_size_with_unit(gc_size, "MB", 1UL << 20);
#ifdef MEASURE_STACK_COST
    print_size_with_unit(global_stack_cost, "MB (stack cost)", 1UL << 20);
#endif
#ifdef MEASURE_TSC_STACK
    print_size_with_unit(global_tsc_stack_cost, "MB (tsc stack cost)",
                         1UL << 20);
#endif
  } else {
    print_size_with_unit(gc_size, "KB", 1UL << 10);
#ifdef MEASURE_STACK_COST
    print_size_with_unit(global_stack_cost, "KB (stack cost)", 1UL << 10);
#endif
#ifdef MEASURE_TSC_STACK
    print_size_with_unit(global_tsc_stack_cost, "KB (tsc stack cost)",
                         1UL << 10);
#endif
  }
}
}  // namespace ZebraGRAM