#pragma once
#include <cxxabi.h>
#include <execinfo.h>
#include <omp.h>
#include <stdlib.h>
#include <unistd.h>

#include <array>
#include <chrono>
#include <csignal>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <numeric>
#include <sstream>
#include <unordered_set>
#include <vector>

#include "global.hpp"
#include "lock.hpp"

// Assert
#ifndef NDEBUG
#define Assert(expr)                                                   \
  do {                                                                 \
    if (!(expr)) {                                                     \
      std::cerr << "Assertion failed: (" #expr << "), "                \
                << "file " << __FILE__ << ", line " << __LINE__ << "." \
                << std::endl;                                          \
      print_stack_trace();                                             \
      std::raise(SIGTRAP); /* Triggers a breakpoint */                 \
      std::abort();                                                    \
    }                                                                  \
  } while (false)

#define Assert_eq(expr1, expr2)                                      \
  do {                                                               \
    if ((expr1) != (expr2)) {                                        \
      std::stringstream ss;                                          \
      ss << "Assertion failed: (" #expr1 " == " #expr2 ") "          \
         << "evaluated as (" << (expr1) << " != " << (expr2) << ") " \
         << "in file " << __FILE__ << ", line " << __LINE__ << ".";  \
      std::cerr << ss.str() << std::endl;                            \
      print_stack_trace();                                           \
      std::raise(SIGTRAP); /* Triggers a breakpoint */               \
      std::abort();                                                  \
    }                                                                \
  } while (false)

#define Assert_less(expr1, expr2)                                              \
  do {                                                                         \
    if (!((expr1) < (expr2))) {                                                \
      std::stringstream ss;                                                    \
      ss << "Assertion failed: (" #expr1 " < " #expr2 ") " << "evaluated as (" \
         << (expr1) << " >= " << (expr2) << ") "                               \
         << "in file " << __FILE__ << ", line " << __LINE__ << ".";            \
      std::cerr << ss.str() << std::endl;                                      \
      print_stack_trace();                                                     \
      std::raise(SIGTRAP); /* Triggers a breakpoint */                         \
      std::abort();                                                            \
    }                                                                          \
  } while (false)
#else
#define Assert(expr) \
  do {               \
    (void)(expr);    \
  } while (0)
#define Assert_eq(expr1, expr2) \
  do {                          \
    (void)(expr1);              \
    (void)(expr2);              \
  } while (0)

#define Assert_less(expr1, expr2) \
  do {                            \
    (void)(expr1);                \
    (void)(expr2);                \
  } while (0)
#endif

#define STRR(X) #X
#define STR(X) STRR(X)

extern void raise_sigusr1();

extern void flint_cleanup_all_threads();

namespace ZebraGRAM {

uint log2(uint64_t n);
uint log2ceil(uint64_t n);
uint bit_width(uint64_t n);

double chernoff_upper_bound(double N, double M, double K,
                            int epsilon_exp = -60);

/**
 * @brief Reverse the count least significant bits of a 64-bit word
 *
 * @param word
 * @param count
 * @return uint64_t
 */
uint64_t reverse_bits(uint64_t word, uint count);

/**
 * @brief Generate an insecure 64 bit random number for testing
 *
 * @return uint64_t
 */
uint64_t rand64();

/**
 * @brief Clear and release the memory of a container
 *
 * @tparam Container
 * @param vec
 */
template <typename Container>
void clear_and_release(Container& vec) {
  Container().swap(vec);
}
}  // namespace ZebraGRAM

void print_stack_trace();

// extern bool print_malloc_flag;

// void* operator new(size_t size);

/**
 * @brief A resource pool that manages the allocation and deallocation of
 * repeatedly used objects.
 *
 * @tparam T the type of the object
 */
template <typename T>
struct ResourcePool {
 private:
  std::vector<T*> pool;
#ifndef NDEBUG
  // for checking if the pointer is already in the pool
  std::unordered_set<T*> pool_set;
#endif
  Lock lock;
  // custom new that returns T*
  T* (*new_func)();

  // custom deleter
  void (*deleter)(T*);

 public:
  ResourcePool() = delete;

  /**
   * @brief Construct a new Resource Pool
   *
   * @param new_func the function to create a new object
   * @param deleter the function to delete an object
   */
  ResourcePool(T* (*new_func)(), void (*deleter)(T*))
      : new_func(new_func), deleter(deleter) {}

  ~ResourcePool() {
    clear();
    Assert(pool.empty());
  }

  /**
   * @brief Get an object from the pool. If the pool is empty, create a new
   * object. The caller is responsible for releasing the object back to the
   * pool. The function is not thread-safe.
   *
   * @return T*
   */
  T* get() {
    if (pool.empty()) {
      return new_func();
    }
    T* ptr = pool.back();
    pool.pop_back();
#ifndef NDEBUG
    Assert(pool_set.erase(ptr));
#endif
    return ptr;
  }

  /**
   * @brief Get an object from the pool. If the pool is empty, create a new
   * object. The caller is responsible for releasing the object back to the
   * pool. The function is thread-safe.
   *
   * @return T*
   */
  T* get_thread_safe() {
    lock.lock();
    T* ptr = get();
    lock.unlock();
    return ptr;
  }

  /**
   * @brief Release an object back to the pool. The function is not thread-safe.
   *
   * @param ptr
   */
  void release(T* ptr) {
    pool.push_back(ptr);
#ifndef NDEBUG
    Assert(pool_set.insert(ptr).second);
#endif
  }

  /**
   * @brief Release an object back to the pool. The function is thread-safe.
   *
   * @param ptr
   */
  void release_thread_safe(T* ptr) {
    lock.lock();
    release(ptr);
    lock.unlock();
  }

  /**
   * @brief Clear the pool and release the memory of the objects. The function
   * is thread-safe.
   *
   */
  void clear() {
    lock.lock();
    for (T*& ptr : pool) {
      deleter(ptr);
      ptr = nullptr;
    }
    pool.clear();
#ifndef NDEBUG
    pool_set.clear();
#endif
    lock.unlock();
  }
};

/**
 * @brief A timer that measures the accumulated time elapsed between tic and
 * toc.
 *
 */
struct Timer {
 private:
  double start;
  double end;
  double total_time;
  std::vector<double> tic_lap_times;
  std::vector<double> toc_lap_times;
#ifndef NDEBUG
  bool running = false;
#endif

 public:
  Timer() : start(0.0), end(0.0), total_time(0.0) {}

  /**
   * @brief Start the timer
   *
   */
  void tic() {
#ifndef NDEBUG
    Assert(!running);
    running = true;
#endif
    start = omp_get_wtime();
  }

  /**
   * @brief Start the timer and record the time elapsed since the last toc
   *
   * @param id the record id
   */
  void tic(uint16_t id) {
#ifndef NDEBUG
    Assert(!running);
    running = true;
#endif
    start = omp_get_wtime();
    if (id >= tic_lap_times.size()) {
      tic_lap_times.resize(id + 1, 0.0);
    }
    if (end != 0.0) {
      tic_lap_times[id] += start - end;
    }
  }

  /**
   * @brief Stop the timer and accumulate the time elapsed since the last tic to
   * the total time
   *
   */
  void toc() {
#ifndef NDEBUG
    Assert(running);
    running = false;
#endif
    end = omp_get_wtime();
    total_time += end - start;
  }

  /**
   * @brief Stop the timer and record the time elapsed since the last tic to
   * both the total time and the lap time for the record id
   *
   * @param id the record id
   */
  void toc(uint16_t id) {
#ifndef NDEBUG
    Assert(running);
    running = false;
#endif
    end = omp_get_wtime();
    if (id >= toc_lap_times.size()) {
      toc_lap_times.resize(id + 1, 0.0);
    }
    double lap_time = end - start;
    toc_lap_times[id] += lap_time;
    total_time += lap_time;
  }

  /**
   * @brief Reset the timer
   *
   */
  void reset() {
    total_time = 0.0;
    start = 0.0;
    end = 0.0;
    tic_lap_times.clear();
    toc_lap_times.clear();
#ifndef NDEBUG
    running = false;
#endif
  }

  /**
   * @brief Get the total time elapsed in seconds
   *
   * @return double
   */
  double get_time() const { return total_time; }

  /**
   * @brief Print the lap times for each id.
   *
   */
  void print_laps() const {
    for (uint16_t i = 0; i < tic_lap_times.size(); ++i) {
      if (tic_lap_times[i] > 0) {
        std::cout << "Tic Lap " << i << ": " << tic_lap_times[i] * 1000 << " ms"
                  << std::endl;
      }
    }
    for (uint16_t i = 0; i < toc_lap_times.size(); ++i) {
      if (toc_lap_times[i] > 0) {
        std::cout << "Toc Lap " << i << ": " << toc_lap_times[i] * 1000 << " ms"
                  << std::endl;
      }
    }
  }
};

// a global timer for testing
extern Timer test_timer;

/**
 * @brief Keep track of whether the program is in parallel execution and measure
 * the time elapsed in parallel execution.
 *
 */
class ParTracker {
 private:
  static Timer par_timer;
  static bool parallel_enabled;

 public:
  static void start() {
    Assert(!parallel_enabled);
    par_timer.tic();
    parallel_enabled = true;
  }

  static void stop() {
    Assert(parallel_enabled);
    par_timer.toc();
    parallel_enabled = false;
  }

  static double get_time() { return par_timer.get_time(); }

  static void reset() { par_timer.reset(); }

  static bool is_in_parallel() { return parallel_enabled; }
};

#ifdef MEASURE_STACK_COST
extern uint64_t global_stack_cost;
extern bool measure_stack_flag;
#endif

#ifdef MEASURE_TSC_STACK
extern uint64_t global_tsc_stack_cost;
#endif

#ifdef PRINT_PROGRESS_GRANULARITY
extern uint64_t next_gc_ckpt;
extern void print_curr_gc_offset(uint64_t curr_gc_offset);
#endif