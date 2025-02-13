#pragma once

#include <atomic>
#include <mutex>
#include <thread>
#if defined(__x86_64__) || defined(_M_X64) || defined(_M_IX86)
#include <immintrin.h>
#define SPIN_WAIT() _mm_pause()
#else
#define SPIN_WAIT() std::this_thread::yield()
#endif

struct BusyThread {
  std::thread thread;
  std::atomic<bool> stopped{true};
  std::atomic<bool> has_job{false};

  void start() {
    if (!stopped.load(std::memory_order_acquire)) {
      return;
    }
    stopped.store(false, std::memory_order_release);
    thread = std::thread(&BusyThread::exec, this);
  }

  void wait() {
    while (has_job.load(std::memory_order_acquire)) {
      // busy waiting
    }
  }

  void stop() {
    if (stopped.load(std::memory_order_acquire)) {
      return;
    }
    stopped.store(true, std::memory_order_release);
    has_job.store(true, std::memory_order_release);
    if (thread.joinable()) {
      thread.join();
    }
  }

  void exec() {
    while (true) {
      while (!has_job.load(std::memory_order_acquire)) {
        SPIN_WAIT();
        // busy waiting
      }
      if (stopped.load(std::memory_order_acquire)) {
        has_job.store(false, std::memory_order_release);
        break;
      }
      exec_step();
      has_job.store(false, std::memory_order_release);
    }
  }

  virtual void exec_step() = 0;
};

// TTASSpinlock
struct Lock {
 public:
  Lock() : flag(false) {}
  Lock(const Lock&) {}

  void lock() {
    // First, test in a non-atomic way
    while (flag.load(std::memory_order_relaxed)) {
      // Spin without attempting to set the flag
      // Optionally, use std::this_thread::yield() to reduce contention
    }
    // Now attempt to set the flag atomically
    while (flag.exchange(true, std::memory_order_acquire)) {
      // If exchange returns true, the flag was already set, so keep spinning
    }
  }

  void unlock() { flag.store(false, std::memory_order_release); }

 private:
  std::atomic<bool> flag;
};

struct Mutex {
  void lock() { _lock.lock(); }
  void unlock() { _lock.unlock(); }
  Mutex() {}
  Mutex(const Mutex&) {}

 private:
  std::mutex _lock;
};

template <typename LockType>
struct Critical {
  LockType& _lock;
  explicit Critical(LockType& _lock) : _lock(_lock) { _lock.lock(); }
  ~Critical() { _lock.unlock(); }
};