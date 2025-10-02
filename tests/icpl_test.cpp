#include <flint/flint.h>
#include <flint/fmpz.h>
#include <gtest/gtest.h>

#include <chrono>
#include <iomanip>
#include <iostream>
#include <ipcl/ipcl.hpp>
#include <ipcl/mod_exp.hpp>
#include <random>
#include <vector>

#include "fmpz_helpers.hpp"
#include "paillier.hpp"

class ModExpBenchmark {
 public:
  ModExpBenchmark() {}

  // Generate random FLINT integer with specified bit length
  void generateRandomFmpz(fmpz_t result, int bit_length) {
    secure_random_fmpz(result, bit_length);
  }

  // Generate odd random number for modulus
  void generateRandomOddFmpz(fmpz_t result, int bit_length) {
    generateRandomFmpz(result, bit_length);
    if (fmpz_is_even(result)) {
      fmpz_add_ui(result, result, 1);
    }
  }

  void runBenchmark(int vector_size, int bit_length, int num_iterations = 10) {
    std::cout << "\n=== Benchmark: ippModExp ===" << std::endl;
    std::cout << "Vector size: " << vector_size << std::endl;
    std::cout << "Bit length: " << bit_length << std::endl;
    std::cout << "Iterations: " << num_iterations << std::endl;
    std::cout << "=============================" << std::endl;

    // Initialize IPCL context
    ipcl::initializeContext("QAT");

    // Generate test data using FLINT
    std::vector<fmpz_t> base_fmpz(vector_size);
    std::vector<fmpz_t> exp_fmpz(vector_size);
    std::vector<fmpz_t> mod_fmpz(vector_size);

    // Initialize FLINT integers
    for (int i = 0; i < vector_size; ++i) {
      fmpz_init(base_fmpz[i]);
      fmpz_init(exp_fmpz[i]);
      fmpz_init(mod_fmpz[i]);
    }

    // Generate random test data
    std::cout << "Generating random test data..." << std::endl;
    for (int i = 0; i < vector_size; ++i) {
      generateRandomFmpz(base_fmpz[i], bit_length);
      generateRandomFmpz(exp_fmpz[i], bit_length / 2);  // Smaller exponent
      generateRandomOddFmpz(mod_fmpz[i], bit_length);

      // Ensure base < mod
      fmpz_mod(base_fmpz[i], base_fmpz[i], mod_fmpz[i]);
    }

    // Convert to BigNumber (conversion time not measured)
    std::cout << "Converting to BigNumber format..." << std::endl;
    std::vector<BigNumber> base_bn =
        flint_ipcl::fmpzVectorToBigNumberVector(base_fmpz);
    std::vector<BigNumber> exp_bn =
        flint_ipcl::fmpzVectorToBigNumberVector(exp_fmpz);
    std::vector<BigNumber> mod_bn =
        flint_ipcl::fmpzVectorToBigNumberVector(mod_fmpz);

    // Run benchmark
    std::vector<double> times;
    times.reserve(num_iterations);

    std::cout << "Running benchmark..." << std::endl;
    for (int iter = 0; iter < num_iterations; ++iter) {
      auto start = std::chrono::high_resolution_clock::now();

      // This is what we're benchmarking
      std::vector<BigNumber> result = ipcl::modExp(base_bn, exp_bn, mod_bn);

      auto end = std::chrono::high_resolution_clock::now();

      auto duration =
          std::chrono::duration_cast<std::chrono::microseconds>(end - start)
              .count();
      times.push_back(duration / 1000.0);  // Convert to milliseconds

      std::cout << "Iteration " << (iter + 1) << ": " << std::fixed
                << std::setprecision(3) << times.back() << " ms" << std::endl;
    }

    // Calculate statistics
    double total_time = 0.0;
    double min_time = times[0];
    double max_time = times[0];

    for (double time : times) {
      total_time += time;
      min_time = std::min(min_time, time);
      max_time = std::max(max_time, time);
    }

    double avg_time = total_time / num_iterations;

    // Calculate standard deviation
    double variance = 0.0;
    for (double time : times) {
      variance += (time - avg_time) * (time - avg_time);
    }
    double std_dev = std::sqrt(variance / num_iterations);

    // Print results
    std::cout << "\n=== Results ===" << std::endl;
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "Average time: " << avg_time << " ms" << std::endl;
    std::cout << "Min time: " << min_time << " ms" << std::endl;
    std::cout << "Max time: " << max_time << " ms" << std::endl;
    std::cout << "Std deviation: " << std_dev << " ms" << std::endl;
    std::cout << "Throughput: " << std::setprecision(1)
              << (vector_size / (avg_time / 1000.0)) << " ops/sec" << std::endl;
    std::cout << "Time per operation: " << std::setprecision(6)
              << (avg_time / vector_size) << " ms/op" << std::endl;

    // Verify correctness with FLINT (optional, for first few elements)
    if (vector_size > 0) {
      std::cout << "\n=== Verification ===" << std::endl;
      std::vector<BigNumber> ipcl_result =
          ipcl::ippModExp(base_bn, exp_bn, mod_bn);

      // Verify first element
      fmpz_t flint_result;
      fmpz_init(flint_result);
      fmpz_powm(flint_result, base_fmpz[0], exp_fmpz[0], mod_fmpz[0]);

      // convert ipcl_result[0] to fmpz_t
      fmpz_t ipcl_fmpz;
      fmpz_init(ipcl_fmpz);
      flint_ipcl::bigNumberToFmpz(ipcl_result[0], ipcl_fmpz);
      if (fmpz_cmp(flint_result, ipcl_fmpz) == 0) {
        std::cout << "First element verification passed." << std::endl;
      } else {
        std::cout << "First element verification FAILED!" << std::endl;
      }
      fmpz_clear(ipcl_fmpz);
      fmpz_clear(flint_result);
    }

    // Clean up FLINT integers
    for (int i = 0; i < vector_size; ++i) {
      fmpz_clear(base_fmpz[i]);
      fmpz_clear(exp_fmpz[i]);
      fmpz_clear(mod_fmpz[i]);
    }
  }

  // New method for speed comparison between IPCL and FLINT
  void runSpeedComparison(int vector_size, int bit_length,
                          int num_iterations = 10) {
    std::cout << "\n=== Speed Comparison: IPCL vs FLINT ===" << std::endl;
    std::cout << "Vector size: " << vector_size << std::endl;
    std::cout << "Bit length: " << bit_length << std::endl;
    std::cout << "Iterations: " << num_iterations << std::endl;
    std::cout << "========================================" << std::endl;

    // Initialize IPCL context
    ipcl::initializeContext("QAT");

    // Generate test data using FLINT
    std::vector<fmpz_t> base_fmpz(vector_size);
    std::vector<fmpz_t> exp_fmpz(vector_size);
    std::vector<fmpz_t> mod_fmpz(vector_size);

    // Initialize FLINT integers
    for (int i = 0; i < vector_size; ++i) {
      fmpz_init(base_fmpz[i]);
      fmpz_init(exp_fmpz[i]);
      fmpz_init(mod_fmpz[i]);
    }

    // Generate random test data
    for (int i = 0; i < vector_size; ++i) {
      generateRandomFmpz(base_fmpz[i], bit_length);
      generateRandomFmpz(exp_fmpz[i], bit_length / 2);  // Smaller exponent
      generateRandomOddFmpz(mod_fmpz[i], bit_length);

      // Ensure base < mod
      fmpz_mod(base_fmpz[i], base_fmpz[i], mod_fmpz[i]);
    }

    // Convert to BigNumber for IPCL
    std::vector<BigNumber> base_bn =
        flint_ipcl::fmpzVectorToBigNumberVector(base_fmpz);
    std::vector<BigNumber> exp_bn =
        flint_ipcl::fmpzVectorToBigNumberVector(exp_fmpz);
    std::vector<BigNumber> mod_bn =
        flint_ipcl::fmpzVectorToBigNumberVector(mod_fmpz);

    // Benchmark IPCL
    std::vector<double> ipcl_times;
    ipcl_times.reserve(num_iterations);

    std::cout << "\nBenchmarking IPCL ippModExp..." << std::endl;
    for (int iter = 0; iter < num_iterations; ++iter) {
      auto start = std::chrono::high_resolution_clock::now();

      std::vector<BigNumber> ipcl_result =
          ipcl::ippModExp(base_bn, exp_bn, mod_bn);

      auto end = std::chrono::high_resolution_clock::now();

      auto duration =
          std::chrono::duration_cast<std::chrono::microseconds>(end - start)
              .count();
      ipcl_times.push_back(duration / 1000.0);  // Convert to milliseconds
    }

    // Benchmark FLINT
    std::vector<double> flint_times;
    flint_times.reserve(num_iterations);

    std::cout << "Benchmarking FLINT fmpz_powm..." << std::endl;
    for (int iter = 0; iter < num_iterations; ++iter) {
      auto start = std::chrono::high_resolution_clock::now();

      // Compute modular exponentiation for all elements using FLINT
      std::vector<fmpz_t> flint_results(vector_size);
      for (int i = 0; i < vector_size; ++i) {
        fmpz_init(flint_results[i]);
        fmpz_powm(flint_results[i], base_fmpz[i], exp_fmpz[i], mod_fmpz[i]);
      }

      auto end = std::chrono::high_resolution_clock::now();

      auto duration =
          std::chrono::duration_cast<std::chrono::microseconds>(end - start)
              .count();
      flint_times.push_back(duration / 1000.0);  // Convert to milliseconds

      // Clean up results
      for (int i = 0; i < vector_size; ++i) {
        fmpz_clear(flint_results[i]);
      }
    }

    // Calculate statistics for IPCL
    double ipcl_total = 0.0;
    for (double time : ipcl_times) {
      ipcl_total += time;
    }
    double ipcl_avg = ipcl_total / num_iterations;

    // Calculate statistics for FLINT
    double flint_total = 0.0;
    for (double time : flint_times) {
      flint_total += time;
    }
    double flint_avg = flint_total / num_iterations;

    // Print comparison results
    std::cout << "\n=== Performance Comparison ===" << std::endl;
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "IPCL average time:  " << ipcl_avg << " ms" << std::endl;
    std::cout << "FLINT average time: " << flint_avg << " ms" << std::endl;
    std::cout << "IPCL time per op:   " << std::setprecision(6)
              << (ipcl_avg / vector_size) << " ms/op" << std::endl;
    std::cout << "FLINT time per op:  " << std::setprecision(6)
              << (flint_avg / vector_size) << " ms/op" << std::endl;

    double speedup = flint_avg / ipcl_avg;
    std::cout << std::setprecision(2);
    if (speedup > 1.0) {
      std::cout << "IPCL is " << speedup << "x faster than FLINT" << std::endl;
    } else {
      std::cout << "FLINT is " << (1.0 / speedup) << "x faster than IPCL"
                << std::endl;
    }

    std::cout << "IPCL throughput:    " << std::setprecision(1)
              << (vector_size / (ipcl_avg / 1000.0)) << " ops/sec" << std::endl;
    std::cout << "FLINT throughput:   " << std::setprecision(1)
              << (vector_size / (flint_avg / 1000.0)) << " ops/sec"
              << std::endl;

    // Verify correctness - compare first element
    std::vector<BigNumber> ipcl_result =
        ipcl::ippModExp(base_bn, exp_bn, mod_bn);

    fmpz_t flint_result;
    fmpz_init(flint_result);
    fmpz_powm(flint_result, base_fmpz[0], exp_fmpz[0], mod_fmpz[0]);

    fmpz_t ipcl_fmpz;
    fmpz_init(ipcl_fmpz);
    flint_ipcl::bigNumberToFmpz(ipcl_result[0], ipcl_fmpz);

    bool results_match = (fmpz_cmp(flint_result, ipcl_fmpz) == 0);
    std::cout << "\nCorrectness check: " << (results_match ? "PASS" : "FAIL")
              << std::endl;

    fmpz_clear(ipcl_fmpz);
    fmpz_clear(flint_result);

    // Clean up FLINT integers
    for (int i = 0; i < vector_size; ++i) {
      fmpz_clear(base_fmpz[i]);
      fmpz_clear(exp_fmpz[i]);
      fmpz_clear(mod_fmpz[i]);
    }
  }
};

TEST(TestModExpBenchmark, BasicRun) {
  // Parse command line arguments
  int vector_size = 100;
  int bit_length = 4096;
  int iterations = 10;

  std::cout << "IPCL ippModExp Benchmark" << std::endl;

  ModExpBenchmark benchmark;

  // Run main benchmark
  benchmark.runBenchmark(vector_size, bit_length, iterations);
}

TEST(TestModExpBenchmark, SpeedComparison) {
  std::cout << "IPCL vs FLINT Speed Comparison" << std::endl;

  ModExpBenchmark benchmark;

  // Test with different vector sizes and bit lengths
  std::vector<std::pair<int, int>> test_cases = {
      {10, 1024},   // Small vectors, moderate precision
      {50, 1024},   // Medium vectors, moderate precision
      {100, 2048},  // Large vectors, high precision
      {25, 4096}    // Small vectors, very high precision
  };

  for (const auto& test_case : test_cases) {
    int vector_size = test_case.first;
    int bit_length = test_case.second;

    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "Test case: " << vector_size << " elements, " << bit_length
              << " bits" << std::endl;

    benchmark.runSpeedComparison(vector_size, bit_length, 5);
  }
}