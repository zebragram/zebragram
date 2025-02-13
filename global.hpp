#pragma once

// Computational security parameter
#define LAMBDA 128
#define LAMBDA_BYTES (LAMBDA / 8)

// Statistical security parameter, i.e. the length of the string in the
// Ungroup gate that lets the evaluator determine the correct row to decrypt
// while we target at a failure probability of 2^{-40}, we set the length to 64
// bits due to the union bound
#define SIGMA 64
#define SIGMA_BYTES (SIGMA / 8)

// set the default number of threads for the omp parallel for loop
// the NUM_THREADS here do not include the worker threads for computing elliptic
// curve multiplication
#define NUM_THREADS 1

// controls how the garbler and evaluator store the garbled circuit
#define STORAGE_TYPE SHARED_MEMORY

// enable this flag to use the Three-Halves-Make-A-Whole technique to garble AND
// gate
#define USE_THREE_HALVES 1

// enable this flag to enable the use of channels in EMP-Toolkit to transmit the
// garbled circuit
// #define USE_EMP_CHANNEL 1

#ifndef USE_EMP_CHANNEL

// enable this flag accelerates the measure mode
// however, correctness is guaranteed only for our use cases and its not compatible
// with the EMP-Toolkit
#define FAST_MEASURE 1

// enable this flag to measure the communication cost from stack in measure mode
#define MEASURE_STACK_COST 1

// if this flag is enabled, the measure mode will also simulate the case where
// we replace all the non-simd stack with the stack in the work of TSC #define
// #define MEASURE_TSC_STACK 1
#endif

namespace PicoGRAM {
/**
 * @brief The mode of the program
 *
 */
enum Mode {
  DEFAULT,  // unspecified mode
  GARBLE,   // garbler mode
  EVAL,     // evaluator mode
  MEASURE,  // measure mode, for performance evaluation
  DEBUG     // debug mode, skip crypto operations
};

enum LinkType { NONE, DIRECT, COND, SIMD_COND };

enum IOType {
  INVALID_IO,
  MEM_IO,
  NET_IO,
  FILE_IO,
  HIGH_SPEED_NET_IO,
};

// how garbled circuits are stored
enum StorageType {
  SHARED_MEMORY,  // for single process
  MEMORY,         // works for multiple processes that share file descriptors
  DISK,
};
}  // namespace PicoGRAM