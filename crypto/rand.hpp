#pragma once

#include <openssl/rand.h>

#include <iomanip>
#include <iostream>
#include <stdexcept>

#include "util.hpp"
namespace PicoGRAM {

/**
 * @brief Generate a 128 bit random number
 *
 * @param random_bytes the output byte buffer, must have length at least 16
 * bytes
 */
void secure_random_128bit(unsigned char* random_bytes);

/**
 * @brief Fill the output buffer with random bytes
 *
 * @param random_bytes the output buffer
 * @param len the output length
 */
void secure_random(uint8_t* random_bytes, size_t len);

/**
 * @brief Generate a random uint64_t
 *
 * @return uint64_t
 */
uint64_t secure_random_uint64();

uint8_t secure_random_byte();

// generate random bits in batch
struct RandBitPool {
  static const size_t pool_size_bytes = 1024;
  static const size_t pool_size_bits = pool_size_bytes * 8;
  uint8_t pool[pool_size_bytes];
  size_t pool_pos = pool_size_bits;

  RandBitPool() {}

  uint8_t get_rand_bits(uint8_t n) {
    Assert_less(n, 9);
    if (pool_pos + n > pool_size_bits) {
      secure_random(pool, pool_size_bytes);
      pool_pos = 0;
    }
    uint8_t result = 0;
    for (uint8_t i = 0; i < n; ++i) {
      uint8_t byte_pos = pool_pos / 8;
      uint8_t bit_pos = pool_pos % 8;
      result |= ((pool[byte_pos] >> bit_pos) & 1) << i;
      pool_pos++;
    }
    return result;
  }
};

extern RandBitPool rand_bit_pool;

/**
 * @brief Permute the elements in the range [begin, end) according to a random
 * number perm. Requires end - begin <= 8.
 *
 * @tparam Iterator
 * @param begin the begin iterator
 * @param end the end iterator
 * @param perm a 64 bit random number
 */
template <class Iterator>
void permute_small(Iterator begin, Iterator end, uint64_t perm) {
  Assert_less(end - begin, 9);
  for (Iterator i = end - 1; i > begin; --i) {
    uint64_t r_lim = i - begin + 1;
    uint64_t r = perm % r_lim;
    perm /= r_lim;
    std::iter_swap(i, begin + r % (i - begin + 1));
  }
}

/**
 * @brief Randomly permute the elements in the range [begin, end) using
 * Fisher-Yate shuffle. To reduce the overhead of generating randomness,
 * consider using permute_small if the size of the range is known to be small.
 *
 * @tparam Iterator
 * @param begin the begin iterator
 * @param end the end iterator
 */
template <class Iterator>
void secure_permute(Iterator begin, Iterator end) {
  // fisher-yates shuffle
  uint64_t size = end - begin;
  if (size > 8 || size <= 2) {
    for (Iterator i = end - 1; i > begin; --i) {
      std::iter_swap(i, begin + secure_random_uint64() % (i - begin + 1));
    }
  } else {
    uint64_t perm = secure_random_uint64();
    permute_small(begin, end, perm);
  }
}

}  // namespace PicoGRAM