#pragma once
#include <immintrin.h>
#include <omp.h>

#include <cstring>

#include "global.hpp"
#include "rand.hpp"
namespace PicoGRAM {

#if LAMBDA_BYTES == 16
struct Label {
 private:
  __m128i bits;

 public:
  Label() : bits(_mm_setzero_si128()) {}

  void* get_ptr() { return reinterpret_cast<void*>(&bits); }

  const void* get_ptr() const { return reinterpret_cast<const void*>(&bits); }

  explicit Label(const uint8_t* bytes)
      : bits(_mm_loadu_si128(reinterpret_cast<const __m128i*>(bytes))) {}

  explicit Label(const __m128i& bits) : bits(bits) {}

  explicit Label(const __m128i&& bits) : bits(bits) {}

  Label operator^(const Label& other) const {
    return Label(_mm_xor_si128(bits, other.bits));
  }

  Label& operator^=(const Label& other) {
    bits = _mm_xor_si128(bits, other.bits);
    return *this;
  }

  bool operator==(const Label& other) const {
    __m128i cmp = _mm_cmpeq_epi8(bits, other.bits);
    return _mm_movemask_epi8(cmp) == 0xFFFF;
  }

  bool operator!=(const Label& other) const { return !(*this == other); }

  /**
   * @brief Retrieve the last significant bit of the label.
   *
   * @return uint8_t either 0 or 1
   */
  uint8_t LSB() const { return bits[0] & 1; }

  /**
   * @brief Generate a random label.
   *
   * @return Label
   */
  static Label random() {
    uint8_t bytes[LAMBDA_BYTES];
    secure_random_128bit(bytes);
    return Label(_mm_loadu_si128(reinterpret_cast<__m128i*>(bytes)));
  }

  /**
   * @brief Generate a random label with last siginificant bit = 1.
   *
   * @return Label
   */
  static Label random_odd() {
    Label label = random();
    label.bits = _mm_or_si128(label.bits, _mm_set_epi64x(0, 1));
    return label;
  }

  // Overload << operator for printing
  friend std::ostream& operator<<(std::ostream& os, const Label& label) {
    uint8_t bytes[LAMBDA_BYTES];
    _mm_storeu_si128(reinterpret_cast<__m128i*>(bytes), label.bits);
    for (int i = 0; i < LAMBDA_BYTES; i++) {
      os << (uint)bytes[i] << " ";
    }
    return os;
  }

  /**
   * @brief Get the internal m128i object
   *
   * @return const __m128i&
   */
  const __m128i& get_m128i() const { return bits; }

  /**
   * @brief Set the internal m128i object
   *
   * @param value
   */
  void set_m128i(__m128i value) { bits = value; }

  /**
   * @brief Check if the least siginificant bit is 1
   *
   * @return true
   * @return false
   */
  bool is_odd() const { return _mm_extract_epi8(bits, 0) & 1; }
};
#else
/**
 * @brief A label for symmetric encryption
 *
 */
struct Label {
  uint8_t bits[LAMBDA_BYTES];

  Label() { memset(bits, 0, LAMBDA_BYTES); }
  explicit Label(const uint8_t* bytes) { memcpy(bits, bytes, LAMBDA_BYTES); }

  void* get_ptr() { return reinterpret_cast<void*>(bits); }

  const void* get_ptr() const { return reinterpret_cast<const void*>(bits); }

  Label operator^(const Label& other) const {
    Label result;
#pragma omp simd
    for (int i = 0; i < LAMBDA_BYTES; i++) {
      result.bits[i] = bits[i] ^ other.bits[i];
    }
    return result;
  }

  Label& operator^=(const Label& other) {
#pragma omp simd
    for (int i = 0; i < LAMBDA_BYTES; i++) {
      bits[i] ^= other.bits[i];
    }
    return *this;
  }

  bool operator==(const Label& other) const {
    for (int i = 0; i < LAMBDA_BYTES; i++) {
      if (bits[i] != other.bits[i]) {
        return false;
      }
    }
    return true;
  }

  bool operator!=(const Label& other) const { return !(*this == other); }

  uint8_t LSB() const { return bits[0] & 1; }

  static Label random() { return Label(secure_random_128bit(label.bits)); }

  static Label random_odd() {
    Label label = random();
    label.bits[0] |= 1;
    return label;
  }

  // Overload << operator for printing
  friend std::ostream& operator<<(std::ostream& os, const Label& label) {
    for (int i = 0; i < LAMBDA_BYTES; i++) {
      os << (int)label.bits[i] << " ";
    }
    return os;
  }
};
#endif

/**
 * @brief A message authentication code (MAC) to check if the correct row is
 * decrypted
 *
 */
struct MAC {
  uint8_t bits[SIGMA_BYTES];

  MAC() { memset(bits, 0, SIGMA_BYTES); }
  explicit MAC(const uint8_t* bytes) { memcpy(bits, bytes, SIGMA_BYTES); }
  explicit MAC(const Label& label) {
    // retrieve the first 5 bytes
    memcpy(bits, label.get_ptr(), SIGMA_BYTES);
  }

  bool operator==(const MAC& other) const {
    for (int i = 0; i < SIGMA_BYTES; i++) {
      if (bits[i] != other.bits[i]) {
        return false;
      }
    }
    return true;
  }

  bool operator!=(const MAC& other) const { return !(*this == other); }
};
}  // namespace PicoGRAM