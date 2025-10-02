#pragma once
#include <cstring>

#include "fmpz_helpers.hpp"
#include "global.hpp"
#include "paillier.hpp"
#include "rand.hpp"
// Helper functions for adding/subtracting powers of 2
// These functions are not in FLINT but are commonly needed

static inline void fmpz_add_2exp(fmpz_t result, const fmpz_t x, ulong exp) {
  fmpz_t power_of_two;
  fmpz_init(power_of_two);
  fmpz_one_2exp(power_of_two, exp);   // power_of_two = 2^exp
  fmpz_add(result, x, power_of_two);  // result = x + 2^exp
  fmpz_clear(power_of_two);
}

static inline void fmpz_sub_2exp(fmpz_t result, const fmpz_t x, ulong exp) {
  fmpz_t power_of_two;
  fmpz_init(power_of_two);
  fmpz_one_2exp(power_of_two, exp);   // power_of_two = 2^exp
  fmpz_sub(result, x, power_of_two);  // result = x - 2^exp
  fmpz_clear(power_of_two);
}

namespace ZebraGRAM {
struct ArithLabel {
  fmpz_t auth;
  fmpz_t raw;

  // maximum length in bytes
  static constexpr uint byte_length_auth = (LEN_AUTH_SHARE + 63) / 64 * 8;
  static constexpr uint byte_length_raw = (LEN_RAW_SHARE + 63) / 64 * 8;
  static constexpr uint byte_length = byte_length_auth + byte_length_raw;

  ArithLabel() {
    fmpz_init(auth);
    fmpz_init(raw);
  }
  ~ArithLabel() {
    fmpz_clear(auth);
    fmpz_clear(raw);
  }

  // copy constructor
  ArithLabel(const ArithLabel &other) {
    fmpz_init(auth);
    fmpz_init(raw);
    fmpz_set(auth, other.auth);
    fmpz_set(raw, other.raw);
  }

  // copy assignment operator
  ArithLabel &operator=(const ArithLabel &other) {
    if (this != &other) {
      fmpz_set(auth, other.auth);
      fmpz_set(raw, other.raw);
    }
    return *this;
  }

  // move constructor
  ArithLabel(ArithLabel &&other) noexcept {
    fmpz_init(auth);
    fmpz_init(raw);
    fmpz_swap(auth, other.auth);
    fmpz_swap(raw, other.raw);
  }

  // move assignment operator
  ArithLabel &operator=(ArithLabel &&other) noexcept {
    if (this != &other) {
      fmpz_swap(auth, other.auth);
      fmpz_swap(raw, other.raw);
    }
    return *this;
  }

  void merge(fmpz_t result) const {
    // result = auth << (LEN_STAT_PARAM + LEN_RAW_SHARE) + raw
    fmpz_mul_2exp(result, auth, LEN_RAW_SHARE);
    fmpz_add(result, result, raw);
  }

  void split(const fmpz_t input) {
    // auth = input >> (LEN_STAT_PARAM + LEN_RAW_SHARE)
    fmpz_fdiv_q_2exp(auth, input, LEN_RAW_SHARE);
    // raw = input mod 2^{LEN_STAT_PARAM + LEN_RAW_SHARE}
    fmpz_fdiv_r_2exp(raw, input, LEN_RAW_SHARE);
  }

  uint from_bytes(const uint8_t *input) {
    fmpz_set_ui_array(auth, (ulong *)input, byte_length_auth / 8);
    fmpz_set_ui_array(raw, (ulong *)(input + byte_length_auth),
                      byte_length_raw / 8);
    return byte_length;
  }

  uint to_bytes(uint8_t *out) const {
    fmpz_get_ui_array((ulong *)out, byte_length_auth / 8, auth);
    fmpz_get_ui_array((ulong *)(out + byte_length_auth), byte_length_raw / 8,
                      raw);
    return byte_length;
  }

  // overload addition and subtraction
  // for the arith part the addition is mod 2^LEN_AUTH_SHARE
  // for the raw part the addition is mod 2^LEN_RAW_SHARE
  ArithLabel operator+(const ArithLabel &other) const {
    ArithLabel result;
    fmpz_add(result.auth, auth, other.auth);
    fmpz_fdiv_r_2exp(result.auth, result.auth, LEN_AUTH_SHARE);
    fmpz_add(result.raw, raw, other.raw);
    fmpz_fdiv_r_2exp(result.raw, result.raw, LEN_RAW_SHARE);
    return result;
  }

  ArithLabel operator-(const ArithLabel &other) const {
    ArithLabel result;
    fmpz_sub(result.auth, auth, other.auth);
    // Handle negative results: if result < 0, add 2^LEN_AUTH_SHARE
    if (fmpz_sgn(result.auth) < 0) {
      fmpz_add_2exp(result.auth, result.auth, LEN_AUTH_SHARE);
    }

    fmpz_sub(result.raw, raw, other.raw);
    // Handle negative results: if result < 0, add 2^LEN_RAW_SHARE
    if (fmpz_sgn(result.raw) < 0) {
      fmpz_add_2exp(result.raw, result.raw, LEN_RAW_SHARE);
    }
    return result;
  }

  ArithLabel &operator+=(const ArithLabel &other) {
    fmpz_add(auth, auth, other.auth);
    fmpz_fdiv_r_2exp(auth, auth, LEN_AUTH_SHARE);
    fmpz_add(raw, raw, other.raw);
    fmpz_fdiv_r_2exp(raw, raw, LEN_RAW_SHARE);
    return *this;
  }

  ArithLabel &operator-=(const ArithLabel &other) {
    fmpz_sub(auth, auth, other.auth);
    // Handle negative results: if auth < 0, add 2^LEN_AUTH_SHARE
    if (fmpz_sgn(auth) < 0) {
      fmpz_add_2exp(auth, auth, LEN_AUTH_SHARE);
    }

    fmpz_sub(raw, raw, other.raw);
    // Handle negative results: if raw < 0, add 2^LEN_RAW_SHARE
    if (fmpz_sgn(raw) < 0) {
      fmpz_add_2exp(raw, raw, LEN_RAW_SHARE);
    }
    return *this;
  }

  // override negation
  ArithLabel operator-() const {
    ArithLabel result;
    fmpz_neg(result.auth, auth);
    // Handle negative results: if result < 0, add 2^LEN_AUTH_SHARE
    if (fmpz_sgn(result.auth) < 0) {
      fmpz_add_2exp(result.auth, result.auth, LEN_AUTH_SHARE);
    }

    fmpz_neg(result.raw, raw);
    // Handle negative results: if result < 0, add 2^LEN_RAW_SHARE
    if (fmpz_sgn(result.raw) < 0) {
      fmpz_add_2exp(result.raw, result.raw, LEN_RAW_SHARE);
    }
    return result;
  }

  // overload output stream
  friend std::ostream &operator<<(std::ostream &os, const ArithLabel &label) {
    char *auth_str = fmpz_get_str(NULL, 10, label.auth);
    char *raw_str = fmpz_get_str(NULL, 10, label.raw);
    os << "ArithLabel(auth=" << auth_str << ", raw=" << raw_str << ")";
    flint_free(auth_str);
    flint_free(raw_str);
    return os;
  }

  bool check() const {
    if (fmpz_sgn(auth) < 0 || fmpz_sizeinbase(auth, 2) > LEN_AUTH_SHARE) {
      return false;
    }
    if (fmpz_sgn(raw) < 0 || fmpz_sizeinbase(raw, 2) > LEN_RAW_SHARE) {
      return false;
    }
    return true;
  }
};

void secret_share(ArithLabel &share1, ArithLabel &share2, const fmpz_t secret,
                  const PaillierPrivKey &sk);

void secret_share(ArithLabel &share1, ArithLabel &share2, const fmpz_t secret,
                  size_t max_bit, const PaillierPrivKey &sk);

void open_share(fmpz_t result, const ArithLabel &share1,
                const ArithLabel &share2);
}  // namespace ZebraGRAM