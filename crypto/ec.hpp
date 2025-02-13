#pragma once
#include <openssl/bn.h>
#include <openssl/ec.h>
#include <openssl/obj_mac.h>

#include <iostream>

#include "label.hpp"
#include "rand.hpp"
#include "util.hpp"

namespace PicoGRAM {

/**
 * @brief A finite field element optimized for modular multiplication and
 * division with 256-bit prime modulus. Most operations are not thread-safe.
 *
 */
struct BigInt {
  static BIGNUM* q;  // The global prime modulus
  // global cache for converting operands to montgomery form
  static BIGNUM* mont_cache;
  // global context for BN operations
  static BN_CTX* global_ctx;
  // global montgomery context
  static BN_MONT_CTX* global_mont_ctx;
  // 256-bit prime modulus
  static constexpr int byte_length = 32;
  // A pool of BIGNUM objects
  static ResourcePool<BIGNUM> big_num_pool;

 private:
  BIGNUM* value = NULL;     // value either in standard or montgomery form
  bool montgomery = false;  // whether the value is in montgomery form

 public:
  // Default Constructor, defer initialization of value
  BigInt() {}

  // Copy constructor
  BigInt(const BigInt& other) : montgomery(other.montgomery) {
    if (other.is_set()) {
      value = big_num_pool.get();
      BN_copy(value, other.value);
    }
  }

  // Move constructor
  BigInt(BigInt&& other) : value(other.value), montgomery(other.montgomery) {
    other.value = nullptr;
  }

  /**
   * @brief Construct a new BigInt object from a byte array
   * @param bytes the byte array
   * If the number is larger than q, the first bit is clamped to 0
   */
  explicit BigInt(const uint8_t* bytes) {
    uint32_t first_4_bytes = *(uint32_t*)bytes;
    if (first_4_bytes != UINT32_MAX) {
      // number is smaller than q for P-256
      from_bytes(bytes);
    } else {
      // set first bit to 0
      uint8_t bytes_copy[byte_length];
      memcpy(bytes_copy, bytes, byte_length);
      bytes_copy[0] &= 0x7F;
      from_bytes(bytes_copy);
    }
  }

  /**
   * @brief Check if the value is set
   *
   */
  bool is_set() const { return value != nullptr; }

  void unset() {
    if (value) {
      big_num_pool.release(value);
    }
    value = nullptr;
    montgomery = false;
  }

  /**
   * @brief Convert this to montgomery form, accelerating modular multiplication
   *
   * @return BigInt& *this
   */
  BigInt& to_montgomery() {
    Assert(is_set());
    if (!montgomery) {
      Assert(!ParTracker::is_in_parallel());
      if (!BN_to_montgomery(value, value, global_mont_ctx, global_ctx)) {
        std::cerr << "Error converting to montgomery form!" << std::endl;
      }
      montgomery = true;
    }
    return *this;
  }

  /**
   * @brief Convert this from montgomery form to standard form.
   *
   * @return BigInt& *this
   */
  BigInt& from_montgomery() {
    Assert(is_set());
    if (montgomery) {
      Assert(!ParTracker::is_in_parallel());
      if (!BN_from_montgomery(value, value, global_mont_ctx, global_ctx)) {
        std::cerr << "Error converting from montgomery form!" << std::endl;
      }
      montgomery = false;
    }
    return *this;
  }

  bool is_montgomery() const { return montgomery; }

  /**
   * @brief Get the internal BIGNUM value (read-only)
   *
   * @return const BIGNUM* the internal BIGNUM value
   */
  const BIGNUM* get_value() const {
    Assert(is_set());
    Assert(!montgomery);
    return value;
  }

  // Destructor
  ~BigInt() { unset(); }

  // Overload assignment operator
  BigInt& operator=(const BigInt& other) {
    if (this != &other) {
      if (!other.is_set()) {
        unset();
      } else if (is_set()) {
        BN_copy(value, other.value);
      } else {
        value = BN_dup(other.value);
      }
    }
    montgomery = other.montgomery;
    return *this;
  }

  // Overload move assignment operator
  BigInt& operator=(BigInt&& other) {
    if (this != &other) {
      if (value) {
        big_num_pool.release(value);
      }
      value = other.value;
      other.value = nullptr;
      montgomery = other.montgomery;
    }
    return *this;
  }

  /**
   * @brief Overload multiplication (may change the value of this to montgomery
   * form). The result is in montgomery form. The operation is not thread-safe.
   *
   * @param other the other BigInt
   * @return BigInt the result of multiplication (in montgomery form)
   */
  BigInt operator*(const BigInt& other) const {
    BigInt result = *this;
    result *= other;
    return result;
  }

  /**
   * @brief Compute the modular inverse of this. The operation is
   * time-consuming. So try to minimize the use or call inv_batch instead.
   * The operation is not thread-safe.
   *
   * @return BigInt the modular inverse (in standard form)
   */
  BigInt inv() {
    Assert(is_set());
    BigInt result;
    result.value = big_num_pool.get();
    inv(result);
    return result;
  }

  /**
   * @brief Compute the modular inverse of this. The operation is
   * time-consuming. So try to minimize the use or call inv_batch instead.
   * The operation is not thread-safe.
   *
   * @param result the result of the modular inverse (in standard form)
   */
  void inv(BigInt& result) {
    from_montgomery();
    Assert(is_set());
    Assert(result.is_set());
    Assert(!ParTracker::is_in_parallel());
    if (BN_mod_inverse(result.value, value, q, global_ctx) == nullptr) {
      std::cerr << "Error computing modular inverse!" << std::endl;
    }
  }

  /**
   * @brief Compute the modular inverse of a batch of numbers more efficiently.
   * Only incurs one expensive modular inversion operation. The operation is not
   * thread-safe.
   *
   * @param nums the array of BigInt pointers (both input and output, the output
   * is in montgomery form)
   * @param count the number of BigInts
   */
  static void inv_batch(BigInt** nums, uint64_t count) {
    Assert(count > 0);
    std::vector<BigInt> prefix_products(count);
    prefix_products[0] = *nums[0];
    for (size_t i = 1; i < count; ++i) {
      prefix_products[i] = prefix_products[i - 1] * *nums[i];
    }
    BigInt& r = prefix_products.back();
    r.from_montgomery();
    r.inv(r);
    for (size_t i = count - 1; i > 0; --i) {
      BigInt new_num = r * prefix_products[i - 1];
      std::swap(*nums[i], new_num);
      r *= new_num;  // now new_num is nums[i]
    }
    *nums[0] = r;
  }

  /**
   * @brief Overload *= operator for multiplication. The result is in
   * montgomery. The operation is not thread-safe.
   *
   * @param other
   * @return BigInt&
   */
  BigInt& operator*=(const BigInt& other) {
    Assert(is_set());
    Assert(other.is_set());
    to_montgomery();
    Assert(!ParTracker::is_in_parallel());
    if (other.montgomery) {
      if (!BN_mod_mul_montgomery(value, value, other.value, global_mont_ctx,
                                 global_ctx)) {
        std::cerr << "Error performing montgomery multiplication!" << std::endl;
      }
    } else {
      // convert other value to montgomery form using mont_cache
      if (!BN_to_montgomery(mont_cache, other.value, global_mont_ctx,
                            global_ctx)) {
        std::cerr << "Error converting to montgomery form!" << std::endl;
      }

      if (!BN_mod_mul_montgomery(value, value, mont_cache, global_mont_ctx,
                                 global_ctx)) {
        std::cerr << "Error performing montgomery multiplication!" << std::endl;
      }
    }
    return *this;
  }

 private:
  /**
   * += when both are in or not in montgomery form
   */
  BigInt& inc_helper(const BigInt& other) {
    Assert(!ParTracker::is_in_parallel());
    Assert_eq(montgomery, other.montgomery);
    if (!BN_mod_add(value, value, other.value, q, global_ctx)) {
      std::cerr << "Error performing addition!" << std::endl;
    }
    return *this;
  }

 public:
  BigInt& operator+=(const BigInt& other) {
    Assert(is_set());
    Assert(other.is_set());
    if (montgomery == other.montgomery) {
      return inc_helper(other);
    }
    BigInt other_copy = other;
    if (other.montgomery) {
      other_copy.from_montgomery();
    } else {
      other_copy.to_montgomery();
    }
    return inc_helper(other_copy);
  }

  BigInt operator+(const BigInt& other) const {
    BigInt result = *this;
    result += other;
    return result;
  }

  /**
   * @brief Overload == operator for comparison. The operation is not
   * thread-safe.
   *
   * @param other
   * @return true
   * @return false
   */
  bool operator==(const BigInt& other) const {
    if (montgomery == other.montgomery) {
      return BN_cmp(value, other.value) == 0;
    } else {
      BigInt this_copy = *this;
      if (other.montgomery) {
        this_copy.to_montgomery();
      } else {
        this_copy.from_montgomery();
      }
      return this_copy == other;
    }
  }

  bool operator!=(const BigInt& other) const { return !(*this == other); }

  /**
   * @brief Generate a random BigInt. The operation is thread-safe, but may not
   * be efficient when run in parallel due to dynamic memory allocation.
   *
   * @return BigInt
   */
  static BigInt random() {
    BigInt result;
    result.value = big_num_pool.get();
    BN_rand_range(result.value, q);
    return result;
  }

  /**
   * @brief Generate a BigInt with value 1. The operation is thread-safe.
   *
   * @return BigInt
   */
  static BigInt unity() {
    BigInt result;
    result.value = big_num_pool.get();
    BN_one(result.value);
    return result;
  }

  static BigInt zero() {
    BigInt result;
    result.value = big_num_pool.get();
    BN_zero(result.value);
    return result;
  }

  // Overload << operator for printing
  friend std::ostream& operator<<(std::ostream& os, const BigInt& bigint) {
    if (!bigint.is_set()) {
      os << "N/A";
      return os;
    }
    os << (bigint.montgomery ? "Montgomery " : "");
    char* str = BN_bn2hex(bigint.value);
    os << str;
    free(str);
    return os;
  }

  /**
   * @brief Convert a byte array to a BigInt. The operation is thread-safe.
   *
   * @param bytes the input byte array
   * @return uint the number of bytes
   */
  uint from_bytes(const uint8_t bytes[byte_length]) {
    if (!is_set()) {
      value = big_num_pool.get();
    }
    BN_bin2bn(bytes, byte_length, value);
    montgomery = true;
    return byte_length;
  }

  /**
   * @brief Convert a BigInt to a byte array. The operation turns value into
   * standard form. The operation is thread-safe only if value is already in
   * standard form.
   *
   * @param bytes the output byte array
   * @return uint the number of bytes
   */
  uint to_bytes(uint8_t bytes[byte_length]) {
    Assert(is_set());
    to_montgomery();  // convert back to standard form
    int len = BN_bn2binpad(value, bytes, byte_length);
    Assert_eq(len, byte_length);
    return byte_length;
  }

  /**
   * @brief Encrypt the bigint using symmetric encryption. The operation turns
   * value into standard form. The operation is thread-safe only if value is
   * already in standard form.
   *
   * @param key the symmetric key
   * @param bytes the output byte array
   * @param nonce1 the first nonce (optional)
   * @param nonce2 the second nonce (optional)
   */
  void enc(const Label& key, uint8_t bytes[byte_length], uint64_t nonce1 = 0,
           uint64_t nonce2 = 0);

  /**
   * @brief Decrypt the bigint using symmetric encryption. The operation is
   * thread-safe.
   *
   * @param key the symmetric key
   * @param bytes the input byte array
   * @param nonce1 the first nonce (optional)
   * @param nonce2 the second nonce (optional)
   */
  void dec(const Label& key, const uint8_t bytes[byte_length],
           uint64_t nonce1 = 0, uint64_t nonce2 = 0);
};

/**
 * @brief A wrapper for OpenSSL's EC_POINT for curve P-256
 *
 */
struct ECPoint {
  static EC_GROUP* group;          // the curve group
  EC_POINT* temp_point = nullptr;  // a buffer for ec operations
  BN_CTX* ctx = nullptr;           // reusable bignum context
  static const uint byte_length = 33;
  uint8_t bytes[byte_length];        // actual data
  bool is_temp_point_fresh = false;  // whether temp_point is sync with bytes
  // a pool of ECPoints to reduce cost of initialization
  static ResourcePool<EC_POINT> ec_point_pool;
  static ResourcePool<BN_CTX> bn_ctx_pool;

  ECPoint() {}

  /**
   * @brief Copy constructor. The temp_point buffer is not copied.
   *
   * @param other the other ECPoint
   */
  ECPoint(const ECPoint& other) : ECPoint() {
    memcpy(bytes, other.bytes, byte_length);
  }

  /**
   * @brief Initialize the temp_point buffer
   *
   */
  void initialize_temp_point() {
    if (temp_point == nullptr) {
      temp_point = ec_point_pool.get_thread_safe();
      if (temp_point == nullptr) {
        std::cerr << "Error creating EC point" << std::endl;
      }
    }
    if (ctx == nullptr) {
      ctx = bn_ctx_pool.get_thread_safe();
      if (ctx == nullptr) {
        std::cerr << "Error creating BN context" << std::endl;
      }
    }
  }

  /**
   * @brief Move constructor
   *
   * @param other the other ECPoint
   */
  ECPoint(ECPoint&& other) : temp_point(other.temp_point), ctx(other.ctx) {
    other.temp_point = nullptr;
    other.ctx = nullptr;
    memcpy(bytes, other.bytes, byte_length);
    is_temp_point_fresh = other.is_temp_point_fresh;
  }

  /**
   * @brief Set the value of the ECPoint from bytes
   *
   */
  void sync_temp_point_from_bytes() {
    Assert(!is_temp_point_fresh);
    is_temp_point_fresh = true;
    initialize_temp_point();
    EC_POINT_oct2point(group, temp_point, bytes, byte_length, ctx);
  }

  /**
   * @brief Set the value of the bytes from the temp_point
   *
   */
  void sync_bytes_from_temp_point() {
    Assert(temp_point != nullptr);
    Assert(ctx);
    is_temp_point_fresh = true;
    size_t len =
        EC_POINT_point2oct(group, temp_point, POINT_CONVERSION_COMPRESSED,
                           bytes, byte_length, ctx);
    Assert_eq(len, byte_length);
  }

  /**
   * @brief Destructor
   *
   */
  ~ECPoint() {
    if (temp_point) {
      // EC_POINT_free(temp_point);
      ec_point_pool.release_thread_safe(temp_point);
    }
    temp_point = nullptr;
    if (ctx) {
      // BN_CTX_free(ctx);
      bn_ctx_pool.release_thread_safe(ctx);
    }
    ctx = nullptr;
  }

  /**
   * @brief Compute scalar multiplication of the ECPoint. The operation is
   * thread-safe.
   *
   * @param exp the scalar, must be in standard form
   * @return ECPoint the result of the scalar multiplication
   */
  ECPoint pow(const BigInt& exp) {
    ECPoint result;
    pow(exp, result);
    return result;
  }

  /**
   * @brief Compute scalar multiplication of the ECPoint. The operation is
   * thread-safe.
   *
   * @param exp the scalar, must be in standard form
   * @param result the result of the scalar multiplication
   */
  void pow(const BigInt& exp, ECPoint& result) {
    Assert(exp.is_set());
    if (!is_temp_point_fresh) {
      sync_temp_point_from_bytes();
    }
    result.initialize_temp_point();
    if (!EC_POINT_mul(group, result.temp_point, nullptr, temp_point,
                      exp.get_value(), ctx)) {
      std::cerr << "Error performing point multiplication for point"
                << std::endl;
    }
    result.sync_bytes_from_temp_point();
  }
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
  /**
   * @brief Compute scalar multiplication of multiple ECPoints. The operation is
   * thread-safe.
   *
   * @param exps the scalars, must be in standard form
   * @param points the ECPoints
   * @return ECPoint the result of the scalar multiplication
   */
  static void aggr_pow(std::vector<BigInt>& exps, std::vector<ECPoint>& points,
                       ECPoint& result) {
    size_t len = exps.size();
    Assert(len == points.size());
    for (size_t i = 0; i < len; ++i) {
      Assert(exps[i].is_set());
      Assert(!exps[i].is_montgomery());
      points[i].sync_temp_point_from_bytes();
    }
    result.initialize_temp_point();
    static constexpr size_t threshold_len = 32;

    if (len <= threshold_len) {
      // avoid dynamic memory allocation in parallel
      const BIGNUM* exps_bn[threshold_len];
      const EC_POINT* points_ec[threshold_len];
      for (size_t i = 0; i < len; ++i) {
        exps_bn[i] = exps[i].get_value();
        points_ec[i] = points[i].temp_point;
      }
      if (!EC_POINTs_mul(group, result.temp_point, nullptr, len, points_ec,
                         exps_bn, result.ctx)) {
        std::cerr << "Error performing point aggregation" << std::endl;
      }
    } else {
      std::vector<const BIGNUM*> exps_bn(exps.size());
      std::vector<const EC_POINT*> points_ec(exps.size());
      for (size_t i = 0; i < exps.size(); ++i) {
        exps_bn[i] = exps[i].get_value();
      }
      for (size_t i = 0; i < points.size(); ++i) {
        points_ec[i] = points[i].temp_point;
      }
      if (!EC_POINTs_mul(group, result.temp_point, nullptr, len,
                         points_ec.data(), exps_bn.data(), result.ctx)) {
        std::cerr << "Error performing point aggregation" << std::endl;
      }
    }
    result.sync_bytes_from_temp_point();
  }
#pragma GCC diagnostic pop

  /**
   * @brief Compute scalar multiplication of the generator. The operation is
   * thread-safe.
   *
   * @param exp the scalar, must be in standard form
   * @return ECPoint the result of the scalar multiplication
   */
  static ECPoint generator_pow(const BigInt& exp) {
    ECPoint result;
    generator_pow(exp, result);
    return result;
  }

  /**
   * @brief Compute scalar multiplication of the generator. The operation is
   * thread-safe.
   *
   * @param exp the scalar, must be in standard form
   * @param result the result of the scalar multiplication
   */
  static void generator_pow(const BigInt& exp, ECPoint& result) {
    Assert(!exp.is_montgomery());
    result.initialize_temp_point();
    if (!EC_POINT_mul(group, result.temp_point, exp.get_value(), nullptr,
                      nullptr, result.ctx)) {
      std::cerr << "Error performing point multiplication for generator"
                << std::endl;
    }
    result.sync_bytes_from_temp_point();
  }

  /**
   * @brief Overload == operator.
   *
   * @param other
   */
  bool operator==(const ECPoint& other) const {
    // return EC_POINT_cmp(group, point, other.point, nullptr) == 0;
    return memcmp(bytes, other.bytes, byte_length) == 0;
  }

  /**
   * @brief Overload != operator.
   *
   * @param other
   */
  bool operator!=(const ECPoint& other) const { return !(*this == other); }

  /**
   * @brief Encrypt the ECPoint using symmetric encryption. The operation is
   * thread-safe.
   *
   * @param key the symmetric key
   * @param ciphertext the output byte array
   * @param nonce1 the first nonce (optional)
   * @param nonce2 the second nonce (optional)
   */
  void enc(const Label& key, uint8_t ciphertext[byte_length],
           uint64_t nonce1 = 0, uint64_t nonce2 = 0) const;

  /**
   * @brief Decrypt the ECPoint using symmetric encryption. The operation is
   * thread-safe.
   *
   * @param key the symmetric key
   * @param ciphertext the input byte array
   * @param nonce1 the first nonce (optional)
   * @param nonce2 the second nonce (optional)
   */
  void dec(const Label& key, const uint8_t ciphertext[byte_length],
           uint64_t nonce1 = 0, uint64_t nonce2 = 0);

  // Overload << operator for printing
  friend std::ostream& operator<<(std::ostream& os, const ECPoint& point) {
    // print hex
    for (uint i = 0; i < byte_length; ++i) {
      os << std::hex << (int)point.bytes[i];
    }
    // reset to dec
    os << std::dec;
    return os;
  }
};

/**
 * @brief A struct for initializing the elliptic curve and bigint parameters as
 * well as global buffers
 *
 */
struct ECInitializer {
  ECInitializer();
  ~ECInitializer();
};
}  // namespace PicoGRAM