#pragma once

#include <openssl/aes.h>

#include <cstring>
#include <iostream>

#include "ec.hpp"
#include "label.hpp"

namespace ZebraGRAM {

/**
 * @brief Compute the hash of a label using fixed-key AES
 *
 * @param label the label to hash
 * @return Label the hash of the label
 */
Label Hash(const Label& label);

/**
 * @brief Compute the hash of a label with a nonce using fixed-key AES
 *
 * @param label the label to hash
 * @param nonce
 * @return Label
 */
Label Hash(const Label& label, uint64_t nonce);

/**
 * @brief Compute the hash of a label with two nonces using fixed-key AES
 *
 * @param label
 * @param nonce1
 * @param nonce2
 * @return Label
 */
Label Hash(const Label& label, uint64_t nonce1, uint64_t nonce2);

/**
 * @brief Compute the hash of an EC point using fixed-key AES
 *
 * @param point the EC point to hash
 * @return Label the hash of the EC point
 */
Label Hash(const ECPoint& point);

/**
 * @brief Compute the hash of a label using fixed-key AES and return a BigInt
 *
 * @param label the label to hash
 * @return BigInt the hash of the label as a 256-bit integer
 */
BigInt HashZ(const Label& label);

/**
 * @brief Compute the hash of a label with a nonce using fixed-key AES and
 * return a BigInt
 *
 * @param label the label to hash
 * @param nonce the nonce
 * @return BigInt the hash of the label as a 256-bit integer
 */
BigInt HashZ(const Label& label, uint64_t nonce1, uint64_t nonce2 = 0);

/**
 * @brief Encrypt a BigInt using symmetric encryption
 *
 * @param plaintext the BigInt to encrypt, the operation turns the BigInt into
 * standard form. The operation is thread-safe only if the BigInt is already in
 * standard form.
 * @param key the symmetric key
 * @param ciphertext the output byte array
 * @param nonce1 the first nonce (optional)
 * @param nonce2 the second nonce (optional)
 */
void enc_big_int(BigInt& plaintext, const Label& key,
                 uint8_t ciphertext[BigInt::byte_length], uint64_t nonce1 = 0,
                 uint64_t nonce2 = 0);

/**
 * @brief Decrypt a BigInt using symmetric encryption
 *
 * @param plaintext the BigInt to decrypt to, the output is in standard form
 * @param key the symmetric key
 * @param ciphertext the input byte array
 * @param nonce1 the first nonce (optional)
 * @param nonce2 the second nonce (optional)
 */
void dec_big_int(BigInt& plaintext, const Label& key,
                 const uint8_t ciphertext[BigInt::byte_length],
                 uint64_t nonce1 = 0, uint64_t nonce2 = 0);

/**
 * @brief Encrypt an ECPoint using symmetric encryption. The operation is
 * thread-safe.
 *
 * @param point the EC point to encrypt
 * @param key the symmetric key
 * @param ciphertext the output byte array
 * @param nonce1 the first nonce (optional)
 * @param nonce2 the second nonce (optional)
 */
void enc_ec_point(const ECPoint& point, const Label& key,
                  uint8_t ciphertext[ECPoint::byte_length], uint64_t nonce1 = 0,
                  uint64_t nonce2 = 0);

/**
 * @brief Decrypt an ECPoint using symmetric encryption. The operation is
 * thread-safe.
 *
 * @param point the EC point to decrypt to
 * @param key the symmetric key
 * @param ciphertext the input byte array
 * @param nonce1 the first nonce (optional)
 * @param nonce2 the second nonce (optional)
 */
void dec_ec_point(ECPoint& point, const Label& key,
                  const uint8_t ciphertext[ECPoint::byte_length],
                  uint64_t nonce1 = 0, uint64_t nonce2 = 0);

/**
 * @brief Compute a sigma-bit hash of a label using fixed-key AES
 *
 * @param label the label to hash
 * @return MAC
 */
MAC Mac(const Label& label);

/**
 * @brief Pseudorandom number generator using AES in counter mode
 *
 * @param seed the seed label
 * @param output the output byte array
 * @param length the length of the output byte array
 */
void prng(const Label& seed, uint8_t* output, size_t length);
}  // namespace ZebraGRAM