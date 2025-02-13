#include "hash.hpp"

#include <omp.h>
#define USE_EMP_TOOL 1
#if USE_EMP_TOOL
#include <emp-tool/utils/ccrh.h>
#endif
namespace PicoGRAM {

#if USE_EMP_TOOL
emp::CCRH ccrh;
#else
// Fixed 128-bit key for AES (16 bytes)
const unsigned char fixed_key[LAMBDA_BYTES] = {
    0x2b, 0x7e, 0x15, 0x16, 0x28, 0xae, 0xd2, 0xa6,
    0xab, 0xf7, 0x43, 0x34, 0x4f, 0x6d, 0x3b, 0x5d};
const unsigned char tweak[LAMBDA_BYTES] = {0x23, 0x44, 0x75, 0x34, 0x34, 0x25,
                                           0xa2, 0x54, 0x12, 0x45, 0x68, 0x23,
                                           0x45, 0x23, 0x45, 0x23};
AES_KEY encryptKey;
struct AES_KEY_initializer {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
  AES_KEY_initializer() { AES_set_encrypt_key(fixed_key, 128, &encryptKey); }
#pragma GCC diagnostic pop
};

AES_KEY_initializer aes_key_initializer;

// Function to XOR two 128-bit blocks (16 bytes)
void xor_128bit(const unsigned char* a, const unsigned char* b,
                unsigned char* result) {
#pragma omp simd
  for (int i = 0; i < LAMBDA_BYTES; ++i) {
    result[i] = a[i] ^ b[i];
  }
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
// Function to perform AES encryption with the fixed key
void aes_encrypt(const unsigned char* input, unsigned char* output) {
  AES_ecb_encrypt(input, output, &encryptKey, AES_ENCRYPT);
}
#pragma GCC diagnostic pop
// TCCR Hash implementation
void tccr_hash(const unsigned char* tweak, const __m128i& x, __m128i& hash) {
  unsigned char pi_x[LAMBDA_BYTES];             // π(x)
  unsigned char intermediate[LAMBDA_BYTES];     // π(x) ⊕ i
  unsigned char pi_intermediate[LAMBDA_BYTES];  // π(π(x) ⊕ i)

  // load x into pi_x
  _mm_storeu_si128((__m128i*)pi_x, x);

  // Step 1: Compute π(x) = AES(fixed_key, x)
  aes_encrypt(pi_x, pi_x);

  // Step 2: Compute π(x) ⊕ i
  xor_128bit(pi_x, tweak, intermediate);

  // Step 3: Compute π(π(x) ⊕ i) = AES(fixed_key, intermediate)
  aes_encrypt(intermediate, pi_intermediate);

  // Step 4: Compute H(i, x) = π(π(x) ⊕ i) ⊕ π(x)
  xor_128bit(pi_intermediate, pi_x, intermediate);

  // store the result in hash
  hash = _mm_loadu_si128((__m128i*)intermediate);
}
#endif

Label Hash(const Label& label) {
  Label result;
// tccr_hash(tweak, label.bits, result.bits);
#if USE_EMP_TOOL
  result.set_m128i(ccrh.H(label.get_m128i()));
#else
  __m128i result_m128i;
  tccr_hash(tweak, label.get_m128i(), result_m128i);
  result.set_m128i(result_m128i);
#endif
  return result;
}

Label add_nonce(const Label& label, uint64_t nonce1, uint64_t nonce2 = 0) {
  Label label_with_nonce;
#if LAMBDA_BYTES == 16
  // set nonce1 concatenated with nonce2 as mask
  __m128i mask = _mm_set_epi64x(nonce2, nonce1);
  label_with_nonce.set_m128i(_mm_xor_si128(label.get_m128i(), mask));

#else
  label_with_nonce = label;
  for (uint64_t i = 0; i < sizeof(uint64_t); ++i) {
    label_with_nonce.bits[i] ^= (nonce1 >> (8 * i)) & 0xFF;
  }
  for (uint64_t i = 0; i < sizeof(uint64_t); ++i) {
    label_with_nonce.bits[i + sizeof(uint64_t)] ^= (nonce2 >> (8 * i)) & 0xFF;
  }
#endif
  return label_with_nonce;
}

Label Hash(const Label& label, uint64_t nonce1, uint64_t nonce2) {
  Label label_with_nonce = add_nonce(label, nonce1, nonce2);
  return Hash(label_with_nonce);
}

Label Hash(const Label& label, uint64_t nonce) { return Hash(label, nonce, 0); }

void prng(const Label& seed, uint8_t* output, size_t length) {
  Label state = seed;
  for (size_t i = 0; i < length; i += LAMBDA_BYTES) {
    state = Hash(state);
    memcpy(output + i, state.get_ptr(),
           std::min((size_t)LAMBDA_BYTES, length - i));
  }
}

BigInt HashZ(const Label& label) {
  uint8_t bytes[2 * LAMBDA_BYTES];
  prng(label, bytes, 2 * LAMBDA_BYTES);
  return BigInt(bytes);
}

BigInt HashZ(const Label& label, uint64_t nonce1, uint64_t nonce2) {
  Label label_with_nonce = add_nonce(label, nonce1, nonce2);
  return HashZ(label_with_nonce);
}

Label Hash(const ECPoint& point) {
  return Hash(Label(point.bytes)) ^ Hash(Label(point.bytes + LAMBDA_BYTES));
}

MAC Mac(const Label& label) {
  Label result_extended = Hash(label);
  return MAC(result_extended);
}

void xor_bytes(uint8_t* a, const uint8_t* b, size_t len) {
  // a ^= b
  // TODO: optimize?
#pragma omp simd
  for (size_t i = 0; i < len; ++i) {
    a[i] ^= b[i];
  }
}

void enc_big_int(BigInt& plaintext, const Label& key,
                 uint8_t ciphertext[BigInt::byte_length], uint64_t nonce1,
                 uint64_t nonce2) {
  Label key_with_nonce = add_nonce(key, nonce1, nonce2);
  prng(key_with_nonce, ciphertext, BigInt::byte_length);
  uint8_t plaintext_bytes[BigInt::byte_length];
  plaintext.to_bytes(plaintext_bytes);
  xor_bytes(ciphertext, plaintext_bytes, BigInt::byte_length);
}

void dec_big_int(BigInt& plaintext, const Label& key,
                 const uint8_t ciphertext[BigInt::byte_length], uint64_t nonce1,
                 uint64_t nonce2) {
  Label key_with_nonce = add_nonce(key, nonce1, nonce2);
  uint8_t plaintext_bytes[BigInt::byte_length];
  prng(key_with_nonce, plaintext_bytes, BigInt::byte_length);
  xor_bytes(plaintext_bytes, ciphertext, BigInt::byte_length);
  plaintext.from_bytes(plaintext_bytes);
}

void enc_ec_point(const ECPoint& point, const Label& key,
                  uint8_t ciphertext[ECPoint::byte_length], uint64_t nonce1,
                  uint64_t nonce2) {
  Label key_with_nonce = add_nonce(key, nonce1, nonce2);
  prng(key_with_nonce, ciphertext, ECPoint::byte_length);
  xor_bytes(ciphertext, point.bytes, ECPoint::byte_length);
}

void dec_ec_point(ECPoint& point, const Label& key,
                  const uint8_t ciphertext[ECPoint::byte_length],
                  uint64_t nonce1, uint64_t nonce2) {
  point.is_temp_point_fresh = false;
  Label key_with_nonce = add_nonce(key, nonce1, nonce2);
  prng(key_with_nonce, point.bytes, ECPoint::byte_length);
  xor_bytes(point.bytes, ciphertext, ECPoint::byte_length);
}

}  // namespace PicoGRAM