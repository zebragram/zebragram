#include "ec.hpp"

#include "hash.hpp"

namespace PicoGRAM {

BIGNUM* BigInt::q;
BIGNUM* BigInt::mont_cache;
BN_CTX* BigInt::global_ctx;
BN_MONT_CTX* BigInt::global_mont_ctx;
EC_GROUP* ECPoint::group = nullptr;
ResourcePool<BIGNUM> BigInt::big_num_pool(BN_new, BN_free);

EC_POINT* new_ec_point() { return EC_POINT_new(ECPoint::group); }
ResourcePool<EC_POINT> ECPoint::ec_point_pool(new_ec_point, EC_POINT_free);
ResourcePool<BN_CTX> ECPoint::bn_ctx_pool(BN_CTX_new, BN_CTX_free);

ECInitializer::ECInitializer() {
  ECPoint::group = EC_GROUP_new_by_curve_name(NID_X9_62_prime256v1);
  if (ECPoint::group == nullptr) {
    std::cerr << "Error creating EC group for P-256" << std::endl;
  }

  BigInt::q = BN_new();
  if (!EC_GROUP_get_order(ECPoint::group, BigInt::q, nullptr)) {
    std::cerr << "Error retrieving curve prime" << std::endl;
  }
  BigInt::mont_cache = BN_new();
  BigInt::global_ctx = BN_CTX_new();
  BigInt::global_mont_ctx = BN_MONT_CTX_new();
  if (!BN_MONT_CTX_set(BigInt::global_mont_ctx, BigInt::q,
                       BigInt::global_ctx)) {
    std::cerr << "Error setting up montgomery context" << std::endl;
  }
}

ECInitializer::~ECInitializer() {
  EC_GROUP_free(ECPoint::group);
  BN_free(BigInt::q);
  BN_free(BigInt::mont_cache);
  BN_CTX_free(BigInt::global_ctx);
  BN_MONT_CTX_free(BigInt::global_mont_ctx);
}
struct ECInitializer ECInitializer_instance;

void BigInt::enc(const Label& key, uint8_t bytes[BigInt::byte_length],
                 uint64_t nonce1, uint64_t nonce2) {
  enc_big_int(*this, key, bytes, nonce1, nonce2);
}

void BigInt::dec(const Label& key, const uint8_t bytes[BigInt::byte_length],
                 uint64_t nonce1, uint64_t nonce2) {
  dec_big_int(*this, key, bytes, nonce1, nonce2);
}

void ECPoint::enc(const Label& key, uint8_t ciphertext[byte_length],
                  uint64_t nonce1, uint64_t nonce2) const {
  enc_ec_point(*this, key, ciphertext, nonce1, nonce2);
}

void ECPoint::dec(const Label& key, const uint8_t ciphertext[byte_length],
                  uint64_t nonce1, uint64_t nonce2) {
  dec_ec_point(*this, key, ciphertext, nonce1, nonce2);
}

}  // namespace PicoGRAM