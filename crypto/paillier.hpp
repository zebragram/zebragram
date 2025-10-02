#pragma once
#include <cstring>
#include <memory>

#include "powm.hpp"
#ifdef IPCL_AVAILABLE
#include "fmpz_helpers.hpp"
#endif

void sample_prime_prob(fmpz_t p, int bits);

// sample prime
void sample_prime(fmpz_t p, int bits);

#define LEN_PRIME 1024
#define LEN_N (2 * LEN_PRIME)
#define LEN_N2 (2 * LEN_N)
#define LEN_K 256
#define LEN_STAT_PARAM 60
#define LEN_RAW_SHARE ((LEN_N - LEN_K) / 2)
#define LEN_AUTH_SHARE (LEN_RAW_SHARE + LEN_K)
#define LEN_ARITH_DIGIT (LEN_RAW_SHARE - LEN_STAT_PARAM)

#define powm powm_cpu
#define powm_batch powm_cpu_batch
#define powm_precomputed powm_precomputed_cpu
#define powm_precomputed_batch powm_precomputed_cpu_batch

struct PaillierPubKey {
  fmpz_t N;                        // modulus p * q
  fmpz_t N2;                       // modulus squared
  fmpz_t g;                        // the hard group generator
  fmpz_t pk;                       // g^\K mod N^2 for damgard jurik elgamal
  fmpz_mod_ctx_struct mod_N2_ctx;  // context for mod N^2 operations
  fmpz_mod_ctx_struct mod_N_ctx;   // context for mod N operations

  PaillierPubKey() {
    fmpz_init(N);
    fmpz_init(N2);
    fmpz_init(g);
    fmpz_init(pk);
  }

  ~PaillierPubKey() {
    fmpz_clear(N);
    fmpz_clear(N2);
    fmpz_clear(g);
    fmpz_clear(pk);
    fmpz_mod_ctx_clear(&mod_N2_ctx);
    fmpz_mod_ctx_clear(&mod_N_ctx);
  }
};

struct PaillierPrivKey {
  fmpz_t g;       // the hard group generator
  fmpz_t p2, q2;  // primes squared
  fmpz_t p2inv;   // p2^{-1} mod q2
  fmpz_t K;       // LEN_K bit long secret key for damgard jurik elgamal
  fmpz_mod_ctx_struct mod_p2_ctx;       // context for mod p^2 operations
  fmpz_mod_ctx_struct mod_q2_ctx;       // context for mod q^2 operations
  fmpz_mod_ctx_struct mod_p_prime_ctx;  // context for mod (p-1)/2 operations
  fmpz_mod_ctx_struct mod_q_prime_ctx;  // context for mod (q-1)/2 operations
  PowmPrecomputeTable table_g_p2;  // precomputation table for g^{x} mod p^2
  PowmPrecomputeTable table_g_q2;  // precomputation table for g^{x} mod q^2
  PowmPrecomputeTable
      table_pk_p2;  // precomputation table for (g^-\K)^{x} mod p^2
  PowmPrecomputeTable
      table_pk_q2;  // precomputation table for (g^-\K)^{x} mod q^2
  PaillierPrivKey() {
    fmpz_init(g);
    fmpz_init(p2);
    fmpz_init(q2);
    fmpz_init(p2inv);
    fmpz_init(K);
  }

  ~PaillierPrivKey() {
    fmpz_clear(g);
    fmpz_clear(p2);
    fmpz_clear(q2);
    fmpz_clear(p2inv);
    fmpz_clear(K);
    fmpz_mod_ctx_clear(&mod_p2_ctx);
    fmpz_mod_ctx_clear(&mod_q2_ctx);
    fmpz_mod_ctx_clear(&mod_p_prime_ctx);
    fmpz_mod_ctx_clear(&mod_q_prime_ctx);
  }
};

extern PaillierPubKey default_pk;
extern PaillierPrivKey default_sk;
extern fmpz_t default_random_mask_raw;
extern fmpz_t default_random_mask_auth;

// Generate a Paillier keypair with modulus size of bits
void paillier_keygen(PaillierPubKey &pk, PaillierPrivKey &sk,
                     size_t precompute_mb = 16);

// 64-byte aligned
#define PAILLIER_CIPHER_TEXT_SIZE (((LEN_N2 + 63) / 64 * 8) * 2)

struct PaillierCiphertext {
  fmpz_t ct[2];
  PaillierCiphertext() {
    fmpz_init(ct[0]);
    fmpz_init(ct[1]);
  }
  ~PaillierCiphertext() {
    fmpz_clear(ct[0]);
    fmpz_clear(ct[1]);
  }
  // delete copy constructor and assignment
  PaillierCiphertext(const PaillierCiphertext &) = delete;
  PaillierCiphertext &operator=(const PaillierCiphertext &) = delete;

  void serialize(uint8_t out[PAILLIER_CIPHER_TEXT_SIZE]) {
    const size_t bytes_per_element = PAILLIER_CIPHER_TEXT_SIZE / 2;
    const size_t limbs_per_element =
        (bytes_per_element + sizeof(ulong) - 1) / sizeof(ulong);

    for (int i = 0; i < 2; i++) {
      ulong *limbs = (ulong *)malloc(limbs_per_element * sizeof(ulong));
      memset(limbs, 0, limbs_per_element * sizeof(ulong));
      fmpz_get_ui_array(limbs, limbs_per_element, ct[i]);
      memcpy(out + i * bytes_per_element, limbs, bytes_per_element);
      free(limbs);
    }
  }

  static void deserialize(fmpz_t ct0, fmpz_t ct1,
                          const uint8_t in[PAILLIER_CIPHER_TEXT_SIZE]) {
    const size_t bytes_per_element = PAILLIER_CIPHER_TEXT_SIZE / 2;
    const size_t limbs_per_element =
        (bytes_per_element + sizeof(ulong) - 1) / sizeof(ulong);

    for (int i = 0; i < 2; i++) {
      ulong *limbs = (ulong *)malloc(limbs_per_element * sizeof(ulong));
      memcpy(limbs, in + i * bytes_per_element, bytes_per_element);
      // Zero-pad any remaining limb bytes
      if (bytes_per_element < limbs_per_element * sizeof(ulong)) {
        memset((uint8_t *)limbs + bytes_per_element, 0,
               limbs_per_element * sizeof(ulong) - bytes_per_element);
      }
      if (i == 0)
        fmpz_set_ui_array(ct0, limbs, limbs_per_element);
      else
        fmpz_set_ui_array(ct1, limbs, limbs_per_element);
      free(limbs);
    }
  }

  void deserialize(const uint8_t in[PAILLIER_CIPHER_TEXT_SIZE]) {
    deserialize(ct[0], ct[1], in);
  }

  void print() const {
    printf("PaillierCiphertext:\n");
    printf("ct[0]: ");
    fmpz_print(ct[0]);
    printf("\nct[1]: ");
    fmpz_print(ct[1]);
    printf("\n");
  }
};

// CRT reconstruction: combines values mod p^2 and mod q^2 into value mod N^2
void crt_reconstruct(fmpz_t result, const fmpz_t val_p2, const fmpz_t val_q2,
                     const PaillierPrivKey &sk);

// Fast inverse mod p^2 using Hensel lifting
void fast_inv_mod_prime2(fmpz_t result, const fmpz_t a,
                         const fmpz_mod_ctx_struct *mod_p_ctx,
                         const fmpz_mod_ctx_struct *mod_p2_ctx);

struct PaillierEncryptCtx {
  fmpz_t *ct_p2[2], *ct_q2[2], *one_plus_mN;

  // size_t total_iter_1 = 0;
  // size_t total_iter_2 = 0;
  size_t batch_size = 0;

  void init(size_t batch_size) {
    // total_iter_1 = batch_size * 4;
    // total_iter_2 = batch_size;
    this->batch_size = batch_size;
    one_plus_mN = new fmpz_t[batch_size];
    for (int i = 0; i < 2; i++) {
      ct_p2[i] = new fmpz_t[batch_size];
      ct_q2[i] = new fmpz_t[batch_size];
      for (size_t j = 0; j < batch_size; j++) {
        fmpz_init(ct_p2[i][j]);
        fmpz_init(ct_q2[i][j]);
      }
    }
    for (size_t j = 0; j < batch_size; j++) {
      fmpz_init(one_plus_mN[j]);
    }
  }

  // void iter_1(size_t i, const PaillierPrivKey &sk, const fmpz_t *r);

  void push_tasks(const fmpz_t *r, const PaillierPrivKey &sk,
                  std::vector<const fmpz_t *> &all_exps,
                  std::vector<const PowmPrecomputeTable *> &all_tables);

  void retrieve_results(size_t i, PaillierCiphertext *ct, const fmpz_t *m,
                        const PaillierPubKey &pk, const PaillierPrivKey &sk,
                        const fmpz_t &result_p2_0, const fmpz_t &result_q2_0,
                        const fmpz_t &result_p2_1, const fmpz_t &result_q2_1);

  // void iter_2(size_t index, PaillierCiphertext *ct, const fmpz_t *m, const
  // PaillierPubKey &pk, const PaillierPrivKey &sk);

  ~PaillierEncryptCtx() {
    for (int i = 0; i < 2; i++) {
      for (size_t j = 0; j < batch_size; j++) {
        fmpz_clear(ct_p2[i][j]);
        fmpz_clear(ct_q2[i][j]);
      }
      delete[] ct_p2[i];
      delete[] ct_q2[i];
    }
    for (size_t j = 0; j < batch_size; j++) {
      fmpz_clear(one_plus_mN[j]);
    }
    delete[] one_plus_mN;
  }
};

void paillier_encrypt_with_sk(PaillierCiphertext &ct, const fmpz_t m,
                              const PaillierPubKey &pk,
                              const PaillierPrivKey &sk, const fmpz_t r);

void paillier_generator_pow(fmpz_t result, const fmpz_t exp,
                            const PaillierPrivKey &sk);

void paillier_encrypt_with_sk(PaillierCiphertext &ct, const fmpz_t m,
                              const PaillierPubKey &pk,
                              const PaillierPrivKey &sk);

void paillier_encrypt_with_sk_batch(PaillierCiphertext *ct, const fmpz_t *m,
                                    const PaillierPubKey &pk,
                                    const PaillierPrivKey &sk, const fmpz_t *r,
                                    size_t batch_size);

void paillier_decrypt(fmpz_t m, const PaillierCiphertext &ct,
                      const PaillierPubKey &pk, const PaillierPrivKey &sk);

void ddlog(fmpz_t result, const fmpz_t x, const PaillierPubKey &pk);
