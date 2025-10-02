#pragma once
#include "arith_label.hpp"
#include "paillier.hpp"
#define GARBLE_MULT_SIZE (2 * PAILLIER_CIPHER_TEXT_SIZE)
namespace ZebraGRAM {

void garble_mul(ArithLabel &z, uint8_t material[GARBLE_MULT_SIZE],
                const ArithLabel &x, const ArithLabel &y,
                const PaillierPubKey &pk, const PaillierPrivKey &sk);

// z_vec = x * y_vec
size_t garble_mul_vec(ArithLabel *z_vec, uint8_t *material, const ArithLabel &x,
                      const ArithLabel *y_vec, size_t vec_size,
                      const PaillierPubKey &pk, const PaillierPrivKey &sk);

void eval_mul(ArithLabel &z, const uint8_t material[GARBLE_MULT_SIZE],
              const ArithLabel &x, const ArithLabel &y,
              const PaillierPubKey &pk);

void eval_mul_vec(ArithLabel *z_vec, const uint8_t *material,
                  const ArithLabel &x, const ArithLabel *y_vec, size_t vec_size,
                  const PaillierPubKey &pk);

// Optimized version that collects all operations into a single large batch
void eval_mul_vec_batched(ArithLabel *z_vec, const uint8_t *material,
                          const ArithLabel &x, const ArithLabel *y_vec,
                          size_t vec_size, const PaillierPubKey &pk);

struct GarbleMulVecCtx {
  fmpz_t *r, *m;
  fmpz_t x_auth_minus_Kx_raw;
  fmpz_t *g_powers, *g_powers_mod_p_prime, *g_powers_mod_q_prime, *tmps;
  PaillierCiphertext *ct;
  PaillierEncryptCtx enc_ctx;
  size_t total_iter_1;
  size_t total_iter_2;
  // size_t total_iter_3;
  // size_t batch_size;

  void init_with_randomness(const ArithLabel &x, const ArithLabel *y_vec,
                            size_t batch_size, const PaillierPrivKey &sk,
                            const uint8_t *randomness) {
    this->total_iter_1 = batch_size + 1;
    this->total_iter_2 = batch_size;
    r = new fmpz_t[total_iter_1];
    m = new fmpz_t[total_iter_1];
    fmpz_init(x_auth_minus_Kx_raw);

    // Initialize arrays
    for (size_t i = 0; i < total_iter_1; i++) {
      fmpz_init(r[i]);
      fmpz_init(m[i]);
      fmpz_set_ui_array(r[i], (ulong *)(randomness + i * (LEN_K / 8)),
                        LEN_K / 64);
    }

    // Compute m[0] = (x.raw * K) << LEN_RAW_SHARE + x.raw (only depends on x)
    fmpz_mul(x_auth_minus_Kx_raw, x.raw, sk.K);
    fmpz_mul_2exp(m[0], x_auth_minus_Kx_raw, LEN_RAW_SHARE);
    fmpz_add(m[0], m[0], x.raw);

    // Precompute x_auth - K * x_raw (only depends on x)
    fmpz_sub(x_auth_minus_Kx_raw, x.auth, x_auth_minus_Kx_raw);

    // Compute m[1..vec_size] for each y_vec element
    for (size_t i = 0; i < batch_size; i++) {
      y_vec[i].merge(m[i + 1]);
    }

    // Encrypt all messages: m[0] once, then m[1..batch_size]
    ct = new PaillierCiphertext[total_iter_1];

    // Process each multiplication
    g_powers = new fmpz_t[batch_size];
    g_powers_mod_p_prime = new fmpz_t[batch_size];
    g_powers_mod_q_prime = new fmpz_t[batch_size];
    tmps = new fmpz_t[batch_size];
    for (size_t i = 0; i < batch_size; i++) {
      fmpz_init(g_powers[i]);
      fmpz_init(tmps[i]);
      fmpz_init(g_powers_mod_p_prime[i]);
      fmpz_init(g_powers_mod_q_prime[i]);
    }
    enc_ctx.init(total_iter_1);
  }

  void init(const ArithLabel &x, const ArithLabel *y_vec, size_t batch_size,
            const PaillierPrivKey &sk) {
    this->total_iter_1 = batch_size + 1;
    uint8_t *randomness = new uint8_t[total_iter_1 * (LEN_K / 8)];
    RAND_bytes(randomness, total_iter_1 * (LEN_K / 8));
    init_with_randomness(x, y_vec, batch_size, sk, randomness);
    delete[] randomness;
  }

  void push_iter_1_tasks(const PaillierPrivKey &sk,
                         std::vector<const fmpz_t *> &all_exps,
                         std::vector<const PowmPrecomputeTable *> &all_tables) {
    // Push exponents for iter_1
    enc_ctx.push_tasks(r, sk, all_exps, all_tables);
  }

  void retrieve_iter_1_results(size_t i, uint8_t *material,
                               const PaillierPubKey &pk,
                               const PaillierPrivKey &sk,
                               const fmpz_t &result_p2_0,
                               const fmpz_t &result_q2_0,
                               const fmpz_t &result_p2_1,
                               const fmpz_t &result_q2_1) {
    enc_ctx.retrieve_results(i, ct, m, pk, sk, result_p2_0, result_q2_0,
                             result_p2_1, result_q2_1);
    ct[i].serialize(material + i * PAILLIER_CIPHER_TEXT_SIZE);
  }

  void set_iter_2_tasks(size_t i, const ArithLabel *y_vec,
                        const PaillierPrivKey &sk, const fmpz_t *&exp1,
                        const fmpz_t *&exp2, const PowmPrecomputeTable *&table1,
                        const PowmPrecomputeTable *&table2) {
    // g's power: r[0] * (y_auth - K * y_raw) + r[i+1] * (x_auth - K * x_raw)
    fmpz_mul(tmps[i], sk.K, y_vec[i].raw);      // tmp = K * y_raw
    fmpz_sub(tmps[i], y_vec[i].auth, tmps[i]);  // tmps[i] = y_auth - K * y_raw
    fmpz_mul(g_powers[i], r[0],
             tmps[i]);  // g_powers[i] = r[0] * (y_auth - K * y_raw)

    fmpz_mul(tmps[i], r[i + 1],
             x_auth_minus_Kx_raw);  // tmps[i] = r[i+1] * (x_auth - K * x_raw)
    fmpz_add(g_powers[i], g_powers[i], tmps[i]);  // g_powers[i] += tmps[i]
    fmpz_mod_set_fmpz(g_powers_mod_p_prime[i], g_powers[i],
                      &sk.mod_p_prime_ctx);
    fmpz_mod_set_fmpz(g_powers_mod_q_prime[i], g_powers[i],
                      &sk.mod_q_prime_ctx);
    exp1 = &g_powers_mod_p_prime[i];
    exp2 = &g_powers_mod_q_prime[i];
    table1 = &sk.table_g_p2;
    table2 = &sk.table_g_q2;
  }

  // result is the output of paillier_generator_pow
  void retrieve_iter_2_results(size_t i, ArithLabel *z_vec,
                               const ArithLabel *y_vec,
                               const PaillierPubKey &pk,
                               const PaillierPrivKey &sk,
                               const fmpz_t &result_p2,
                               const fmpz_t &result_q2) {
    crt_reconstruct(tmps[i], result_p2, result_q2, sk);
    ddlog(tmps[i], tmps[i], pk);

    // z_combine = -(tmps[i] + m[0] * y_raw) mod N
    fmpz_mul(g_powers[i], m[0], y_vec[i].raw);  // reuse g_powers[i] as temp
    fmpz_add(tmps[i], tmps[i], g_powers[i]);
    fmpz_neg(tmps[i], tmps[i]);
    fmpz_mod_set_fmpz(tmps[i], tmps[i], &pk.mod_N_ctx);

    z_vec[i].split(tmps[i]);
  }

  // void run_iter_3(size_t i, uint8_t *material, const PaillierPubKey &pk,
  // const PaillierPrivKey &sk) {
  //     enc_ctx.iter_2(i, ct, m, pk, sk);
  //     ct[i].serialize(material + i* PAILLIER_CIPHER_TEXT_SIZE);
  // }

  ~GarbleMulVecCtx() {
    // Cleanup
    for (size_t i = 0; i < total_iter_2 + 1; i++) {
      fmpz_clear(r[i]);
      fmpz_clear(m[i]);
    }
    for (size_t i = 0; i < total_iter_2; i++) {
      fmpz_clear(g_powers[i]);
      fmpz_clear(g_powers_mod_p_prime[i]);
      fmpz_clear(g_powers_mod_q_prime[i]);
      fmpz_clear(tmps[i]);
    }
    delete[] r;
    delete[] m;
    delete[] ct;
    delete[] g_powers;
    delete[] g_powers_mod_p_prime;
    delete[] g_powers_mod_q_prime;
    delete[] tmps;

    fmpz_clear(x_auth_minus_Kx_raw);
  }
};
}  // namespace ZebraGRAM