#pragma once
#include "arith_label.hpp"
#include "global.hpp"
#include "key_manager.hpp"
#include "label.hpp"
#include "paillier.hpp"

namespace ZebraGRAM {
class Mod2Ctx {
  static constexpr uint sec_param = LAMBDA;
  uint8_t Delta[sec_param] = {0};  // the global difference for free-XOR
  PaillierCiphertext cts[sec_param];
  fmpz_t rands[sec_param];                    // garbler's random exponents
  fmpz_t rands_times_K[sec_param];            // rands[i] * K
  PowmPrecomputeTable *tables[sec_param][2];  // powers of ciphertexts of cts
 public:
  static Label bits_to_label(const uint8_t bits[sec_param]) {
    uint8_t bytes[LAMBDA_BYTES] = {0};
    for (size_t i = 0; i < sec_param; ++i) {
      bytes[i / 8] |= (bits[i] & 1) << (i % 8);
    }
    return Label(bytes);
  }

  void garble_init(const PaillierPrivKey &sk, const PaillierPubKey &pk) {
    // randomly sample a 128-bit string R
    fmpz_t Delta_fmpz[sec_param];
    Label Delta_label = key_manager.get_Delta();
    uint8_t *Delta_bytes = (uint8_t *)Delta_label.get_ptr();
    for (size_t i = 0; i < sec_param; ++i) {
      // extract bit i of Delta_label
      // Note: little-endian bit order to match bits_to_label
      Delta[i] = (Delta_bytes[i >> 3] >> (i & 7)) & 0x1;
      fmpz_init(Delta_fmpz[i]);
      fmpz_set_ui(Delta_fmpz[i], Delta[i]);
      // sample random rands
      fmpz_init(rands[i]);
      fmpz_init(rands_times_K[i]);
      secure_random_fmpz(rands[i], LEN_K);
      fmpz_mul(rands_times_K[i], rands[i], sk.K);
      // encrypt Delta[i] with sk
    }
    // encrypt with sk batch
    paillier_encrypt_with_sk_batch(cts, Delta_fmpz, pk, sk, rands, sec_param);
  }

  void eval_set_cts(const PaillierCiphertext cts_in[sec_param]) {
    for (size_t i = 0; i < sec_param; ++i) {
      fmpz_set(cts[i].ct[0], cts_in[i].ct[0]);
      fmpz_set(cts[i].ct[1], cts_in[i].ct[1]);
    }
  }

  void eval_precompute(const PaillierPubKey &pk,
                       const double precompute_mb = 1024.0) {
#pragma omp parallel for collapse(2)
    for (size_t i = 0; i < sec_param; ++i) {
      for (int b = 0; b < 2; ++b) {
        tables[i][b] = new PowmPrecomputeTable();
        if (b == 0) {
          tables[i][b]->init(cts[i].ct[b], pk.N2, LEN_AUTH_SHARE,
                             precompute_mb / sec_param / 2);
        } else {
          tables[i][b]->init(cts[i].ct[b], pk.N2, LEN_RAW_SHARE,
                             precompute_mb / sec_param / 2);
        }
      }
    }
  }

  Label garble_mod2(const ArithLabel &x, const PaillierPrivKey &sk,
                    const PaillierPubKey &pk) const {
    // g^(x.auth * rands[i] - x.raw * rands_times_K[i]) mod N^2
    fmpz_t exps[sec_param * 2];
    const PowmPrecomputeTable *exp_tables[sec_param * 2];
    fmpz_t results[sec_param * 2];
    for (size_t i = 0; i < sec_param; ++i) {
      fmpz_init(exps[2 * i]);
      fmpz_init(exps[2 * i + 1]);
      fmpz_t tmp1, tmp2;
      fmpz_init(tmp1);
      fmpz_init(tmp2);
      // x.auth * rands[i]
      fmpz_mul(tmp1, x.auth, rands[i]);
      // x.raw * rands_times_K[i]
      fmpz_mul(tmp2, x.raw, rands_times_K[i]);
      fmpz_sub(exps[2 * i], tmp1, tmp2);
      fmpz_mod_set_fmpz(exps[2 * i + 1], exps[2 * i], &sk.mod_q_prime_ctx);
      fmpz_mod_set_fmpz(exps[2 * i], exps[2 * i], &sk.mod_p_prime_ctx);
      fmpz_clear(tmp1);
      fmpz_clear(tmp2);
      fmpz_init(results[2 * i]);
      fmpz_init(results[2 * i + 1]);
      exp_tables[2 * i] = &sk.table_g_p2;
      exp_tables[2 * i + 1] = &sk.table_g_q2;
    }
    // batch powm_precomputed
    // std::cout << "before powm_precomputed_batch" << std::endl;
    powm_precomputed_batch(results, (const fmpz_t *)exps,
                           (const PowmPrecomputeTable **)exp_tables,
                           sec_param * 2);
    // std::cout << "after powm_precomputed_batch" << std::endl;
    // combine with CRT
    fmpz_t reconstructed[sec_param];
    for (size_t i = 0; i < sec_param; ++i) {
      fmpz_init(reconstructed[i]);
    }
    uint8_t y_bool[sec_param];
#pragma omp parallel for
    for (size_t i = 0; i < sec_param; ++i) {
      crt_reconstruct(reconstructed[i], results[2 * i], results[2 * i + 1], sk);
      if (Delta[i]) {
        // *= (1 + N * x.raw)
        fmpz_t one_plus_Nx;
        fmpz_init(one_plus_Nx);
        fmpz_mul(one_plus_Nx, x.raw, pk.N);
        fmpz_add_ui(one_plus_Nx, one_plus_Nx, 1);
        fmpz_mod_mul(reconstructed[i], reconstructed[i], one_plus_Nx,
                     &pk.mod_N2_ctx);
        fmpz_clear(one_plus_Nx);
      }  // not constant time
      // ddlog
      ddlog(reconstructed[i], reconstructed[i], pk);
      // set y_bool = reconstructed[i] mod 2
      y_bool[i] = fmpz_fdiv_ui(reconstructed[i], 2);
    }

    for (size_t i = 0; i < sec_param * 2; ++i) {
      fmpz_clear(exps[i]);
      fmpz_clear(results[i]);
    }
    for (size_t i = 0; i < sec_param; ++i) {
      fmpz_clear(reconstructed[i]);
    }
    Label y = bits_to_label(y_bool);
    return y;
  }

  // Batch version: computes mod 2 over a batch in a single
  // powm_precomputed_batch and a single omp parallel for
  void garble_mod2_batch(const ArithLabel *xs, size_t batch, Label *ys,
                         const PaillierPrivKey &sk,
                         const PaillierPubKey &pk) const {
    const size_t total = batch * sec_param;

    // Allocate exponent and result arrays for the entire batch
    fmpz_t *exps = new fmpz_t[2 * total];
    fmpz_t *results = new fmpz_t[2 * total];
    const PowmPrecomputeTable **exp_tables =
        new const PowmPrecomputeTable *[2 * total];

    // Build exponents and tables for all (batch, bit) pairs
    for (size_t j = 0; j < batch; ++j) {
      const ArithLabel &x = xs[j];
      for (size_t i = 0; i < sec_param; ++i) {
        const size_t t = j * sec_param + i;

        // exps[2*t] = (x.auth * rands[i] - x.raw * rands_times_K[i]) mod p'
        // exps[2*t+1] = same mod q'
        fmpz_init(exps[2 * t]);
        fmpz_init(exps[2 * t + 1]);

        fmpz_t tmp1, tmp2, diff;
        fmpz_init(tmp1);
        fmpz_init(tmp2);
        fmpz_init(diff);

        fmpz_mul(tmp1, x.auth, rands[i]);         // x.auth * rands[i]
        fmpz_mul(tmp2, x.raw, rands_times_K[i]);  // x.raw * rands_times_K[i]
        fmpz_sub(diff, tmp1, tmp2);               // diff

        fmpz_mod_set_fmpz(exps[2 * t], diff, &sk.mod_p_prime_ctx);
        fmpz_mod_set_fmpz(exps[2 * t + 1], diff, &sk.mod_q_prime_ctx);

        fmpz_clear(tmp1);
        fmpz_clear(tmp2);
        fmpz_clear(diff);

        fmpz_init(results[2 * t]);
        fmpz_init(results[2 * t + 1]);

        exp_tables[2 * t] = &sk.table_g_p2;
        exp_tables[2 * t + 1] = &sk.table_g_q2;
      }
    }

    // Single batch exponentiation
    powm_precomputed_batch(results, (const fmpz_t *)exps,
                           (const PowmPrecomputeTable **)exp_tables, 2 * total);

    // Reconstruct and extract bits in a single parallel-for
    fmpz_t *reconstructed = new fmpz_t[total];
    for (size_t t = 0; t < total; ++t) {
      fmpz_init(reconstructed[t]);
    }

    uint8_t *y_bits = new uint8_t[total];

#pragma omp parallel for
    for (long long t = 0; t < (long long)total; ++t) {
      const size_t j = (size_t)t / sec_param;
      const size_t i = (size_t)t % sec_param;

      // CRT reconstruct from p- and q- branches
      crt_reconstruct(reconstructed[t], results[2 * t], results[2 * t + 1], sk);

      if (Delta[i]) {
        // *= (1 + N * x.raw)
        fmpz_t one_plus_Nx;
        fmpz_init(one_plus_Nx);
        fmpz_mul(one_plus_Nx, xs[j].raw, pk.N);
        fmpz_add_ui(one_plus_Nx, one_plus_Nx, 1);
        fmpz_mod_mul(reconstructed[t], reconstructed[t], one_plus_Nx,
                     &pk.mod_N2_ctx);
        fmpz_clear(one_plus_Nx);
      }

      // ddlog and take mod 2
      ddlog(reconstructed[t], reconstructed[t], pk);
      y_bits[t] = (uint8_t)fmpz_fdiv_ui(reconstructed[t], 2);
    }

    // Convert bits to Labels per batch element
    for (size_t j = 0; j < batch; ++j) {
      ys[j] = bits_to_label(&y_bits[j * sec_param]);
    }

    // Cleanup
    for (size_t t = 0; t < 2 * total; ++t) {
      fmpz_clear(exps[t]);
      fmpz_clear(results[t]);
    }
    for (size_t t = 0; t < total; ++t) {
      fmpz_clear(reconstructed[t]);
    }

    delete[] exps;
    delete[] results;
    delete[] exp_tables;
    delete[] reconstructed;
    delete[] y_bits;
  }

  Label eval_mod2(const ArithLabel &x, const PaillierPubKey &pk) const {
    //
    fmpz_t exps[sec_param * 2];
    const PowmPrecomputeTable *exp_tables[sec_param * 2];
    fmpz_t results[sec_param * 2];
    for (size_t i = 0; i < sec_param; ++i) {
      fmpz_init_set(exps[2 * i], x.auth);
      fmpz_init_set(exps[2 * i + 1], x.raw);
      fmpz_init(results[2 * i]);
      fmpz_init(results[2 * i + 1]);
      exp_tables[2 * i] = tables[i][0];
      exp_tables[2 * i + 1] = tables[i][1];
    }
    // batch powm_precomputed
    powm_precomputed_batch(results, (const fmpz_t *)exps,
                           (const PowmPrecomputeTable **)exp_tables,
                           sec_param * 2);
    fmpz_t reconstructed[sec_param];
    for (size_t i = 0; i < sec_param; ++i) {
      fmpz_init(reconstructed[i]);
    }
    uint8_t y_bool[sec_param];
#pragma omp parallel for
    for (size_t i = 0; i < sec_param; ++i) {
      fmpz_mod_mul(reconstructed[i], results[2 * i], results[2 * i + 1],
                   &pk.mod_N2_ctx);
      // ddlog
      ddlog(reconstructed[i], reconstructed[i], pk);
      // set y_bool = reconstructed[i] mod 2
      y_bool[i] = fmpz_fdiv_ui(reconstructed[i], 2);
    }
    for (size_t i = 0; i < sec_param * 2; ++i) {
      fmpz_clear(exps[i]);
      fmpz_clear(results[i]);
    }
    for (size_t i = 0; i < sec_param; ++i) {
      fmpz_clear(reconstructed[i]);
    }
    Label y = bits_to_label(y_bool);
    return y;
  }

  // Batch version: single powm_precomputed_batch and single omp parallel for
  void eval_mod2_batch(const ArithLabel *xs, size_t batch, Label *ys,
                       const PaillierPubKey &pk) const {
    const size_t total = batch * sec_param;

    // Allocate exponent and result arrays for the entire batch
    fmpz_t *exps = new fmpz_t[2 * total];
    fmpz_t *results = new fmpz_t[2 * total];
    const PowmPrecomputeTable **exp_tables =
        new const PowmPrecomputeTable *[2 * total];

    // Build exps and tables: reuse tables[i][0/1] across batch
    for (size_t j = 0; j < batch; ++j) {
      const ArithLabel &x = xs[j];
      for (size_t i = 0; i < sec_param; ++i) {
        const size_t t = j * sec_param + i;

        fmpz_init_set(exps[2 * t], x.auth);     // exponent for cts[i].ct[0]
        fmpz_init_set(exps[2 * t + 1], x.raw);  // exponent for cts[i].ct[1]

        fmpz_init(results[2 * t]);
        fmpz_init(results[2 * t + 1]);

        exp_tables[2 * t] = tables[i][0];
        exp_tables[2 * t + 1] = tables[i][1];
      }
    }

    // Single batch exponentiation
    powm_precomputed_batch(results, (const fmpz_t *)exps,
                           (const PowmPrecomputeTable **)exp_tables, 2 * total);

    // Reconstruct and extract bits in a single parallel-for
    fmpz_t *reconstructed = new fmpz_t[total];
    for (size_t t = 0; t < total; ++t) {
      fmpz_init(reconstructed[t]);
    }

    uint8_t *y_bits = new uint8_t[total];

#pragma omp parallel for
    for (long long t = 0; t < (long long)total; ++t) {
      // Multiply the two branches modulo N^2
      fmpz_mod_mul(reconstructed[t], results[2 * t], results[2 * t + 1],
                   &pk.mod_N2_ctx);

      // ddlog and take mod 2
      ddlog(reconstructed[t], reconstructed[t], pk);
      y_bits[t] = (uint8_t)fmpz_fdiv_ui(reconstructed[t], 2);
    }

    // Pack per-item bits into Labels
    for (size_t j = 0; j < batch; ++j) {
      ys[j] = bits_to_label(&y_bits[j * sec_param]);
    }

    // Cleanup
    for (size_t t = 0; t < 2 * total; ++t) {
      fmpz_clear(exps[t]);
      fmpz_clear(results[t]);
    }
    for (size_t t = 0; t < total; ++t) {
      fmpz_clear(reconstructed[t]);
    }
    delete[] exps;
    delete[] results;
    delete[] exp_tables;
    delete[] reconstructed;
    delete[] y_bits;
  }

  ~Mod2Ctx() {
    for (size_t i = 0; i < sec_param; ++i) {
      fmpz_clear(rands[i]);
      fmpz_clear(rands_times_K[i]);
      if (tables[i][0]) {
        tables[i][0]->cleanup();
        delete tables[i][0];
      }
      if (tables[i][1]) {
        tables[i][1]->cleanup();
        delete tables[i][1];
      }
    }
  }
};
}  // namespace ZebraGRAM