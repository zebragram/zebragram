#include "arith_mul.hpp"

namespace ZebraGRAM {

void garble_mul(ArithLabel &z, uint8_t material[GARBLE_MULT_SIZE],
                const ArithLabel &x, const ArithLabel &y,
                const PaillierPubKey &pk, const PaillierPrivKey &sk) {
  // sample r1 and r2
  fmpz_t r[2], m[2];
  fmpz_t g_power, tmp, z_combine;
  fmpz_init(g_power);
  fmpz_init(tmp);
  fmpz_init(z_combine);
  for (int i = 0; i < 2; i++) {
    fmpz_init(r[i]);
    fmpz_init(m[i]);
    secure_random_fmpz(r[i], LEN_K);
  }
  PaillierCiphertext ct[2];
  fmpz_mul(tmp, x.raw, sk.K);
  // shift m[0] left by LEN_RAW_SHARE
  fmpz_mul_2exp(m[0], tmp, LEN_RAW_SHARE);
  // m[0] += x.raw
  fmpz_add(m[0], m[0], x.raw);
  y.merge(m[1]);

  paillier_encrypt_with_sk_batch(ct, m, pk, sk, r, 2);
  // serialize
  ct[0].serialize(material);
  ct[1].serialize(material + PAILLIER_CIPHER_TEXT_SIZE);

  // g's power: r[0] * (y_auth - K * y_raw) + r[1] * (x_auth - K * x_raw)
  fmpz_sub(tmp, x.auth, tmp);    // tmp = x_auth - K * x_raw
  fmpz_mul(g_power, r[1], tmp);  // g_power = r[1] * (x_auth - K * x_raw)
  fmpz_mul(tmp, sk.K, y.raw);    // tmp = K * y_raw
  fmpz_sub(tmp, y.auth, tmp);    // tmp = y_auth - K * y_raw
  fmpz_mul(tmp, r[0], tmp);      // tmp = r[0] * (y_auth - K * y_raw)

  fmpz_add(g_power, g_power, tmp);  // g_power += tmp

  paillier_generator_pow(z_combine, g_power, sk);

  ddlog(z_combine, z_combine, pk);
  // z_combine = - (z_combine + m[0] * Ry_raw) mod N where m[0] = Rx.raw * Delta
  fmpz_mul(m[0], m[0], y.raw);
  fmpz_add(z_combine, z_combine, m[0]);
  fmpz_neg(z_combine, z_combine);
  fmpz_mod_set_fmpz(z_combine, z_combine, &pk.mod_N_ctx);

  z.split(z_combine);

  // clear
  for (int i = 0; i < 2; ++i) {
    fmpz_clear(r[i]);
    fmpz_clear(m[i]);
  }
  fmpz_clear(g_power);
  fmpz_clear(tmp);
  fmpz_clear(z_combine);
}

// z_vec = x * y_vec
size_t garble_mul_vec(ArithLabel *z_vec, uint8_t *material, const ArithLabel &x,
                      const ArithLabel *y_vec, size_t vec_size,
                      const PaillierPubKey &pk, const PaillierPrivKey &sk) {
  // Allocate arrays: m[0] is shared, m[1..vec_size] are for each y_vec element
  GarbleMulVecCtx context;
  context.init(x, y_vec, vec_size, sk);
  std::vector<const fmpz_t *> all_exps;
  std::vector<const PowmPrecomputeTable *> all_tables;
  size_t inner_iteration_1 = context.total_iter_1;
  size_t inner_iteration_2 = context.total_iter_2;
  context.push_iter_1_tasks(sk, all_exps, all_tables);

  size_t iter_2_begin_index = all_exps.size();
  all_exps.resize(inner_iteration_2 * 2 + iter_2_begin_index);
  all_tables.resize(inner_iteration_2 * 2 + iter_2_begin_index);

  for (size_t j = 0; j < inner_iteration_2; ++j) {
    size_t idx = iter_2_begin_index + j * 2;
    context.set_iter_2_tasks(j, y_vec, default_sk, all_exps[idx],
                             all_exps[idx + 1], all_tables[idx],
                             all_tables[idx + 1]);
  }

  size_t total_tasks = all_exps.size();
  fmpz_t *all_powm_results = new fmpz_t[total_tasks];
  for (size_t i = 0; i < total_tasks; i++) {
    fmpz_init(all_powm_results[i]);
  }
  powm_precomputed_batch(all_powm_results, all_exps.data(), all_tables.data(),
                         total_tasks);
  size_t iter_1_total_ops = context.total_iter_1 * 4;
  // Retrieve results for iter_1

#pragma omp parallel for
  for (size_t j = 0; j < inner_iteration_1 + inner_iteration_2; j++) {
    if (j < inner_iteration_1) {
      size_t res_start_index = j * 4;
      context.retrieve_iter_1_results(j, material, pk, sk,
                                      all_powm_results[res_start_index],
                                      all_powm_results[res_start_index + 1],
                                      all_powm_results[res_start_index + 2],
                                      all_powm_results[res_start_index + 3]);
    } else {
      size_t res_start_index = iter_1_total_ops + (j - inner_iteration_1) * 2;
      context.retrieve_iter_2_results(j - inner_iteration_1, z_vec, y_vec, pk,
                                      sk, all_powm_results[res_start_index],
                                      all_powm_results[res_start_index + 1]);
    }
  }
  // clear
  for (size_t i = 0; i < total_tasks; i++) {
    fmpz_clear(all_powm_results[i]);
  }
  delete[] all_powm_results;

  return (vec_size + 1) * PAILLIER_CIPHER_TEXT_SIZE;
}

void eval_mul(ArithLabel &z, const uint8_t material[GARBLE_MULT_SIZE],
              const ArithLabel &x, const ArithLabel &y,
              const PaillierPubKey &pk) {
  PaillierCiphertext ct[2];
  ct[0].deserialize(material);
  ct[1].deserialize(material + PAILLIER_CIPHER_TEXT_SIZE);

  const fmpz_t *exps[4];
  exps[0] = &y.auth;
  exps[1] = &y.raw;
  exps[2] = &x.auth;
  exps[3] = &x.raw;

  const fmpz_t *bases[4];
  bases[0] = &ct[0].ct[0];
  bases[1] = &ct[0].ct[1];
  bases[2] = &ct[1].ct[0];
  bases[3] = &ct[1].ct[1];

  fmpz_t powm_result[4];
  for (int i = 0; i < 4; i++) {
    fmpz_init(powm_result[i]);
  }
  // batch powerings
  powm_batch(powm_result, bases, exps, pk.N2, 4);

  // combine results to powm_result[0]
  for (int i = 1; i < 4; i++) {
    fmpz_mod_mul(powm_result[0], powm_result[0], powm_result[i],
                 &pk.mod_N2_ctx);
  }
  // store ddlog result in powm_result[1]
  ddlog(powm_result[1], powm_result[0], pk);
  // combine y to powm_result[2]
  y.merge(powm_result[2]);
  // multiply x.raw * y_combine
  fmpz_mul(powm_result[2], x.raw, powm_result[2]);
  // powm_result[3] = powm_result[2] - powm_result[1]
  fmpz_sub(powm_result[3], powm_result[2], powm_result[1]);
  fmpz_mod_set_fmpz(powm_result[3], powm_result[3], &pk.mod_N_ctx);
  // split to z
  z.split(powm_result[3]);
}

// Optimized eval_mul_vec that collects all operations into a single large batch
void eval_mul_vec_batched(ArithLabel *z_vec, const uint8_t *material,
                          const ArithLabel &x, const ArithLabel *y_vec,
                          size_t vec_size, const PaillierPubKey &pk) {
  // Calculate total operations: 4 * vec_size (4 operations per vector element)
  size_t total_ops = 4 * vec_size;

  // Allocate arrays for the large batch operation
  const fmpz_t **all_bases = new const fmpz_t *[total_ops];
  const fmpz_t **all_exps = new const fmpz_t *[total_ops];
  fmpz_t *all_powm_results = new fmpz_t[total_ops];

  // Temporary storage for deserialized ciphertexts
  fmpz_t *ct_bases_1 = new fmpz_t[vec_size];
  fmpz_t *ct_bases_2 = new fmpz_t[vec_size];
  fmpz_t base_3, base_4;
  fmpz_init(base_3);
  fmpz_init(base_4);

  // Initialize result arrays
  for (size_t i = 0; i < total_ops; i++) {
    fmpz_init(all_powm_results[i]);
  }

  // Initialize and deserialize ciphertext bases
  for (size_t i = 0; i < vec_size; i++) {
    fmpz_init(ct_bases_1[i]);
    fmpz_init(ct_bases_2[i]);
    PaillierCiphertext::deserialize(
        ct_bases_1[i], ct_bases_2[i],
        material + (i + 1) * PAILLIER_CIPHER_TEXT_SIZE);
  }
  PaillierCiphertext::deserialize(base_3, base_4, material);

  // // print out all the bases
  // std::cout << "Deserialized bases:" << std::endl;
  // for (size_t i = 0; i < vec_size; i++) {
  //     std::cout << "ct_bases_1[" << i << "]: ";
  //     fmpz_print(ct_bases_1[i]);
  //     std::cout << std::endl;
  //     std::cout << "ct_bases_2[" << i << "]: ";
  //     fmpz_print(ct_bases_2[i]);
  //     std::cout << std::endl;
  // }
  // std::cout << "base_3: ";
  // fmpz_print(base_3);
  // std::cout << std::endl;
  // std::cout << "base_4: ";
  // fmpz_print(base_4);
  // std::cout << std::endl;
  // Populate the batch arrays
  for (size_t i = 0; i < vec_size; i++) {
    size_t base_idx = i * 4;

    // Operation 1: ct_bases_1[i] ^ x.auth
    all_bases[base_idx] = &ct_bases_1[i];
    all_exps[base_idx] = &x.auth;

    // Operation 2: ct_bases_2[i] ^ x.raw
    all_bases[base_idx + 1] = &ct_bases_2[i];
    all_exps[base_idx + 1] = &x.raw;

    // Operation 3: base_3 ^ y_vec[i].auth
    all_bases[base_idx + 2] = &base_3;
    all_exps[base_idx + 2] = &y_vec[i].auth;

    // Operation 4: base_4 ^ y_vec[i].raw
    all_bases[base_idx + 3] = &base_4;
    all_exps[base_idx + 3] = &y_vec[i].raw;
  }

  // Single large batch operation
  powm_batch(all_powm_results, all_bases, all_exps, pk.N2, total_ops);

  // Process results for each vector element
  for (size_t i = 0; i < vec_size; i++) {
    size_t base_idx = i * 4;

    // Combine the 4 results for this element
    fmpz_t combined_result, ddlog_result, y_combined, final_result;
    fmpz_init(combined_result);
    fmpz_init(ddlog_result);
    fmpz_init(y_combined);
    fmpz_init(final_result);

    // combined_result = result1 * result2 * result3 * result4
    fmpz_mod_mul(combined_result, all_powm_results[base_idx],
                 all_powm_results[base_idx + 1], &pk.mod_N2_ctx);
    fmpz_mod_mul(combined_result, combined_result,
                 all_powm_results[base_idx + 2], &pk.mod_N2_ctx);
    fmpz_mod_mul(combined_result, combined_result,
                 all_powm_results[base_idx + 3], &pk.mod_N2_ctx);

    // Compute discrete log
    ddlog(ddlog_result, combined_result, pk);

    // Combine y_vec[i] and multiply by x.raw
    y_vec[i].merge(y_combined);
    fmpz_mul(y_combined, x.raw, y_combined);

    // final_result = y_combined - ddlog_result
    fmpz_sub(final_result, y_combined, ddlog_result);
    fmpz_mod_set_fmpz(final_result, final_result, &pk.mod_N_ctx);

    // Split result into z_vec[i]
    z_vec[i].split(final_result);

    // Clean up temporary variables
    fmpz_clear(combined_result);
    fmpz_clear(ddlog_result);
    fmpz_clear(y_combined);
    fmpz_clear(final_result);
  }

  // Clean up arrays
  for (size_t i = 0; i < total_ops; i++) {
    fmpz_clear(all_powm_results[i]);
  }
  for (size_t i = 0; i < vec_size; i++) {
    fmpz_clear(ct_bases_1[i]);
    fmpz_clear(ct_bases_2[i]);
  }
  fmpz_clear(base_3);
  fmpz_clear(base_4);

  delete[] all_bases;
  delete[] all_exps;
  delete[] all_powm_results;
  delete[] ct_bases_1;
  delete[] ct_bases_2;
}

void eval_mul_vec(ArithLabel *z_vec, const uint8_t *material,
                  const ArithLabel &x, const ArithLabel *y_vec, size_t vec_size,
                  const PaillierPubKey &pk) {
  // Use the optimized batched version instead of the old implementation
  eval_mul_vec_batched(z_vec, material, x, y_vec, vec_size, pk);
}

}  // namespace ZebraGRAM