#include "paillier.hpp"

#include "util.hpp"

PaillierPrivKey default_sk;
PaillierPubKey default_pk;
fmpz_t default_random_mask_raw;
fmpz_t default_random_mask_auth;
bool default_keypair_initialized = false;

struct DefaultKeypairInitializer {
  DefaultKeypairInitializer() {
    if (!default_keypair_initialized) {
      paillier_keygen(default_pk, default_sk, 256);  // 256 MB precompute tables
      secure_random_fmpz(default_random_mask_raw, LEN_RAW_SHARE);
      secure_random_fmpz(default_random_mask_auth, LEN_AUTH_SHARE);
      default_keypair_initialized = true;
    }
  }
};

DefaultKeypairInitializer default_keypair_initializer;

void sample_prime_prob(fmpz_t p, int bits) {
  do {
    secure_random_fmpz(p, bits);
    // make odd
    if (fmpz_is_even(p)) fmpz_add_ui(p, p, 1);
  } while (!fmpz_is_probabprime(p));
}

// sample prime
void sample_prime(fmpz_t p, int bits) {
  // generate random odd number of specified bits and test for primality
  do {
    sample_prime_prob(p, bits);
  } while (!fmpz_is_prime(p));  // ensure definitely prime
}

void sample_safe_prime_prob(fmpz_t p, int bits) {
  fmpz_t q;
  fmpz_init(q);
  do {
    secure_random_fmpz(q, bits - 1);
    // make odd
    if (fmpz_is_even(q)) fmpz_add_ui(q, q, 1);
    // p = 2q + 1
    fmpz_mul_2exp(p, q, 1);
    fmpz_add_ui(p, p, 1);
  } while (!fmpz_is_probabprime(p) || !fmpz_is_probabprime(q));
  fmpz_clear(q);
}

// Generate a Paillier keypair with modulus size of bits
void paillier_keygen(PaillierPubKey &pk, PaillierPrivKey &sk,
                     size_t precompute_mb) {
  fmpz_t p, q;
  fmpz_init(p);
  fmpz_init(q);
#if defined(TESTING)
#if (LEN_PRIME == 1560)
  // fmpz_set_str(sk.p,
  // "21019291990590200087044852906748028682122556773934469331050211920737553651275991370523419465549507203861227823404904409194594374021033509251972644106258933455543160367342109253097739437571297040126572251587566354248542778485174322898201017886929491962198403535885288422627892254859282090215307200661096504100762663759729019190039619459750020091517313046319931235065778878189205239367879815591474332713118908161515450313430255691698087351489868188817609528450562227325563",
  // 10); fmpz_set_str(sk.q,
  // "25803698197120525590032358582620562168217222651133212206134966726942204474015782687727008855129105742019356860374078509392406844300748372466746152736153325407650494843617450865199229167080128619164800512727341788337003788886790677268864935809641238471426409358569131284987635718193106038764559968997647606826747219900498314422749841641044863369363721824698666155557633613469825149615158360494372847829403632341731960005847600779460707802173537498460614794008182039420383",
  // 10);
  static const ulong p_limbs[] = {0x54c9ec78586e4a7b,
                                  0x87b305b74cdecf89,
                                  0xe004ad853dcae8b,
                                  0xafb3eab70ebb9ea5,
                                  0xc143d92ed85c0019,
                                  0x8ef3a51d12b4c2bd,
                                  0x5f51f37568cb13c1,
                                  0x253d202068148be6,
                                  0xb6f62eb8296c4ecb,
                                  0x263c5e8d96bae757,
                                  0x4d15d07e78d79217,
                                  0x45942b6410cc18ee,
                                  0x997df10e63db7e27,
                                  0x40004745f7a87bed,
                                  0x56c267c373da0054,
                                  0x7e0d03914c2b1c6f,
                                  0xea5478a9b944cab7,
                                  0x78979e73453ec984,
                                  0x3d1c20146edfea61,
                                  0x8e23d2d7f824300f,
                                  0x4a3d621f5ef1ccf8,
                                  0x56c862438f4c861f,
                                  0x44b3b4f41d13d545,
                                  0x6828dfecac4bc6e1,
                                  0x8510b7};  // 25 limbs
  static const ulong q_limbs[] = {0x82817a82b28e9ddf,
                                  0xeef2ba6a86af75f3,
                                  0x381b07f97a9f77ba,
                                  0x1ec15f08868e1ca,
                                  0x198b513a1cb18e4e,
                                  0xfc1f423d03b3af97,
                                  0x8fec9c6d9def2be8,
                                  0xa5db5ceb7c2c6b11,
                                  0x8d8bef49669febc5,
                                  0xaac5c5f926418acd,
                                  0xb8ea28a655d54a4,
                                  0x6ce5e955e01e54e1,
                                  0xdbc35aebbb2f6a4d,
                                  0x6f15a945199c8860,
                                  0xcbed5df33c1ce6e2,
                                  0x56c4ab6b1b5f256d,
                                  0xd32ce58825527eae,
                                  0x5fdb1e78e430b883,
                                  0x97b1863a5c33193f,
                                  0xad5da2242dd8febf,
                                  0x4699f1ef26bf430e,
                                  0x4afe1d0570ff0831,
                                  0x8dbadf151e586949,
                                  0xdeda4462ecbd5e6b,
                                  0xa35a84};  // 25 limbs
  fmpz_set_ui_array(p, p_limbs, 25);
  fmpz_set_ui_array(q, q_limbs, 25);
#elif (LEN_PRIME == 1024)
  // Hardcoded primes as ulong arrays:
  static const ulong p_limbs[] = {
      0xc5e5135fc5fae9db, 0x99f3eda40616208a, 0x4d727c79760a359b,
      0x85d11bd180ed2980, 0xa46b73cd31fd03b,  0xc91697dd4eb7ee2e,
      0x2cf172d3b230179b, 0xc5d8780be84cf62b, 0x8050b257546b9892,
      0x5b65fecf7017b910, 0x8412036981d29a22, 0x7334e5867772275f,
      0xfb09969d4be5f14e, 0x18a8fdc5f53f8e24, 0x811121951da60b51,
      0xe6678d7d6293cfaa};  // 16 limbs
  static const ulong q_limbs[] = {
      0x6577ded9c92e4f3f, 0xa3477fee60b24aac, 0xff3c932c88e96973,
      0xe9f8c5d6799467c7, 0x97762b21f3329326, 0x4c47959c39ec4697,
      0xfac44085d33ddc91, 0x81932c2f9a068ce1, 0x2375c633b8701040,
      0x429262da28c6f52b, 0x1e3c563a8655563c, 0x7651a05ef5af629,
      0xdc895ea4ade06aa2, 0x8a70593d6fb449cc, 0x48fa26182d6a31c6,
      0xec853123b7481847};  // 16 limbs
  fmpz_set_ui_array(p, p_limbs, 16);
  fmpz_set_ui_array(q, q_limbs, 16);
#endif
#else
  // generate two distinct primes p,q each of size bits/2
  do {
    sample_safe_prime_prob(p, LEN_PRIME);
    sample_safe_prime_prob(q, LEN_PRIME);
  } while (fmpz_cmp(p, q) == 0);
#ifdef PRINT_HARDCODED_PRIMES
  printf("// Hardcoded primes as ulong arrays:\n");

  // Print p as ulong array
  size_t p_limbs_size = fmpz_size(sk.p);
  ulong *p_limb_array = (ulong *)malloc(p_limbs_size * sizeof(ulong));
  fmpz_get_ui_array(p_limb_array, p_limbs_size, sk.p);
  printf("static const ulong p_limbs[] = {");
  for (size_t i = 0; i < p_limbs_size; i++) {
    if (i > 0) printf(", ");
    printf("0x%lx", p_limb_array[i]);
  }
  printf("}; // %zu limbs\n", p_limbs_size);
  free(p_limb_array);

  // Print q as ulong array
  size_t q_limbs_size = fmpz_size(sk.q);
  ulong *q_limb_array = (ulong *)malloc(q_limbs_size * sizeof(ulong));
  fmpz_get_ui_array(q_limb_array, q_limbs_size, sk.q);
  printf("static const ulong q_limbs[] = {");
  for (size_t i = 0; i < q_limbs_size; i++) {
    if (i > 0) printf(", ");
    printf("0x%lx", q_limb_array[i]);
  }
  printf("}; // %zu limbs\n", q_limbs_size);
  free(q_limb_array);
#endif
#endif
  // print p and q as an unsigned long array here for testing
  // N = p * q
  fmpz_mul(pk.N, p, q);
  // N2 = N^2
  fmpz_mul(pk.N2, pk.N, pk.N);

  // compute p_prime = (p-1)/2, q_prime = (q-1)/2
  fmpz_t p_prime, q_prime;
  fmpz_init(p_prime);
  fmpz_init(q_prime);
  fmpz_sub_ui(p_prime, p, 1);
  fmpz_fdiv_q_2exp(p_prime, p_prime, 1);
  fmpz_sub_ui(q_prime, q, 1);
  fmpz_fdiv_q_2exp(q_prime, q_prime, 1);

  // initialize mod N^2 context
  fmpz_mod_ctx_init(&pk.mod_N2_ctx, pk.N2);
  fmpz_mod_ctx_init(&pk.mod_N_ctx, pk.N);
  // initialize mod p' and q' context
  fmpz_mod_ctx_init(&sk.mod_p_prime_ctx, p_prime);
  fmpz_mod_ctx_init(&sk.mod_q_prime_ctx, q_prime);
  fmpz_clear(p_prime);
  fmpz_clear(q_prime);

  secure_random_fmpz(pk.g, LEN_N);
  // g = r^{2N}
  fmpz_mul(pk.g, pk.g, pk.g);

  // could be optimized with CRT, but only done once in keygen, so not a big
  // deal
  powm(pk.g, pk.g, pk.N, pk.N2);

  fmpz_set(sk.g, pk.g);

  // random K of LEN_K bits
  secure_random_fmpz(sk.K, LEN_K);
  // neg k
  // compute pk = g^\K mod N^2
  powm(pk.pk, pk.g, sk.K, pk.N2);

  // p2 = p^2
  fmpz_mul(sk.p2, p, p);
  // q2 = q^2
  fmpz_mul(sk.q2, q, q);
  fmpz_clear(p);
  fmpz_clear(q);

  // initialize modular contexts
  fmpz_mod_ctx_init(&sk.mod_p2_ctx, sk.p2);
  fmpz_mod_ctx_init(&sk.mod_q2_ctx, sk.q2);
  // p2inv = p^{-1} mod q
  if (fmpz_invmod(sk.p2inv, sk.p2, sk.q2) == 0) {
    fprintf(stderr, "Error: fmpz_invert failed in paillier_keygen\n");
    exit(1);
  }

  fmpz_t inv_pk_p2;
  fmpz_t inv_pk_q2;
  fmpz_init(inv_pk_p2);
  fmpz_init(inv_pk_q2);
  if (fmpz_invmod(inv_pk_p2, pk.pk, sk.p2) == 0) {
    fprintf(stderr, "Error: fmpz_invert failed in paillier_keygen\n");
    exit(1);
  }
  if (fmpz_invmod(inv_pk_q2, pk.pk, sk.q2) == 0) {
    fprintf(stderr, "Error: fmpz_invert failed in paillier_keygen\n");
    exit(1);
  }

  // precompute table for g^{x} mod p^2 and g^{x} mod q^2
  // also compute a small table of (g^K)-th powers for fast computing the
  // ciphertext
  const size_t aux_table_size = precompute_mb * LEN_K / (LEN_K + LEN_PRIME);
  const size_t main_table_size = precompute_mb - aux_table_size;
  // parallelize using omp
#pragma omp parallel sections
  {
#pragma omp section
    {
      sk.table_g_p2.init(pk.g, sk.p2, (size_t)LEN_PRIME, main_table_size / 2);
      flint_cleanup();
    }
#pragma omp section
    {
      sk.table_g_q2.init(pk.g, sk.q2, (size_t)LEN_PRIME, main_table_size / 2);
      flint_cleanup();
    }
#pragma omp section
    {
      sk.table_pk_p2.init(inv_pk_p2, sk.p2, (size_t)LEN_K, aux_table_size / 2);
      flint_cleanup();
    }
#pragma omp section
    {
      sk.table_pk_q2.init(inv_pk_q2, sk.q2, (size_t)LEN_K, aux_table_size / 2);
      flint_cleanup();
    }
  }
  fmpz_clear(inv_pk_p2);
  fmpz_clear(inv_pk_q2);
}

// CRT reconstruction: combines values mod p^2 and mod q^2 into value mod N^2
void crt_reconstruct(fmpz_t result, const fmpz_t val_p2, const fmpz_t val_q2,
                     const PaillierPrivKey &sk) {
  // result = val_p2 + p2 * ((val_q2 - val_p2) * p2inv mod q2)
  fmpz_sub(result, val_q2, val_p2);
  fmpz_mod_mul(result, result, sk.p2inv, &sk.mod_q2_ctx);
  fmpz_mul(result, result, sk.p2);
  fmpz_add(result, val_p2, result);
}

// Fast inverse mod p^2 using Hensel lifting
void fast_inv_mod_prime2(fmpz_t result, const fmpz_t a,
                         const fmpz_mod_ctx_struct *mod_p_ctx,
                         const fmpz_mod_ctx_struct *mod_p2_ctx) {
  const fmpz *p = fmpz_mod_ctx_modulus(mod_p_ctx);

  // First compute a^(-1) mod p using Fermat's little theorem: a^(p-2) mod p
  fmpz_t inv_p, p_minus_2;
  fmpz_init(inv_p);
  fmpz_init(p_minus_2);

  fmpz_sub_ui(p_minus_2, p, 2);
  fmpz_mod_inv(inv_p, a, mod_p_ctx);  // a^(p-2) mod p = a^(-1) mod p

  // Hensel lift: if x ≡ a^(-1) (mod p), then a^(-1) mod p^2 = x(2 - ax) mod p^2
  fmpz_t ax, two_minus_ax;
  fmpz_init(ax);
  fmpz_init(two_minus_ax);

  fmpz_mod_mul(ax, a, inv_p, mod_p2_ctx);  // ax mod p^2
  fmpz_set_ui(two_minus_ax, 2);
  fmpz_mod_sub(two_minus_ax, two_minus_ax, ax, mod_p2_ctx);  // (2 - ax) mod p^2
  fmpz_mod_mul(result, inv_p, two_minus_ax, mod_p2_ctx);  // x(2 - ax) mod p^2

  fmpz_clear(inv_p);
  fmpz_clear(p_minus_2);
  fmpz_clear(ax);
  fmpz_clear(two_minus_ax);
}

void paillier_encrypt_with_sk(PaillierCiphertext &ct, const fmpz_t m,
                              const PaillierPubKey &pk,
                              const PaillierPrivKey &sk, const fmpz_t r) {
  // sample a random LEN_K bit r
  fmpz_t ct_p2[2], ct_q2[2];
  for (int i = 0; i < 2; i++) {
    fmpz_init(ct_p2[i]);
    fmpz_init(ct_q2[i]);
  }
// batch powerings with precomputed table
#pragma omp parallel sections
  {
#pragma omp section
    {
      powm_precomputed(ct_p2[0], r, sk.table_g_p2);
      // powm_batch(ct_p2, pk.g, exps_p2, sk.p2, 2);
    }
#pragma omp section
    {
      powm_precomputed(ct_q2[0], r, sk.table_g_q2);
      // powm_batch(ct_q2, pk.g, exps_q2, sk.q2, 2);
    }
#pragma omp section
    {
      powm_precomputed(ct_p2[1], r, sk.table_pk_p2);
      // powm_batch(ct_p2, pk.pk, exps_p2, sk.p2, 2);
      // let ct_p2[1] = inv(ct_p2[1]) mod p2 using fast Hensel lifting
      // fast_inv_mod_prime2(ct_p2[1], ct_p2[1], &sk.mod_p_ctx, &sk.mod_p2_ctx);
    }
#pragma omp section
    {
      powm_precomputed(ct_q2[1], r, sk.table_pk_q2);
      // powm_batch(ct_q2, pk.pk, exps_q2, sk.q2, 2);
      // let ct_q2[1] = inv(ct_q2[1]) mod q2 using fast Hensel lifting
      // fast_inv_mod_prime2(ct_q2[1], ct_q2[1], &sk.mod_q_ctx, &sk.mod_q2_ctx);
    }
  }
  // combine with CRT: ct.ct[i] = ct_p2[i] + p2 * ((ct_q2[i] - ct_p2[i]) * p2inv
  // mod q2)
  for (int i = 0; i < 2; i++) {
    crt_reconstruct(ct.ct[i], ct_p2[i], ct_q2[i], sk);
  }

  // let ct.ct[1] = ct.ct[1] * (1 + m * N) mod N^2
  fmpz_t one_plus_mN;
  fmpz_init(one_plus_mN);
  fmpz_mul(one_plus_mN, m, pk.N);
  fmpz_add_ui(one_plus_mN, one_plus_mN, 1);
  fmpz_mod_mul(ct.ct[1], ct.ct[1], one_plus_mN, &pk.mod_N2_ctx);
  fmpz_clear(one_plus_mN);
  for (int i = 0; i < 2; i++) {
    fmpz_clear(ct_p2[i]);
    fmpz_clear(ct_q2[i]);
  }
}

// void PaillierEncryptCtx::iter_1(size_t i, const PaillierPrivKey &sk, const
// fmpz_t *r)
// {
//     size_t group = i / (total_iter_1 / 4);
//     size_t index = i % (total_iter_1 / 4);
//     switch (group)
//     {
//     case 0:
//         sk.table_g_p2.powm(ct_p2[0][i], r[index]);
//         break;
//     case 1:
//         sk.table_g_q2.powm(ct_q2[0][index], r[index]);
//         break;
//     case 2:
//         sk.table_pk_p2.powm(ct_p2[1][index], r[index]);
//         fast_inv_mod_prime2(ct_p2[1][index], ct_p2[1][index], &sk.mod_p_ctx,
//         &sk.mod_p2_ctx); break;
//     case 3:
//         sk.table_pk_q2.powm(ct_q2[1][index], r[index]);
//         fast_inv_mod_prime2(ct_q2[1][index], ct_q2[1][index], &sk.mod_q_ctx,
//         &sk.mod_q2_ctx); break;
//     default:
//         break;
//     }
// }

void PaillierEncryptCtx::push_tasks(
    const fmpz_t *r, const PaillierPrivKey &sk,
    std::vector<const fmpz_t *> &all_exps,
    std::vector<const PowmPrecomputeTable *> &all_tables) {
  for (size_t i = 0; i < batch_size; i++) {
    all_exps.push_back(&r[i]);
    all_exps.push_back(&r[i]);
    all_exps.push_back(&r[i]);
    all_exps.push_back(&r[i]);
    all_tables.push_back(&sk.table_g_p2);
    all_tables.push_back(&sk.table_g_q2);
    all_tables.push_back(&sk.table_pk_p2);
    all_tables.push_back(&sk.table_pk_q2);
  }
}

void PaillierEncryptCtx::retrieve_results(
    size_t i, PaillierCiphertext *ct, const fmpz_t *m, const PaillierPubKey &pk,
    const PaillierPrivKey &sk, const fmpz_t &result_p2_0,
    const fmpz_t &result_q2_0, const fmpz_t &result_p2_1,
    const fmpz_t &result_q2_1) {
  fmpz_set(ct_p2[0][i], result_p2_0);
  fmpz_set(ct_q2[0][i], result_q2_0);
  fmpz_set(ct_p2[1][i], result_p2_1);
  fmpz_set(ct_q2[1][i], result_q2_1);
  // fast_inv_mod_prime2(ct_p2[1][i], ct_p2[1][i], &default_sk.mod_p_ctx,
  // &default_sk.mod_p2_ctx); fast_inv_mod_prime2(ct_q2[1][i], ct_q2[1][i],
  // &default_sk.mod_q_ctx, &default_sk.mod_q2_ctx);
  for (int j = 0; j < 2; j++) {
    crt_reconstruct(ct[i].ct[j], ct_p2[j][i], ct_q2[j][i], sk);
  }
  fmpz_mul(one_plus_mN[i], m[i], pk.N);
  fmpz_add_ui(one_plus_mN[i], one_plus_mN[i], 1);
  fmpz_mod_mul(ct[i].ct[1], ct[i].ct[1], one_plus_mN[i], &pk.mod_N2_ctx);
}

// void PaillierEncryptCtx::iter_2(size_t index, PaillierCiphertext *ct, const
// fmpz_t *m, const PaillierPubKey &pk, const PaillierPrivKey &sk) {
//     for (int j = 0; j < 2; j++)
//     {
//         crt_reconstruct(ct[index].ct[j], ct_p2[j][index], ct_q2[j][index],
//         sk);
//     }
//     fmpz_mul(one_plus_mN[index], m[index], pk.N);
//     fmpz_add_ui(one_plus_mN[index], one_plus_mN[index], 1);
//     fmpz_mod_mul(ct[index].ct[1], ct[index].ct[1], one_plus_mN[index],
//     &pk.mod_N2_ctx);
// }

void paillier_encrypt_with_sk_batch(PaillierCiphertext *ct, const fmpz_t *m,
                                    const PaillierPubKey &pk,
                                    const PaillierPrivKey &sk, const fmpz_t *r,
                                    size_t batch_size) {
  PaillierEncryptCtx context;
  context.init(batch_size);
  std::vector<const fmpz_t *> all_exps;
  std::vector<const PowmPrecomputeTable *> all_tables;
  context.push_tasks(r, sk, all_exps, all_tables);
  fmpz_t *all_results = new fmpz_t[all_exps.size()];
  for (size_t i = 0; i < all_exps.size(); i++) {
    fmpz_init(all_results[i]);
  }
  powm_precomputed_batch(all_results, all_exps.data(), all_tables.data(),
                         all_exps.size());
#pragma omp parallel for
  for (size_t i = 0; i < batch_size; i++) {
    context.retrieve_results(i, ct, m, pk, sk, all_results[4 * i],
                             all_results[4 * i + 1], all_results[4 * i + 2],
                             all_results[4 * i + 3]);
  }
  for (size_t i = 0; i < all_exps.size(); i++) {
    fmpz_clear(all_results[i]);
  }
  delete[] all_results;
  // #pragma omp parallel for
  //     for (size_t i = 0; i < context.total_iter_1; i++)
  //     {
  //         context.iter_1(i, sk, r);
  //     }

  //     for (size_t i = 0; i < context.total_iter_2; i++)
  //     {
  //         context.iter_2(i, ct, m, pk, sk);
  //     }
}

void paillier_encrypt_with_sk(PaillierCiphertext &ct, const fmpz_t m,
                              const PaillierPubKey &pk,
                              const PaillierPrivKey &sk) {
  // sample a random LEN_K bit r
  fmpz_t r;
  fmpz_init(r);
  secure_random_fmpz(r, LEN_K);

  paillier_encrypt_with_sk(ct, m, pk, sk, r);

  fmpz_clear(r);
}

void paillier_generator_pow(fmpz_t result, const fmpz_t exp,
                            const PaillierPrivKey &sk) {
  // compute g^{exp} mod N^2 using CRT
  fmpz_t g_exp_p2, g_exp_q2;
  fmpz_init(g_exp_p2);
  fmpz_init(g_exp_q2);
  fmpz_t exp_p_prime, exp_q_prime;
  fmpz_init(exp_p_prime);
  fmpz_init(exp_q_prime);

  fmpz_mod_set_fmpz(exp_p_prime, exp, &sk.mod_p_prime_ctx);
  fmpz_mod_set_fmpz(exp_q_prime, exp, &sk.mod_q_prime_ctx);
  // #pragma omp parallel sections num_threads(2)
  {
    // #pragma omp section
    {
      powm_precomputed(g_exp_p2, exp_p_prime, sk.table_g_p2);
      // powm(g_exp_p2, sk.g, exp_abs, sk.p2);
    }
    // #pragma omp section
    {
      powm_precomputed(g_exp_q2, exp_q_prime, sk.table_g_q2);
      // powm(g_exp_q2, sk.g, exp_abs, sk.q2);
    }
  }

  // combine using CRT
  crt_reconstruct(result, g_exp_p2, g_exp_q2, sk);
  fmpz_clear(g_exp_p2);
  fmpz_clear(g_exp_q2);
  fmpz_clear(exp_p_prime);
  fmpz_clear(exp_q_prime);
}

void paillier_decrypt(fmpz_t m, const PaillierCiphertext &ct,
                      const PaillierPubKey &pk, const PaillierPrivKey &sk) {
  // Compute ct[0]^K mod N^2 using CRT
  fmpz_t ct0_K_p2, ct0_K_q2;
  fmpz_init(ct0_K_p2);
  fmpz_init(ct0_K_q2);

  // Compute ct[0]^K mod p^2 and mod q^2 in parallel
#pragma omp parallel sections
  {
#pragma omp section
    powm(ct0_K_p2, ct.ct[0], sk.K, pk.N2);

#pragma omp section
    powm(ct0_K_q2, ct.ct[0], sk.K, pk.N2);
  }

  // Combine using CRT
  fmpz_t ct0_K;
  fmpz_init(ct0_K);
  crt_reconstruct(ct0_K, ct0_K_p2, ct0_K_q2, sk);

  // Compute ct0_K * ct[1] mod N^2
  fmpz_t numerator;
  fmpz_init(numerator);
  fmpz_mod_mul(numerator, ct0_K, ct.ct[1], &pk.mod_N2_ctx);

  // Subtract 1: numerator = numerator - 1
  fmpz_sub_ui(numerator, numerator, 1);

  // Divide by N: m = (numerator) / N
  fmpz_divexact(m, numerator, pk.N);

  // Clean up
  fmpz_clear(ct0_K_p2);
  fmpz_clear(ct0_K_q2);
  fmpz_clear(ct0_K);
  fmpz_clear(numerator);
}

void ddlog(fmpz_t result, const fmpz_t x, const PaillierPubKey &pk) {
  // result = (x / (x mod N) - 1) / N

  // Compute x mod N
  fmpz_t tmp;
  fmpz_init(tmp);

  fmpz_mod_set_fmpz(tmp, x, &pk.mod_N_ctx);

  // Compute x / (x mod N) mod N^2
  fmpz_mod_divides(tmp, x, tmp, &pk.mod_N2_ctx);

  // Divide by N: result = (quotient - 1) / N
  fmpz_divexact(result, tmp, pk.N);

  // Clean up
  fmpz_clear(tmp);
}
