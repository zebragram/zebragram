#include <gtest/gtest.h>
#include <omp.h>

#include "arith_label.hpp"
#include "paillier.hpp"

using namespace ZebraGRAM;

TEST(TestPaillier, Memory) {
  // Create a vector of arith labels, and see if they get properly destructed
  // with clear_and_release
  const int N = 10000;
  std::cout << "Test Begins" << std::endl;
  _fmpz_cleanup();
  raise_sigusr1();
  // first create an array of fmpz_t and initialize them
  fmpz_t *arr = new fmpz_t[N];
  for (int i = 0; i < N; ++i) {
    fmpz_init(arr[i]);
  }
  std::cout << "After Init" << std::endl;
  raise_sigusr1();
  // now assign some random values
  for (int i = 0; i < N; ++i) {
    secure_random_fmpz(arr[i], LEN_AUTH_SHARE);
  }
  std::cout << "Random values assigned" << std::endl;
  raise_sigusr1();
#pragma omp parallel for
  for (int i = 0; i < N - 1; ++i) {
    fmpz_mul(arr[i], arr[i], arr[i + 1]);
  }
  std::cout << "After multiplications" << std::endl;
  raise_sigusr1();
  // now clear and release
  for (int i = 0; i < N; ++i) {
    fmpz_clear(arr[i]);
  }
  std::cout << "After fmpz_clear" << std::endl;
  raise_sigusr1();
  _fmpz_cleanup();
  std::cout << "After fmpz_cleanup" << std::endl;
  raise_sigusr1();
  delete[] arr;
  std::cout << "Array deleted" << std::endl;
  raise_sigusr1();

  std::cout << "Array cleared and released" << std::endl;
  std::vector<ArithLabel> labels(N);
  std::cout << "Vector of " << N << " ArithLabels created." << std::endl;
  raise_sigusr1();
  // now set these labels to some random values
  for (int i = 0; i < N; ++i) {
    secure_random_fmpz(labels[i].auth, LEN_AUTH_SHARE);
    secure_random_fmpz(labels[i].raw, LEN_RAW_SHARE);
  }
  std::cout << "Random values assigned." << std::endl;
  raise_sigusr1();
  // do some operations
  for (int i = 0; i < N - 1; ++i) {
    labels[i] = labels[i] + labels[i + 1];
    labels[i + 1] = labels[i + 1] + labels[i + 1] - labels[i];
  }

  // clear and release
  clear_and_release(labels);
  std::cout << "Vector cleared and released." << std::endl;
  raise_sigusr1();
  _fmpz_cleanup();
  std::cout << "After _fmpz_cleanup in test." << std::endl;
  raise_sigusr1();
}

TEST(TestPaillier, powm_simple) {
  // test 3^11 mod 17 == 7
  fmpz_t base, exp, mod, result;
  fmpz_init_set_ui(base, 3);
  fmpz_init_set_ui(exp, 11);
  fmpz_init_set_ui(mod, 17);
  fmpz_init(result);
  powm_cpu(result, base, exp, mod);
  ASSERT_EQ(fmpz_get_ui(result), 7);
  fmpz_clear(base);
  fmpz_clear(exp);
  fmpz_clear(mod);
  fmpz_clear(result);
}

void gen_random_base_mod(fmpz_t base, fmpz_t mod, int bits = 2048) {
  // 2048-bit mod
  secure_random_fmpz(mod, bits);
  // make mod odd
  if (fmpz_is_even(mod)) fmpz_add_ui(mod, mod, 1);
  // base in [1, mod-1]
  do {
    secure_random_fmpz(base, bits);
    fmpz_mod(base, base, mod);
  } while (fmpz_sgn(base) <= 0 || fmpz_cmp(base, mod) >= 0);
}

void gen_random_exp(fmpz_t exp, int bits = 256) {
  // exp in [1, 2^256-1]
  do {
    secure_random_fmpz(exp, bits);
  } while (fmpz_sgn(exp) <= 0);
}

TEST(TestPaillier, powm_ddh) {
  // sample large random base, mod, and two exponents a,b
  // check that (base^a)^b == (base^b)^a mod mod
  fmpz_t base, mod, a, b, res1, res2;
  fmpz_init(base);
  fmpz_init(mod);
  fmpz_init(a);
  fmpz_init(b);
  fmpz_init(res1);
  fmpz_init(res2);
  gen_random_base_mod(base, mod);
  gen_random_exp(a);
  gen_random_exp(b);

  // compute res1 = (base^a)^b mod mod
  powm_cpu(res1, base, a, mod);
  powm_cpu(res1, res1, b, mod);
  // compute res2 = (base^b)^a mod mod
  powm_cpu(res2, base, b, mod);
  powm_cpu(res2, res2, a, mod);
  // compare
  ASSERT_EQ(fmpz_cmp(res1, res2), 0);
}

TEST(TestPaillier, powm_precomputed) {
  fmpz_t base, exp, mod, result, result_ref;
  fmpz_init(base);
  fmpz_init(exp);
  fmpz_init(mod);
  fmpz_init(result);
  fmpz_init(result_ref);
  gen_random_base_mod(base, mod);
  gen_random_exp(exp);
  // set exp to 1
  // compute reference result
  powm_cpu(result_ref, base, exp, mod);
  // precompute table
  PowmPrecomputeTable table;  // 4 MB
  table.init(base, mod, (size_t)fmpz_bits(exp), 4);
  // compute with precomputed table
  powm_precomputed_cpu(result, exp, table);
  // compare
  ASSERT_EQ(fmpz_cmp(result, result_ref), 0);
  fmpz_clear(base);
  fmpz_clear(exp);
  fmpz_clear(mod);
  fmpz_clear(result);
  fmpz_clear(result_ref);
}

TEST(TestPaillier, bench_powm_precomputed) {
  static constexpr int round = 128;
  fmpz_t base, exp[round], result[round];
  fmpz_init(base);
  secure_random_fmpz(base, LEN_N);
  for (int i = 0; i < round; i++) {
    fmpz_init(exp[i]);
    fmpz_init(result[i]);

    secure_random_fmpz(exp[i], LEN_AUTH_SHARE);
  }
  // first without precomputation
  double start = omp_get_wtime();
#pragma omp parallel for schedule(static)
  for (int i = 0; i < round; i++) {
    powm_cpu(result[i], base, exp[i], default_pk.N2);
  }
  double end = omp_get_wtime();
  double time_no_precompute = end - start;
  printf("powm without precomputation for %d ops: %f s, avg %f ms/op\n", round,
         time_no_precompute, (time_no_precompute / round) * 1000.0);
  // now with precomputation
  int min_window_size = 1;
  int max_window_size = 12;
  for (int window_size = min_window_size; window_size <= max_window_size;
       window_size++) {
    // precompute table for max exp bits
    PowmPrecomputeTable table;
    start = omp_get_wtime();
    table.init_with_window_size(base, default_pk.N2, LEN_AUTH_SHARE,
                                window_size);
    end = omp_get_wtime();
    double time_precompute_init = end - start;
    printf("Precomputation time (window size %d): %f s\n", window_size,
           time_precompute_init);
    start = omp_get_wtime();
#pragma omp parallel for schedule(static)
    for (int i = 0; i < round; i++) {
      powm_precomputed_cpu(result[i], exp[i], table);
    }
    end = omp_get_wtime();
    double time_precompute = end - start;
    printf(
        "powm with precomputation (window size %d) for %d ops: %f s, avg %f "
        "ms/op\n",
        window_size, round, time_precompute,
        (time_precompute / round) * 1000.0);
    printf("Speedup: %fx\n", time_no_precompute / time_precompute);
  }
  fmpz_clear(base);
  for (int i = 0; i < round; i++) {
    fmpz_clear(exp[i]);
    fmpz_clear(result[i]);
  }
}

TEST(TestPaillier, test_batch) {
  int batch_size = 16;
  fmpz_t base, mod;
  fmpz_init(base);
  fmpz_init(mod);
  gen_random_base_mod(base, mod);

  fmpz_t *exps = new fmpz_t[batch_size];
  fmpz_t *results_ref = new fmpz_t[batch_size];
  fmpz_t *results = new fmpz_t[batch_size];
  for (int i = 0; i < batch_size; i++) {
    fmpz_init(exps[i]);
    fmpz_init(results_ref[i]);
    fmpz_init(results[i]);
    gen_random_exp(exps[i]);
  }

  // compute reference results
  powm_cpu_batch(results_ref, base, exps, mod, batch_size);

  // precompute table with max exp bits
  size_t max_exp_bits = 0;
  for (int i = 0; i < batch_size; i++) {
    size_t bits = (size_t)fmpz_bits(exps[i]);
    if (bits > max_exp_bits) max_exp_bits = bits;
  }
  PowmPrecomputeTable table;  // 4 MB
  table.init(base, mod, max_exp_bits, 4);

  // compute with precomputed table
  powm_precomputed_cpu_batch(results, exps, table, batch_size);

  // compare
  for (int i = 0; i < batch_size; i++) {
    ASSERT_EQ(fmpz_cmp(results[i], results_ref[i]), 0);
    fmpz_clear(exps[i]);
    fmpz_clear(results_ref[i]);
    fmpz_clear(results[i]);
  }
  delete[] exps;
  delete[] results_ref;
  delete[] results;
  fmpz_clear(base);
  fmpz_clear(mod);
}

TEST(TestPaillier, bench_powm) {
  int mod_bits = LEN_N;
  int exp_bits = 128;
  int batch_size = 1024;
  fmpz_t base, mod;  // reuse base, mod for all

  // array of exponents and results
  fmpz_t *exps = new fmpz_t[batch_size];
  fmpz_t *exps_warmup = new fmpz_t[batch_size];
  fmpz_t *results = new fmpz_t[batch_size];
  fmpz_t *results_precompute = new fmpz_t[batch_size];

  fmpz_init(base);
  fmpz_init(mod);
  gen_random_base_mod(base, mod, mod_bits);
  for (int i = 0; i < batch_size; i++) {
    fmpz_init(exps[i]);
    fmpz_init(exps_warmup[i]);
    fmpz_init(results[i]);
    fmpz_init(results_precompute[i]);
    gen_random_exp(exps[i], exp_bits);
    fmpz_set(exps_warmup[i], exps[i]);
  }
  // warm up
  powm_cpu_batch(results, base, exps_warmup, mod, batch_size);

  // use omp_get_wtime() for timing
  double start, end;
  // time normal batch powm
  start = omp_get_wtime();
  powm_cpu_batch(results, base, exps, mod, batch_size);
  end = omp_get_wtime();
  double time_normal = end - start;
  printf("Normal batch powm time for %d ops: %f s, avg %f ms/op\n", batch_size,
         time_normal, (time_normal / batch_size) * 1000.0);
  // precompute table
  PowmPrecomputeTable table;
  // time precomputation
  start = omp_get_wtime();
  table.init(base, mod, (size_t)exp_bits, 128);  // 128 MB
  end = omp_get_wtime();
  double time_precompute = end - start;
  printf("Precomputation time: %f s\n", time_precompute);
  // warm up
  powm_precomputed_cpu_batch(results_precompute, exps_warmup, table,
                             batch_size);
  // time precomputed batch powm
  start = omp_get_wtime();
  powm_precomputed_cpu_batch(results_precompute, exps, table, batch_size);

  end = omp_get_wtime();
  double time_precomputed = end - start;
  printf("Precomputed batch powm time for %d ops: %f s, avg %f ms/op\n",
         batch_size, time_precomputed,
         (time_precomputed / batch_size) * 1000.0);
  printf("Speedup: %fx\n", time_normal / time_precomputed);
  // verify results
  for (int i = 0; i < batch_size; i++) {
    ASSERT_EQ(fmpz_cmp(results[i], results_precompute[i]), 0);
    fmpz_clear(exps[i]);
    fmpz_clear(exps_warmup[i]);
    fmpz_clear(results[i]);
    fmpz_clear(results_precompute[i]);
  }
  delete[] exps;
  delete[] exps_warmup;
  delete[] results;
  delete[] results_precompute;
  fmpz_clear(base);
  fmpz_clear(mod);
}

TEST(TestPaillier, bench_prime_prob_gen) {
  double start, end;
  int bits = 1536;
  fmpz_t p;
  fmpz_init(p);
  start = omp_get_wtime();
  sample_prime_prob(p, bits);
  end = omp_get_wtime();
  double time_prime = end - start;
  printf("Generated %d-bit prime in %f s\n", bits, time_prime);
  fmpz_clear(p);
}

TEST(TestPaillier, bench_prime_gen) {
  // skip
  GTEST_SKIP();
  double start, end;
  int bits = LEN_PRIME;
  fmpz_t p;
  fmpz_init(p);
  start = omp_get_wtime();
  sample_prime(p, bits);
  end = omp_get_wtime();
  double time_prime = end - start;
  printf("Generated %d-bit prime in %f s\n", bits, time_prime);
  fmpz_clear(p);
}

TEST(TestPaillier, generator_pow) {
  fmpz_t exp, result, result_ref;
  fmpz_init(exp);
  fmpz_init(result);
  fmpz_init(result_ref);
  paillier_generator_pow(result, exp, default_sk);
  powm(result_ref, default_pk.g, exp, default_pk.N2);
  ASSERT_EQ(fmpz_cmp(result, result_ref), 0);
  fmpz_clear(exp);
  fmpz_clear(result);
  fmpz_clear(result_ref);
}

TEST(TestPaillier, encrypt_decrypt_correctness) {
  // Generate keypair
  double start = omp_get_wtime();
  paillier_keygen(default_pk, default_sk, 1024);  // 256 MB precompute tables
  double end = omp_get_wtime();
  printf("Paillier keygen time: %f s\n", end - start);

  // Test with several different message sizes
  int test_messages[] = {0, 1, 42, 1000000, (1 << 20), (1 << 30)};
  int num_tests = sizeof(test_messages) / sizeof(test_messages[0]);

  for (int i = 0; i < num_tests; i++) {
    fmpz_t message, decrypted;
    fmpz_init(message);
    fmpz_init(decrypted);

    fmpz_set_ui(message, test_messages[i]);

    // Encrypt
    PaillierCiphertext ct;
    start = omp_get_wtime();
    paillier_encrypt_with_sk(ct, message, default_pk, default_sk);
    end = omp_get_wtime();
    printf("Encrypt message %d time: %f ms\n", test_messages[i],
           (end - start) * 1000.0);

    // Decrypt
    start = omp_get_wtime();
    paillier_decrypt(decrypted, ct, default_pk, default_sk);
    end = omp_get_wtime();
    printf("Decrypt message %d time: %f ms\n", test_messages[i],
           (end - start) * 1000.0);

    // Verify correctness
    ASSERT_EQ(fmpz_cmp(message, decrypted), 0)
        << "Decryption failed for message: " << test_messages[i];

    fmpz_clear(message);
    fmpz_clear(decrypted);
  }

  // Test with a random large message
  fmpz_t large_message, large_decrypted;
  fmpz_init(large_message);
  fmpz_init(large_decrypted);

  // Generate random message smaller than N (to avoid overflow)
  secure_random_fmpz(large_message,
                     LEN_ARITH_DIGIT);  // Use arithmetic digit size

  // Encrypt and decrypt
  PaillierCiphertext large_ct;
  paillier_encrypt_with_sk(large_ct, large_message, default_pk, default_sk);
  paillier_decrypt(large_decrypted, large_ct, default_pk, default_sk);

  // Verify correctness
  ASSERT_EQ(fmpz_cmp(large_message, large_decrypted), 0)
      << "Decryption failed for large random message";

  printf("Encrypt/decrypt test passed for large random message\n");

  fmpz_clear(large_message);
  fmpz_clear(large_decrypted);
}

TEST(TestPaillier, secret_share) {
  paillier_keygen(default_pk, default_sk, 256);  // 256 MB precompute tables

  fmpz_t secret;
  fmpz_init(secret);
  secure_random_fmpz(secret, LEN_ARITH_DIGIT);

  ArithLabel share1, share2;
  secret_share(share1, share2, secret, default_sk);

  // merge shares
  fmpz_t merge1, merge2;
  fmpz_init(merge1);
  fmpz_init(merge2);
  share1.merge(merge1);
  share2.merge(merge2);

  // split merged shares
  ArithLabel split1, split2;
  split1.split(merge1);
  split2.split(merge2);

  fmpz_t reconstructed;
  fmpz_init(reconstructed);
  open_share(reconstructed, split1, split2);

  // check reconstructed == secret
  ASSERT_EQ(fmpz_cmp(reconstructed, secret), 0);

  fmpz_clear(secret);
  fmpz_clear(reconstructed);
  fmpz_clear(merge1);
  fmpz_clear(merge2);
}

TEST(TestPaillier, test_hss) {
  // test homomorphic secret-sharing
  // Let A and B be two LEN_RAW_SHARE / 2 bit integers
  int round = 100;
  for (int i = 0; i < round; i++) {
    fmpz_t A, B, A_times_B;
    fmpz_init(A);
    fmpz_init(B);
    fmpz_init(A_times_B);
    secure_random_fmpz(A, LEN_ARITH_DIGIT / 2);
    secure_random_fmpz(B, LEN_ARITH_DIGIT / 2);
    // compute A * B
    fmpz_mul(A_times_B, A, B);

    // encrypt A
    PaillierCiphertext enc_A;
    paillier_encrypt_with_sk(enc_A, A, default_pk, default_sk);

    // secret-share B into B1, B2
    ArithLabel Bshares[2];
    fmpz_t ABshares[2];
    secret_share(Bshares[0], Bshares[1], B, default_sk);
    // compute x[i] = enc_A.ct[0]^(Bshares[i].auth) * enc_A.ct[1]^Bshares[i].raw
    fmpz_t x;
    fmpz_init(x);

    for (int party = 0; party < 2; ++party) {
      fmpz_t powm_result[2];
      fmpz_init(powm_result[0]);
      fmpz_init(powm_result[1]);
      powm(powm_result[0], enc_A.ct[0], Bshares[party].auth, default_pk.N2);
      powm(powm_result[1], enc_A.ct[1], Bshares[party].raw, default_pk.N2);
      fmpz_mod_mul(x, powm_result[0], powm_result[1], &default_pk.mod_N2_ctx);
      ddlog(ABshares[party], x, default_pk);

      fmpz_clear(powm_result[0]);
      fmpz_clear(powm_result[1]);
    }
    fmpz_t revealed_value;
    fmpz_init(revealed_value);
    // open the result from ABshares
    fmpz_sub(revealed_value, ABshares[1], ABshares[0]);
    // check revealed_value == A * B
    ASSERT_EQ(fmpz_cmp(revealed_value, A_times_B), 0);
    fmpz_clear(A);
    fmpz_clear(B);
    fmpz_clear(x);
    fmpz_clear(A_times_B);
  }
}

TEST(TestPaillier, bench_ddlog) {
  // generate a batch of random x of size LEN_N2

  int batch_size = 1000;
  fmpz_t *x_array = new fmpz_t[batch_size];
  for (int i = 0; i < batch_size; i++) {
    fmpz_init(x_array[i]);
    secure_random_fmpz(x_array[i], LEN_N2);
  }
  fmpz_t *results = new fmpz_t[batch_size];
  for (int i = 0; i < batch_size; i++) {
    fmpz_init(results[i]);
  }
  double start = omp_get_wtime();
  for (int i = 0; i < batch_size; i++) {
    ddlog(results[i], x_array[i], default_pk);
  }
  double end = omp_get_wtime();
  double time_ddlog = end - start;
  printf("ddlog time for %d ops: %f s, avg %f ms/op\n", batch_size, time_ddlog,
         (time_ddlog / batch_size) * 1000.0);
  // now use omp parallel for to speed up
  start = omp_get_wtime();
#pragma omp parallel for schedule(static)
  for (int i = 0; i < batch_size; i++) {
    ddlog(results[i], x_array[i], default_pk);
  }
  end = omp_get_wtime();
  time_ddlog = end - start;
  printf("ddlog (omp) time for %d ops: %f s, avg %f ms/op\n", batch_size,
         time_ddlog, (time_ddlog / batch_size) * 1000.0);
  for (int i = 0; i < batch_size; i++) {
    fmpz_clear(x_array[i]);
    fmpz_clear(results[i]);
  }
}

TEST(TestPaillier, serialize_deserialize_ciphertext) {
  // Test with a few different messages
  int test_messages[] = {0, 42, 1000000};
  int num_tests = sizeof(test_messages) / sizeof(test_messages[0]);

  for (int i = 0; i < num_tests; i++) {
    fmpz_t message;
    fmpz_init(message);
    fmpz_set_ui(message, test_messages[i]);

    // Encrypt to get original ciphertext
    PaillierCiphertext original_ct;
    paillier_encrypt_with_sk(original_ct, message, default_pk, default_sk);

    // Serialize and deserialize
    uint8_t buffer[PAILLIER_CIPHER_TEXT_SIZE];
    original_ct.serialize(buffer);

    PaillierCiphertext deserialized_ct;
    deserialized_ct.deserialize(buffer);

    // Verify both ciphertext components are identical
    ASSERT_EQ(fmpz_cmp(original_ct.ct[0], deserialized_ct.ct[0]), 0);
    ASSERT_EQ(fmpz_cmp(original_ct.ct[1], deserialized_ct.ct[1]), 0);

    fmpz_clear(message);
  }
}

TEST(TestPaillier, EncryptDecrypt) {
  // Test encryption and decryption correctness
  PaillierPubKey pk;
  PaillierPrivKey sk;
  paillier_keygen(pk, sk);

  for (int i = 0; i < 100; i++) {
    fmpz_t m, m_dec;
    fmpz_init(m);
    fmpz_init(m_dec);

    secure_random_fmpz(m, 128);  // Random message

    // Encrypt
    PaillierCiphertext ct;
    paillier_encrypt_with_sk(ct, m, pk, sk);

    // Decrypt
    paillier_decrypt(m_dec, ct, pk, sk);

    // Check if decrypted message matches original
    ASSERT_EQ(fmpz_cmp(m, m_dec), 0);

    fmpz_clear(m);
    fmpz_clear(m_dec);
  }
}

TEST(TestPaillier, test_encrypt_batch) {
  PaillierPubKey pk;
  PaillierPrivKey sk;
  paillier_keygen(pk, sk);

  const int batch_size = 16;
  PaillierCiphertext ct_single[batch_size];
  PaillierCiphertext ct_batch[batch_size];
  fmpz_t m[batch_size];
  fmpz_t r[batch_size];

  for (int i = 0; i < batch_size; i++) {
    fmpz_init(m[i]);
    fmpz_init(r[i]);
    secure_random_fmpz(m[i], 128);
    secure_random_fmpz(r[i], LEN_K);
  }

  // Batch encryption
  paillier_encrypt_with_sk_batch(ct_batch, m, pk, sk, r, batch_size);

  // Single encryption
  for (int i = 0; i < batch_size; i++) {
    paillier_encrypt_with_sk(ct_single[i], m[i], pk, sk, r[i]);
  }

  // Compare results
  for (int i = 0; i < batch_size; i++) {
    EXPECT_EQ(fmpz_equal(ct_single[i].ct[0], ct_batch[i].ct[0]), 1);
    EXPECT_EQ(fmpz_equal(ct_single[i].ct[1], ct_batch[i].ct[1]), 1);
  }

  for (int i = 0; i < batch_size; i++) {
    fmpz_clear(m[i]);
    fmpz_clear(r[i]);
  }
}

TEST(TestPaillier, bench_encrypt) {
  PaillierPubKey pk;
  PaillierPrivKey sk;
  paillier_keygen(pk, sk);

  const int batch_size = 128;
  PaillierCiphertext ct[batch_size];
  fmpz_t m[batch_size];
  fmpz_t r[batch_size];

  for (int i = 0; i < batch_size; i++) {
    fmpz_init(m[i]);
    fmpz_init(r[i]);
    secure_random_fmpz(m[i], 128);
    secure_random_fmpz(r[i], LEN_K);
  }

  // Benchmark single encryption
  double start_single = omp_get_wtime();
  for (int i = 0; i < batch_size; i++) {
    paillier_encrypt_with_sk(ct[i], m[i], pk, sk, r[i]);
  }
  double end_single = omp_get_wtime();
  double single_duration = end_single - start_single;
  printf("Single encryption time for %d items: %fs, average: %fms\n",
         batch_size, single_duration, (single_duration / batch_size) * 1000);

  // Benchmark batch encryption
  double start_batch = omp_get_wtime();
  paillier_encrypt_with_sk_batch(ct, m, pk, sk, r, batch_size);
  double end_batch = omp_get_wtime();
  double batch_duration = end_batch - start_batch;
  printf("Batch encryption time for %d items: %fs, average: %fms\n", batch_size,
         batch_duration, (batch_duration / batch_size) * 1000);

  printf("Speedup: %fx\n", single_duration / batch_duration);

  for (int i = 0; i < batch_size; i++) {
    fmpz_clear(m[i]);
    fmpz_clear(r[i]);
  }
}