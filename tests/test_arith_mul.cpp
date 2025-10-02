#include <gtest/gtest.h>
#include <omp.h>

#include "arith_mul.hpp"
using namespace ZebraGRAM;

TEST(TestArith, garble) {
  int round = 100;
  for (int r = 0; r < round; ++r) {
    printf("round %d\n", r);
    fmpz_t a, b, c;
    fmpz_init(a);
    fmpz_init(b);
    fmpz_init(c);
    secure_random_fmpz(a, LEN_ARITH_DIGIT / 2);
    secure_random_fmpz(b, LEN_ARITH_DIGIT / 2);
    ArithLabel a_share[2];
    ArithLabel b_share[2];
    secret_share(a_share[0], a_share[1], a, default_sk);
    secret_share(b_share[0], b_share[1], b, default_sk);
    uint8_t material[GARBLE_MULT_SIZE];
    ArithLabel c_share[2];
    garble_mul(c_share[0], material, a_share[0], b_share[0], default_pk,
               default_sk);
    // turn material back to ciphertexts
    PaillierCiphertext ct[2];
    ct[0].deserialize(material);
    ct[1].deserialize(material + PAILLIER_CIPHER_TEXT_SIZE);

    fmpz_t powm_result[4];
    for (int i = 0; i < 4; i++) {
      fmpz_init(powm_result[i]);
    }
    powm(powm_result[0], ct[0].ct[0], b_share[0].auth, default_pk.N2);
    powm(powm_result[1], ct[0].ct[1], b_share[0].raw, default_pk.N2);
    powm(powm_result[2], ct[1].ct[0], a_share[0].auth, default_pk.N2);
    powm(powm_result[3], ct[1].ct[1], a_share[0].raw, default_pk.N2);
    // combine results to powm_result[0]
    for (int i = 1; i < 4; i++) {
      fmpz_mod_mul(powm_result[0], powm_result[0], powm_result[i],
                   &default_pk.mod_N2_ctx);
    }
    // store ddlog result in powm_result[1]
    ddlog(powm_result[1], powm_result[0], default_pk);
    // store a.raw * b.combine in powm_result[2]
    b_share[0].merge(powm_result[2]);
    fmpz_mul(powm_result[2], a_share[0].raw, powm_result[2]);
    // powm_result[3] = powm_result[2] - powm_result[1]
    fmpz_sub(powm_result[3], powm_result[2], powm_result[1]);
    fmpz_mod_set_fmpz(powm_result[3], powm_result[3], &default_pk.mod_N_ctx);

    ArithLabel c_share_ref;
    c_share_ref.split(powm_result[3]);
    // compare c_share[0] and c_share_ref
    ASSERT_EQ(fmpz_cmp(c_share[0].raw, c_share_ref.raw), 0);
    ASSERT_EQ(fmpz_cmp(c_share[0].auth, c_share_ref.auth), 0);

    for (int i = 0; i < 4; i++) {
      fmpz_clear(powm_result[i]);
    }
    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(c);
  }
}

TEST(TestArith, mul) {
  int round = 100;
  for (int r = 0; r < round; ++r) {
    printf("round %d\n", r);
    fmpz_t a, b, c, ab;
    fmpz_init(a);
    fmpz_init(b);
    fmpz_init(c);
    fmpz_init(ab);
    secure_random_fmpz(a, LEN_ARITH_DIGIT / 2);
    secure_random_fmpz(b, LEN_ARITH_DIGIT / 2);
    fmpz_mul(ab, a, b);
    ArithLabel a_share[2];
    ArithLabel b_share[2];
    secret_share(a_share[0], a_share[1], a, default_sk);
    secret_share(b_share[0], b_share[1], b, default_sk);
    uint8_t material[GARBLE_MULT_SIZE];
    ArithLabel c_share[2];
    garble_mul(c_share[0], material, a_share[0], b_share[0], default_pk,
               default_sk);
    eval_mul(c_share[1], material, a_share[1], b_share[1], default_pk);
    fmpz_t c_merge[2];
    fmpz_init(c_merge[0]);
    fmpz_init(c_merge[1]);
    c_share[0].merge(c_merge[0]);
    c_share[1].merge(c_merge[1]);
    // check ab * ((K << LEN_RAW_SHARE) + 1) == c_merge[1] - c_merge[0] mod N
    fmpz_t check_left, check_right;
    fmpz_init(check_left);
    fmpz_init(check_right);
    fmpz_mul_2exp(check_left, default_sk.K, LEN_RAW_SHARE);
    fmpz_add_ui(check_left, check_left, 1);
    fmpz_sub(check_right, c_merge[1], c_merge[0]);
    fmpz_mod_set_fmpz(check_right, check_right, &default_pk.mod_N_ctx);

    fmpz_t c_revealed;
    fmpz_init(c_revealed);
    open_share(c_revealed, c_share[0], c_share[1]);
    ASSERT_EQ(fmpz_cmp(c_revealed, ab), 0);
    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(c);
    fmpz_clear(ab);
    fmpz_clear(c_revealed);
  }
}

TEST(TestArith, bench_garble_eval) {
  constexpr int round = 100;
  fmpz_t a[round], b[round];
  ArithLabel a_share[round][2];
  ArithLabel b_share[round][2];
  uint8_t material[round][GARBLE_MULT_SIZE];
  ArithLabel c_share[round][2];
  for (int r = 0; r < round; ++r) {
    fmpz_init(a[r]);
    fmpz_init(b[r]);
    secure_random_fmpz(a[r], LEN_ARITH_DIGIT / 2);
    secure_random_fmpz(b[r], LEN_ARITH_DIGIT / 2);
    secret_share(a_share[r][0], a_share[r][1], a[r], default_sk);
    secret_share(b_share[r][0], b_share[r][1], b[r], default_sk);
  }

  double start = omp_get_wtime();
  for (int r = 0; r < round; ++r) {
    garble_mul(c_share[r][0], material[r], a_share[r][0], b_share[r][0],
               default_pk, default_sk);
  }
  double garble_time = omp_get_wtime() - start;
  printf("garble time for %d mul: %f s, avg %f ms/op\n", round, garble_time,
         (garble_time / round) * 1000.0);

  start = omp_get_wtime();
  for (int r = 0; r < round; ++r) {
    eval_mul(c_share[r][1], material[r], a_share[r][1], b_share[r][1],
             default_pk);
  }
  double eval_time = omp_get_wtime() - start;
  printf("eval time for %d mul: %f s, avg %f ms/op\n", round, eval_time,
         (eval_time / round) * 1000.0);
  for (int r = 0; r < round; ++r) {
    fmpz_clear(a[r]);
    fmpz_clear(b[r]);
  }
}

TEST(TestArith, batch_mul_correctness) {
  // set thread number as 1
  constexpr int vec_size = 20;
  int round = 10;

  for (int r = 0; r < round; ++r) {
    printf("batch round %d\n", r);

    // Generate random values
    fmpz_t a, b[vec_size], ab[vec_size];
    fmpz_init(a);
    for (int i = 0; i < vec_size; i++) {
      fmpz_init(b[i]);
      fmpz_init(ab[i]);
    }

    secure_random_fmpz(a, 1);
    for (int i = 0; i < vec_size; i++) {
      secure_random_fmpz(b[i], LEN_ARITH_DIGIT);
      fmpz_mul(ab[i], a, b[i]);
    }

    // Create shares
    ArithLabel a_share[2];
    ArithLabel b_share_vec0[vec_size], b_share_vec1[vec_size];
    ArithLabel c_share_vec0[vec_size], c_share_vec1[vec_size];
    secret_share(a_share[0], a_share[1], a, 1, default_sk);
    for (int i = 0; i < vec_size; i++) {
      secret_share(b_share_vec0[i], b_share_vec1[i], b[i], default_sk);
    }

    // Prepare batch inputs

    // Batch garble
    size_t material_size = (vec_size + 1) * PAILLIER_CIPHER_TEXT_SIZE;
    uint8_t *material = new uint8_t[material_size];

    size_t actual_size =
        garble_mul_vec(c_share_vec0, material, a_share[0], b_share_vec0,
                       vec_size, default_pk, default_sk);
    ASSERT_EQ(actual_size, material_size);

    eval_mul_vec(c_share_vec1, material, a_share[1], b_share_vec1, vec_size,
                 default_pk);

    // Verify results
    for (int i = 0; i < vec_size; i++) {
      fmpz_t c_revealed;
      fmpz_init(c_revealed);
      open_share(c_revealed, c_share_vec0[i], c_share_vec1[i]);
      ASSERT_EQ(fmpz_cmp(c_revealed, ab[i]), 0)
          << "Multiplication " << i << " failed in round " << r;
      fmpz_clear(c_revealed);
    }

    // Cleanup
    delete[] material;
    fmpz_clear(a);
    for (int i = 0; i < vec_size; i++) {
      fmpz_clear(b[i]);
      fmpz_clear(ab[i]);
    }
  }

  printf("Batch multiplication test passed\n");
}

TEST(TestArith, bench_batch_mul) {
  // only test the batch garble for now
  constexpr int vec_size = 100;
  int round = 10;
  double total_garble_time = 0.0;
  double total_eval_time = 0.0;
  // rerun key generation
  PaillierPrivKey sk;
  PaillierPubKey pk;
  paillier_keygen(pk, sk, 1024);

  for (int r = 0; r < round; ++r) {
    // Generate random values
    fmpz_t a;
    fmpz_t b[vec_size];
    fmpz_init(a);
    for (int i = 0; i < vec_size; i++) {
      fmpz_init(b[i]);
    }

    secure_random_fmpz(a, LEN_ARITH_DIGIT / 2);
    for (int i = 0; i < vec_size; i++) {
      secure_random_fmpz(b[i], LEN_ARITH_DIGIT / 2);
    }

    // Create shares
    ArithLabel a_share[2];
    ArithLabel b_share_vec0[vec_size];
    ArithLabel c_share_vec0[vec_size];
    ArithLabel dummy;

    secret_share(a_share[0], a_share[1], a, 1, sk);
    for (int i = 0; i < vec_size; i++) {
      secret_share(b_share_vec0[i], dummy, b[i], sk);
    }

    // Allocate material
    size_t material_size = (vec_size + 1) * PAILLIER_CIPHER_TEXT_SIZE;
    uint8_t *material = new uint8_t[material_size];

    // Benchmark batch garble
    double start = omp_get_wtime();
    garble_mul_vec(c_share_vec0, material, a_share[0], b_share_vec0, vec_size,
                   pk, sk);
    double end = omp_get_wtime();
    eval_mul_vec(c_share_vec0, material, a_share[1], b_share_vec0, vec_size,
                 pk);
    double eval_end = omp_get_wtime();
    total_eval_time += (eval_end - end);

    total_garble_time += (end - start);
    total_eval_time += (eval_end - end);

    // Cleanup
    delete[] material;
    fmpz_clear(a);
    for (int i = 0; i < vec_size; i++) {
      fmpz_clear(b[i]);
    }
  }

  double avg_garble_time = total_garble_time / round;
  double avg_garble_per_op =
      (total_garble_time / (round * vec_size)) * 1000.0;  // ms per operation

  printf("Batch garble time for %d rounds of %d multiplications: %f s\n", round,
         vec_size, total_garble_time);
  printf("Average batch garble time per round: %f s\n", avg_garble_time);
  printf("Average time per multiplication: %f ms\n", avg_garble_per_op);

  double avg_eval_time = total_eval_time / round;
  double avg_eval_per_op =
      (total_eval_time / (round * vec_size)) * 1000.0;  // ms per operation
  printf("Batch eval time for %d rounds of %d multiplications: %f s\n", round,
         vec_size, total_eval_time);
  printf("Average batch eval time per round: %f s\n", avg_eval_time);
  printf("Average eval time per multiplication: %f ms\n", avg_eval_per_op);
}