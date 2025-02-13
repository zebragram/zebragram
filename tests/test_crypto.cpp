#include <gtest/gtest.h>
#include <omp.h>

#include "bit.hpp"
#include "ec.hpp"
#include "hash.hpp"
#include "three_halves.hpp"
using namespace PicoGRAM;

TEST(TestCrypto, hash) {
  Label label = Label::random();
  Label hash = Hash(label);
  Label hash2 = Hash(label);
  ASSERT_EQ(hash, hash2);
}

TEST(TestCrypto, three_halves) {
  for (int r = 0; r < 10000; ++r) {
    Label A = Label::random();
    Label B = Label::random();
    uint64_t nonce = rand();
    uint8_t material[25];
    Label C = ThreeHalves::Garble(A, B, nonce, material);
    Label C_eval = ThreeHalves::Eval(A, B, nonce, material);
    ASSERT_EQ(C, C_eval);
    C_eval = ThreeHalves::Eval(A ^ key_manager.get_Delta(), B, nonce, material);
    ASSERT_EQ(C, C_eval);
    C_eval = ThreeHalves::Eval(A, B ^ key_manager.get_Delta(), nonce, material);
    ASSERT_EQ(C, C_eval);
    C_eval = ThreeHalves::Eval(A ^ key_manager.get_Delta(),
                               B ^ key_manager.get_Delta(), nonce, material);
    ASSERT_EQ(C ^ key_manager.get_Delta(), C_eval);
  }
}

TEST(TestCrypto, three_halves_perf) {
  Label A = Label::random();
  Label B = Label::random();
  uint64_t nonce = rand();
  uint8_t material[25];
  Label C;
  int round = 100000;
  Timer timer;
  timer.tic();
  for (int r = 0; r < round; ++r) {
    C = ThreeHalves::Garble(A, B, nonce, material);
  }
  timer.toc();
  std::cout << "Garble time: " << timer.get_time() / round * 1e6 << " us\n";
  timer.reset();
  timer.tic();
  for (int r = 0; r < round; ++r) {
    C = ThreeHalves::Eval(A, B, nonce, material);
  }
  timer.toc();
  std::cout << "Eval time: " << timer.get_time() / round * 1e6 << " us\n";
}

TEST(TestCrypto, xor128bit_perf) {
  __m128i a = _mm_set_epi64x(rand(), rand());
  __m128i b = _mm_set_epi64x(rand(), rand());
  int round = 1000000;
  Timer timer;
  timer.tic();
  for (int r = 0; r < round; ++r) {
    a = _mm_xor_si128(a, b);
  }
  timer.toc();
  std::cout << "a = " << a[0] << " " << a[1] << std::endl;
  std::cout << "xor128bit time: " << timer.get_time() / round * 1e6 << " us\n";
}

TEST(TestCrypto, bigint_encode_decode) {
  for (int r = 0; r < 10000; ++r) {
    BigInt a = BigInt::random();
    uint8_t bytes[32];
    a.to_bytes(bytes);
    BigInt b(bytes);
    ASSERT_EQ(a, b);
  }
}

TEST(TestCrypto, bigint_encrypt_decrypt) {
  for (int r = 0; r < 10000; ++r) {
    Label key = Label::random();
    uint64_t nonce1 = rand();
    uint64_t nonce2 = rand();
    BigInt a = BigInt::random();
    uint8_t bytes[BigInt::byte_length];
    a.enc(key, bytes, nonce1, nonce2);
    BigInt b;
    b.dec(key, bytes, nonce1, nonce2);
    ASSERT_EQ(a, b);
    BigInt c;
    uint8_t flipped_bytes[BigInt::byte_length];
    memcpy(flipped_bytes, bytes, BigInt::byte_length);
    flipped_bytes[rand() % BigInt::byte_length] ^= 1 << (rand() % 8);
    c.dec(key, flipped_bytes, nonce1, nonce2);
    ASSERT_NE(a, c);
    BigInt d;
    __m128i flipped_key_m128 =
        _mm_xor_si128(key.get_m128i(), _mm_set_epi64x(0, 1));
    Label flipped_key;
    flipped_key.set_m128i(flipped_key_m128);
    d.dec(flipped_key, bytes, nonce1, nonce2);
    ASSERT_NE(a, d);
    uint64_t nonce1_flipped = nonce1 ^ 1;
    BigInt e;
    e.dec(key, bytes, nonce1_flipped, nonce2);
    ASSERT_NE(a, e);
    uint64_t nonce2_flipped = nonce2 ^ 1;
    BigInt f;
    f.dec(key, bytes, nonce1, nonce2_flipped);
    ASSERT_NE(a, f);
  }
}

TEST(TestCrypto, ec_point_encrypt_decrypt) {
  for (int r = 0; r < 10000; ++r) {
    Label key = Label::random();
    uint64_t nonce1 = rand();
    uint64_t nonce2 = rand();
    ECPoint A = ECPoint::generator_pow(BigInt::random());
    uint8_t bytes[ECPoint::byte_length];
    A.enc(key, bytes, nonce1, nonce2);
    ECPoint B;
    B.dec(key, bytes, nonce1, nonce2);
    ASSERT_EQ(A, B);
  }
}

TEST(TestCrypto, modular) {
  BigInt a = BigInt::random();
  BigInt b = BigInt::random();
  BigInt c = a * b;
  BigInt d = c * a.inv();
  ASSERT_EQ(b, d);
}

TEST(TestCrypto, inv_batch) {
  uint64_t batch_size = 100;
  std::vector<BigInt> a(batch_size);
  std::vector<BigInt> a_inv_ref(batch_size);
  for (uint64_t i = 0; i < batch_size; ++i) {
    a[i] = BigInt::random();
    a_inv_ref[i] = a[i].inv();
  }
  std::vector<BigInt*> a_ptrs(batch_size);
  for (uint64_t i = 0; i < batch_size; ++i) {
    a_ptrs[i] = &a[i];
  }
  BigInt::inv_batch(a_ptrs.data(), batch_size);
  for (uint64_t i = 0; i < batch_size; ++i) {
    ASSERT_EQ(a[i], a_inv_ref[i]);
  }
}

TEST(TestCrypto, exp) {
  BigInt a = BigInt::random();
  BigInt b = BigInt::random();
  ECPoint A = ECPoint::generator_pow(a);
  ECPoint B = ECPoint::generator_pow(b);
  ECPoint C = A.pow(b);
  ECPoint D = B.pow(a);
  ASSERT_EQ(C, D);
}

TEST(TestCrypto, dh) {
  BigInt a = BigInt::random();
  BigInt b = BigInt::random();
  ECPoint A = ECPoint::generator_pow(a);
  ECPoint B = A.pow(b);
  ECPoint C = ECPoint::generator_pow((a * b).from_montgomery());
  ASSERT_EQ(B, C);
}

TEST(TestCrypto, exp_inv) {
  BigInt a = BigInt::random();
  BigInt b = BigInt::random();
  BigInt b_inv = b.inv();
  ECPoint A = ECPoint::generator_pow(a);
  ECPoint B = A.pow(b).pow(b_inv);
  ASSERT_EQ(A, B);
}

TEST(TestCrypto, copy) {
  ECPoint A = ECPoint::generator_pow(BigInt::random());
  ECPoint B = A;
  BigInt a = BigInt::random();
  A.pow(a, A);
  B.pow(a, B);
  ASSERT_EQ(A, B);
}

TEST(TestCrypto, BN_perf) {
  uint64_t r = 1000000;
  Timer timer;
  std::vector<BIGNUM*> rand_a(r);
  std::vector<BIGNUM*> rand_b(r);
  // first test BN_new perf on rand_a
  timer.tic();
  for (uint64_t i = 0; i < r; ++i) {
    rand_a[i] = BN_new();
  }
  timer.toc();
  std::cout << "BN_new time: " << timer.get_time() / r * 1e6 << " us\n";
  timer.reset();
  uint8_t bytes[32];
  for (int i = 0; i < 32; ++i) {
    bytes[i] = rand() % 256;
  }
  // then test BN_bin2bn perf on rand_a
  timer.tic();
  for (uint64_t i = 0; i < r; ++i) {
    BN_bin2bn(bytes, 32, rand_a[i]);
  }
  timer.toc();
  std::cout << "BN_bin2bn time: " << timer.get_time() / r * 1e6 << " us\n";
  timer.reset();
  // then test BN_dup perf by copying rand_a to rand_b
  timer.tic();
  for (uint64_t i = 0; i < r; ++i) {
    rand_b[i] = BN_dup(rand_a[i]);
  }
  timer.toc();
  std::cout << "BN_dup time: " << timer.get_time() / r * 1e6 << " us\n";
  timer.reset();
  // now test BN_copy perf by copying rand_a to rand_b
  timer.tic();
  for (uint64_t i = 0; i < r; ++i) {
    BN_copy(rand_b[i], rand_a[i]);
  }
  timer.toc();
  std::cout << "BN_copy time: " << timer.get_time() / r * 1e6 << " us\n";
  timer.reset();
  // finally test BN_free perf
  timer.tic();
  for (uint64_t i = 0; i < r; ++i) {
    BN_free(rand_a[i]);
    BN_free(rand_b[i]);
  }
  timer.toc();
  std::cout << "BN_free time: " << timer.get_time() / (r * 2) * 1e6 << " us\n";
  timer.reset();

  // now test BN_CTX_new perf
  std::vector<BN_CTX*> ctx(r);
  timer.tic();
  for (uint64_t i = 0; i < r; ++i) {
    ctx[i] = BN_CTX_new();
  }
  timer.toc();
  std::cout << "BN_CTX_new time: " << timer.get_time() / r * 1e6 << " us\n";
  timer.reset();
  // now test BN_CTX_free perf
  timer.tic();
  for (uint64_t i = 0; i < r; ++i) {
    BN_CTX_free(ctx[i]);
  }
  timer.toc();
  std::cout << "BN_CTX_free time: " << timer.get_time() / r * 1e6 << " us\n";
  timer.reset();
}

TEST(TestCrypto, EC_POINT_perf) {
  uint64_t r = 100000;
  Timer timer;
  std::vector<EC_POINT*> rand_a(r);
  std::vector<EC_POINT*> rand_b(r);

  // first test EC_POINT_new perf on rand_a
  timer.tic();
  for (uint64_t i = 0; i < r; ++i) {
    rand_a[i] = EC_POINT_new(ECPoint::group);
    rand_b[i] = EC_POINT_new(ECPoint::group);
  }
  timer.toc();
  std::cout << "EC_POINT_new time: " << timer.get_time() / (2 * r) * 1e6
            << " us\n";
  timer.reset();
  // then test EC_POINT_oct2point perf on rand_a
  ECPoint A = ECPoint::generator_pow(BigInt::random());
  timer.tic();
  for (uint64_t i = 0; i < r; ++i) {
    EC_POINT_oct2point(ECPoint::group, rand_a[i], A.bytes, ECPoint::byte_length,
                       BigInt::global_ctx);
  }
  timer.toc();
  std::cout << "EC_POINT_oct2point time: " << timer.get_time() / r * 1e6
            << " us\n";
  timer.reset();

  // then test EC_POINT_copy perf by copying rand_a to rand_b
  timer.tic();
  for (uint64_t i = 0; i < r; ++i) {
    EC_POINT_copy(rand_b[i], rand_a[i]);
  }
  timer.toc();
  std::cout << "EC_POINT_copy time: " << timer.get_time() / r * 1e6 << " us\n";
  timer.reset();

  // then test EC_POINT_point2oct perf on rand_a
  uint8_t bytes[ECPoint::byte_length];
  timer.tic();
  for (uint64_t i = 0; i < r; ++i) {
    EC_POINT_point2oct(ECPoint::group, rand_a[i], POINT_CONVERSION_COMPRESSED,
                       bytes, ECPoint::byte_length, BigInt::global_ctx);
  }
  timer.toc();
  std::cout << "EC_POINT_point2oct time: " << timer.get_time() / r * 1e6
            << " us\n";
  timer.reset();
  // finally test EC_POINT_free perf
  timer.tic();
  for (uint64_t i = 0; i < r; ++i) {
    EC_POINT_free(rand_a[i]);
    EC_POINT_free(rand_b[i]);
  }
  timer.toc();
  std::cout << "EC_POINT_free time: " << timer.get_time() / (r * 2) * 1e6
            << " us\n";
  timer.reset();
}

#ifdef NDEBUG

TEST(TestCrypto, BigIntMult_perf) {
  uint64_t r = 1000000;
  std::vector<BigInt> rand_a(r);
  std::vector<BigInt> rand_b(r);
  for (uint64_t i = 0; i < r; ++i) {
    rand_a[i] = BigInt::random();
    rand_b[i] = BigInt::random();
  }
  auto start = std::chrono::high_resolution_clock::now();
  for (uint64_t i = 0; i < r; ++i) {
    BigInt c = rand_a[i] * rand_b[i];
  }
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff = end - start;
  std::cout << "BigInt multiply time: " << diff.count() / r * 1e6 << " us\n";
}

TEST(TestCrypto, BigIntMult_repeat_perf) {
  uint64_t r = 1000000;
  BigInt rand_a = BigInt::random();
  std::vector<BigInt> rand_b(r);
  for (uint64_t i = 0; i < r; ++i) {
    rand_b[i] = BigInt::random();
  }
  auto start = std::chrono::high_resolution_clock::now();
  for (uint64_t i = 0; i < r; ++i) {
    rand_a *= rand_b[i];
  }
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff = end - start;
  std::cout << "BigInt multiply time: " << diff.count() / r * 1e6 << " us\n";
}

TEST(TestCrypto, BigInt_convert_montgomery_perf) {
  uint64_t r = 1000000;
  std::vector<BigInt> rand_b(r);
  for (uint64_t i = 0; i < r; ++i) {
    rand_b[i] = BigInt::random();
  }
  auto start = std::chrono::high_resolution_clock::now();
  for (uint64_t i = 0; i < r; ++i) {
    rand_b[i].to_montgomery();
  }
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff = end - start;
  std::cout << "To Montomery time: " << diff.count() / r * 1e6 << " us\n";

  start = std::chrono::high_resolution_clock::now();
  for (uint64_t i = 0; i < r; ++i) {
    rand_b[i].from_montgomery();
  }
  end = std::chrono::high_resolution_clock::now();
  diff = end - start;
  std::cout << "From Montomery time: " << diff.count() / r * 1e6 << " us\n";
}

TEST(TestCrypto, BigInt_MontMult_repeat_perf) {
  uint64_t r = 1000000;
  BigInt rand_a = BigInt::random();
  std::vector<BigInt> rand_b(r);
  for (uint64_t i = 0; i < r; ++i) {
    rand_b[i] = BigInt::random().to_montgomery();
  }
  auto start = std::chrono::high_resolution_clock::now();
  for (uint64_t i = 0; i < r; ++i) {
    rand_a *= rand_b[i];
  }
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff = end - start;
  std::cout << "BigInt multiply time: " << diff.count() / r * 1e6 << " us\n";
}

TEST(TestCrypto, BigIntInv_perf) {
  uint64_t r = 100000;
  std::vector<BigInt> rand_a(r);
  for (uint64_t i = 0; i < r; ++i) {
    rand_a[i] = BigInt::random();
  }
  auto start = std::chrono::high_resolution_clock::now();
  for (uint64_t i = 0; i < r; ++i) {
    BigInt b = rand_a[i].inv();
  }
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff = end - start;
  std::cout << "BigInt inv time: " << diff.count() / r * 1e6 << " us\n";
}

TEST(TestCrypto, BigIntInvBatch_perf) {
  uint64_t r = 100000;
  std::vector<BigInt> rand_a(r);
  for (uint64_t i = 0; i < r; ++i) {
    rand_a[i] = BigInt::random();
  }
  std::vector<BigInt*> rand_a_ptrs(r);
  for (uint64_t i = 0; i < r; ++i) {
    rand_a_ptrs[i] = &rand_a[i];
  }
  auto start = std::chrono::high_resolution_clock::now();
  BigInt::inv_batch(rand_a_ptrs.data(), r);
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff = end - start;
  std::cout << "BigInt inv time: " << diff.count() / r * 1e6 << " us\n";
}

TEST(TestCrypto, hash_perf) {
  Label a = Label::random();
  int round = 10000000;
  auto start = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < round; ++i) {
    a = Hash(a);
  }
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff = end - start;
  std::cout << "Hash time: " << diff.count() / round * 1e6 << " us\n";
}

TEST(TestCrypto, ec_perf) {
  ECPoint A = ECPoint::generator_pow(BigInt::random());
  BigInt exponent = BigInt::random();
  int round = 100000;
  auto start = std::chrono::high_resolution_clock::now();
  for (int r = 0; r < round; ++r) {
    A.pow(exponent);
  }
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff = end - start;
  std::cout << "EC power time: " << diff.count() / round * 1e6 << " us\n";
}

TEST(TestCrypto, ec_perf_gen_pow) {
  uint64_t r = 1000000;
  std::vector<BigInt> rand_exponents(r);
  for (uint64_t i = 0; i < r; ++i) {
    rand_exponents[i] = BigInt::random();
  }
  auto start = std::chrono::high_resolution_clock::now();
  for (uint64_t i = 0; i < r; ++i) {
    ECPoint::generator_pow(rand_exponents[i]);
  }
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff = end - start;
  std::cout << "EC generator power time: " << diff.count() / r * 1e6 << " us\n";
}

TEST(TestCrypto, ec_perf_par) {
  ECPoint A = ECPoint::generator_pow(BigInt::random());
  BigInt exponent = BigInt::random();
#pragma omp parallel
  for (int r = 0; r < 100000; ++r) {
    A.pow(exponent);
  }
}

TEST(TestCrypto, ec_perf_par32) {
#pragma omp parallel for
  for (int i = 0; i < 32; ++i) {
    ECPoint A = ECPoint::generator_pow(BigInt::random());
    BigInt exponent = BigInt::random();
    for (int r = 0; r < 100000; ++r) {
      A.pow(exponent);
    }
  }
}

TEST(TestCrypto, ec_perf_gen_pow_par) {
  uint64_t r = 100000;
  std::vector<BigInt> rand_exponents(r);
  for (uint64_t i = 0; i < r; ++i) {
    rand_exponents[i] = BigInt::random();
  }
#pragma omp parallel for
  for (uint64_t i = 0; i < 100000; ++i) {
    ECPoint::generator_pow(rand_exponents[i]);
  }
}

TEST(TestCrypto, ec_perf_gen_pow_par_inner_loop) {
  uint64_t r = 100000;
  uint64_t inner_loop_size = 10;
  uint64_t outer_loop_size = r / inner_loop_size;
  std::vector<BigInt> rand_exponents(r);
  for (uint64_t i = 0; i < r; ++i) {
    rand_exponents[i] = BigInt::random();
  }

  for (uint64_t i = 0; i < outer_loop_size; ++i) {
#pragma omp parallel for num_threads(4)
    for (uint64_t j = 0; j < inner_loop_size; ++j) {
      ECPoint::generator_pow(rand_exponents[i * inner_loop_size + j]);
    }
  }
}
#endif