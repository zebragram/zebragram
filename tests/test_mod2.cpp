#include <gtest/gtest.h>

#include <vector>

#include "crypto/mod2.hpp"
#include "test_util.hpp"

using namespace ZebraGRAM;

TEST(Mod2CtxTest, BasicTest) {
  Mod2Ctx mod2;
  mod2.garble_init(default_sk, default_pk);
  std::cout << "Garble init done." << std::endl;
  mod2.eval_precompute(
      default_pk);  // for simplicity, use the same struct for eval
  std::cout << "Eval precompute done." << std::endl;
  ArithLabel x;
  secure_random_fmpz(x.auth, LEN_AUTH_SHARE);
  secure_random_fmpz(x.raw, LEN_RAW_SHARE);
  std::cout << "Random arith label sampled." << std::endl;
  Label y_g = mod2.garble_mod2(x, default_sk, default_pk);
  std::cout << "Garble mod2 done." << std::endl;
  Label y_e = mod2.eval_mod2(x, default_pk);
  std::cout << "Eval mod2 done." << std::endl;
  // check y_g == y_e
  EXPECT_EQ(y_g, y_e);
}

TEST(Mod2CtxTest, RandTest) {
  Mod2Ctx mod2;
  mod2.garble_init(default_sk, default_pk);
  std::cout << "Garble init done." << std::endl;
  mod2.eval_precompute(
      default_pk);  // for simplicity, use the same struct for eval
  std::cout << "Eval precompute done." << std::endl;
  ArithLabel x_g, x_e;
  for (int round = 0; round < 20; ++round) {
    int val = rand();
    fmpz_t val_fmpz;
    fmpz_init_set_ui(val_fmpz, val);
    secret_share(x_g, x_e, val_fmpz, default_sk);
    Label y_g = mod2.garble_mod2(x_g, default_sk, default_pk);
    Label y_e = mod2.eval_mod2(x_e, default_pk);
    // check y_g == y_e
    if (val % 2 == 0) {
      EXPECT_EQ(y_g, y_e);
    } else {
      EXPECT_EQ(y_g, y_e ^ key_manager.get_Delta());
    }
  }
}

TEST(Mod2CtxTest, RandBatchTest) {
  Mod2Ctx mod2;
  mod2.garble_init(default_sk, default_pk);
  mod2.eval_precompute(default_pk);

  const size_t batch = 32;
  std::vector<ArithLabel> x_gs(batch), x_es(batch);
  std::vector<int> vals(batch);

  for (size_t j = 0; j < batch; ++j) {
    int val = rand();
    vals[j] = val;
    fmpz_t val_fmpz;
    fmpz_init_set_ui(val_fmpz, val);
    secret_share(x_gs[j], x_es[j], val_fmpz, default_sk);
    fmpz_clear(val_fmpz);
  }

  std::vector<Label> y_gs(batch), y_es(batch);
  mod2.garble_mod2_batch(x_gs.data(), batch, y_gs.data(), default_sk,
                         default_pk);
  mod2.eval_mod2_batch(x_es.data(), batch, y_es.data(), default_pk);

  for (size_t j = 0; j < batch; ++j) {
    if (vals[j] % 2 == 0) {
      EXPECT_EQ(y_gs[j], y_es[j]);
    } else {
      EXPECT_EQ(y_gs[j], y_es[j] ^ key_manager.get_Delta());
    }
  }
}

TEST(Mod2CtxTest, Bench) {
  auto start = std::chrono::high_resolution_clock::now();
  paillier_keygen(default_pk, default_sk,
                  1024);  // 10 GB precompute tables for garbler
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff = end - start;
  std::cout << "Garbler Keygen time: " << diff.count() << " s" << std::endl;
  Mod2Ctx mod2;
  mod2.garble_init(default_sk, default_pk);
  std::cout << "Garble init done." << std::endl;
  start = std::chrono::high_resolution_clock::now();
  mod2.eval_precompute(default_pk,
                       6e3);  // for simplicity, use the same struct for eval
  end = std::chrono::high_resolution_clock::now();
  diff = end - start;
  std::cout << "Eval precompute time: " << diff.count() << " s" << std::endl;
  ArithLabel x;
  secure_random_fmpz(x.auth, LEN_AUTH_SHARE);
  secure_random_fmpz(x.raw, LEN_RAW_SHARE);
  uint64_t rounds = 1000;
  start = std::chrono::high_resolution_clock::now();
  for (uint64_t i = 0; i < rounds; ++i) {
    mod2.garble_mod2(x, default_sk, default_pk);
  }
  end = std::chrono::high_resolution_clock::now();
  diff = end - start;
  std::cout << "Garble mod2 time for " << rounds << " rounds: " << diff.count()
            << " s, avg: " << (diff.count() / rounds * 1e3) << " ms"
            << std::endl;
  start = std::chrono::high_resolution_clock::now();
  for (uint64_t i = 0; i < rounds; ++i) {
    mod2.eval_mod2(x, default_pk);
  }
  end = std::chrono::high_resolution_clock::now();
  diff = end - start;
  std::cout << "Eval mod2 time for " << rounds << " rounds: " << diff.count()
            << " s, avg: " << (diff.count() / rounds * 1e3) << " ms"
            << std::endl;

  // Batch benchmarks (batch size = 16)
  const size_t batch = 16;
  std::vector<ArithLabel> xs_g(batch), xs_e(batch);
  for (size_t j = 0; j < batch; ++j) {
    secure_random_fmpz(xs_g[j].auth, LEN_AUTH_SHARE);
    secure_random_fmpz(xs_g[j].raw, LEN_RAW_SHARE);
    secure_random_fmpz(xs_e[j].auth, LEN_AUTH_SHARE);
    secure_random_fmpz(xs_e[j].raw, LEN_RAW_SHARE);
  }
  std::vector<Label> ys_g(batch), ys_e(batch);

  start = std::chrono::high_resolution_clock::now();
  for (uint64_t i = 0; i < rounds; ++i) {
    mod2.garble_mod2_batch(xs_g.data(), batch, ys_g.data(), default_sk,
                           default_pk);
  }
  end = std::chrono::high_resolution_clock::now();
  diff = end - start;
  std::cout << "Garble mod2 batch size " << batch << " time for " << rounds
            << " batches: " << diff.count()
            << " s, avg per batch: " << (diff.count() / rounds * 1e3) << " ms"
            << ", avg per item: " << (diff.count() / (rounds * batch) * 1e3)
            << " ms" << std::endl;

  start = std::chrono::high_resolution_clock::now();
  for (uint64_t i = 0; i < rounds; ++i) {
    mod2.eval_mod2_batch(xs_e.data(), batch, ys_e.data(), default_pk);
  }
  end = std::chrono::high_resolution_clock::now();
  diff = end - start;
  std::cout << "Eval mod2 batch size " << batch << " time for " << rounds
            << " batches: " << diff.count()
            << " s, avg per batch: " << (diff.count() / rounds * 1e3) << " ms"
            << ", avg per item: " << (diff.count() / (rounds * batch) * 1e3)
            << " ms" << std::endl;
}

Mod2Ctx mod2;

struct Mod2TestMain : Gadget {
  const uint64_t payload_width_4k =
      (4096 * 8 + LEN_ARITH_DIGIT - 1) / LEN_ARITH_DIGIT;
  const uint64_t payload_width_1k =
      (1024 * 8 + LEN_ARITH_DIGIT - 1) / LEN_ARITH_DIGIT;
  // uint word_width = payload_width_1k;
  uint word_width = payload_width_4k;
  uint64_t decompose_bits = LEN_ARITH_DIGIT;
  std::function<FuncOutput(FuncInput)> main_func = [&](FuncInput) {
    std::vector<uint64_t> in_vals(word_width, 0);
    for (uint i = 0; i < word_width; ++i) {
      in_vals[i] = rand64();
    }
    ArithWord in_word(ArithWord::input_dbg(this, in_vals));
    if (get_mode() == GARBLE) {
      paillier_keygen(default_pk, default_sk,
                      10240);  // 10 GB precompute tables for garbler
      mod2.garble_init(default_sk, default_pk);
    } else if (get_mode() == EVAL) {
      mod2.eval_precompute(default_pk,
                           102400);  // 100 GB precompute tables for evaluator
    }
    auto start = std::chrono::high_resolution_clock::now();
    Word out_word =
        in_word.decompose_all(mod2, decompose_bits, default_sk, default_pk);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    std::cout << "Decompose_all time: " << diff.count() << " s" << std::endl;
    return FuncOutput{};
  };
  DEFINE_FUNC(main, {}, main_func);

  Mod2TestMain(uint64_t T, Mode mode) : Gadget(mode, T) {}
};

TEST(Mod2CtxTest, DecomposeAll) { test_gadget<Mod2TestMain>(ANY, {1, 1}); }