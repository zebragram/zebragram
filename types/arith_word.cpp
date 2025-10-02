#include "arith_word.hpp"

#include "fmpz_helpers.hpp"
#include "gadget.hpp"
#include "word.hpp"

namespace ZebraGRAM {
void ArithWord::join(const ArithWord &src, ArithWord &dst) {
  Assert_eq(src.get_owner(), dst.get_owner());

  Mode mode = src.get_mode();
  Gadget *owner = src.get_owner();
  if (mode == GARBLE) {
    Assert_eq(src.width(), dst.width());
    for (uint i = 0; i < src.width(); ++i) {
      ArithLabel diff = dst.payloads[i] - src.payloads[i];
      Assert_eq(GARBLE, dst.get_mode());
      Assert_eq(owner, dst.get_owner());
      owner->get_gc().write_arith_label(diff);
    }
  } else if (mode == EVAL) {
    dst.payloads.resize(src.width());
    for (uint i = 0; i < src.width(); ++i) {
      ArithLabel diff;
      owner->get_gc().read_arith_label(diff);
      dst.payloads[i] = src.payloads[i] + diff;
    }
  } else if (mode == MEASURE) {
    owner->get_gc().skip_arith_label(src.width());
  } else if (mode == DEBUG) {
    // do nothing
  }
}

ArithLabel bit_to_arith_label(const Bit &bit, const PaillierPrivKey &sk,
                              uint8_t *buffer) {
  ArithLabel res_label, gc_label;
  uint8_t lsb = bit.get_label().LSB();
  static constexpr uint len_auth_bytes = (LEN_K + LEN_STAT_PARAM + 63) / 64 * 8;
  // print input bit label
  // std::cout << "input bit label: " << bit.get_label() << " with lsb " <<
  // (int)lsb << std::endl;
  if (bit.get_mode() == GARBLE) {
    ArithLabel h0, h1;
    Label h0_raw = Hash(bit.get_label());
    Label h1_raw = Hash(bit.get_label() ^ key_manager.get_Delta());
    // convert to arith label

    uint8_t h0_auth_bytes[len_auth_bytes] = {0};
    uint8_t h1_auth_bytes[len_auth_bytes] = {0};
    prng(h0_raw, h0_auth_bytes,
         len_auth_bytes);  // expand h0_raw to len_auth_bytes
    prng(h1_raw, h1_auth_bytes, len_auth_bytes);

    fmpz_set_ui_array(h0.raw, (ulong *)h0_raw.get_ptr(), LAMBDA_BYTES / 8);
    fmpz_set_ui_array(h1.raw, (ulong *)h1_raw.get_ptr(), LAMBDA_BYTES / 8);
    fmpz_set_ui_array(h0.auth, (ulong *)h0_auth_bytes, len_auth_bytes / 8);
    fmpz_set_ui_array(h1.auth, (ulong *)h1_auth_bytes, len_auth_bytes / 8);

    if (lsb == 0) {
      fmpz_set(res_label.raw, h0.raw);
      fmpz_sub(gc_label.raw, h0.raw, h1.raw);
      fmpz_add_ui(gc_label.raw, gc_label.raw, 1ULL);
      // if gc_label.raw < 0, add 2^LAMBDA_BYTES
      if (fmpz_sgn(gc_label.raw) < 0) {
        fmpz_add_2exp(gc_label.raw, gc_label.raw, LAMBDA_BYTES * 8);
      }
      fmpz_set(res_label.auth, h0.auth);
      fmpz_sub(gc_label.auth, h0.auth, h1.auth);
      fmpz_add(gc_label.auth, gc_label.auth, sk.K);
      // if gc_label.auth < 0, add 2^LEN_AUTH_SHARE
      if (fmpz_sgn(gc_label.auth) < 0) {
        fmpz_add_2exp(gc_label.auth, gc_label.auth, len_auth_bytes * 8);
      }
    } else {
      fmpz_sub_ui(res_label.raw, h1.raw, 1ULL);
      fmpz_sub(gc_label.raw, h1.raw, h0.raw);
      fmpz_sub_ui(gc_label.raw, gc_label.raw, 1ULL);
      // if gc_label.raw < 0, add 2^LAMBDA_BYTES
      if (fmpz_sgn(gc_label.raw) < 0) {
        fmpz_add_2exp(gc_label.raw, gc_label.raw, LAMBDA_BYTES * 8);
      }
      fmpz_sub(res_label.auth, h1.auth, sk.K);
      fmpz_sub(gc_label.auth, h1.auth, h0.auth);
      fmpz_sub(gc_label.auth, gc_label.auth, sk.K);
      // if gc_label.auth < 0, add 2^LEN_AUTH_SHARE
      if (fmpz_sgn(gc_label.auth) < 0) {
        fmpz_add_2exp(gc_label.auth, gc_label.auth, len_auth_bytes * 8);
      }
    }
    // write gc to buffer
    gc_label.to_bytes(buffer);
  } else if (bit.get_mode() == EVAL) {
    // read gc from buffer
    Label h_raw = Hash(bit.get_label());
    uint8_t h_auth_bytes[len_auth_bytes] = {0};
    prng(h_raw, h_auth_bytes,
         len_auth_bytes);  // expand h_raw to len_auth_bytes
    fmpz_set_ui_array(res_label.raw, (ulong *)h_raw.get_ptr(),
                      LAMBDA_BYTES / 8);
    fmpz_set_ui_array(res_label.auth, (ulong *)h_auth_bytes,
                      len_auth_bytes / 8);
    if (lsb) {
      gc_label.from_bytes(buffer);
      res_label += gc_label;
      // check if res_label.raw > 2^(8 * LAMBDA_BYTES)
      if (fmpz_sizeinbase(res_label.raw, 2) > LAMBDA_BYTES * 8) {
        fmpz_sub_2exp(res_label.raw, res_label.raw, LAMBDA_BYTES * 8);
      }
      // check if res_label.auth > 2^(8 * len_auth_bytes)
      if (fmpz_sizeinbase(res_label.auth, 2) > len_auth_bytes * 8) {
        fmpz_sub_2exp(res_label.auth, res_label.auth, len_auth_bytes * 8);
      }
    }
  }
  // std::cout << "output arith label: " << res_label << std::endl;
  return res_label;
}

// gc_offset is the relative offset of the gc from the current gc pointer of the
// owner gadget
void bit_arith_word_mul(const Bit &bit, const ArithWord &word,
                        ArithWord &result, uint8_t *buffer) {
  Assert_eq(bit.get_owner(), word.get_owner());
  Assert_eq(bit.get_owner(), result.get_owner());
  Assert_eq(bit.get_mode(), word.get_mode());
  Assert_eq(bit.get_mode(), result.get_mode());
  Mode mode = bit.get_mode();
  ArithLabel bit_arith_label = bit_to_arith_label(bit, default_sk, buffer);
  // resize result payloads
  result.payloads.resize(word.width());
  if (mode == GARBLE) {
    // perform batch garble
    garble_mul_vec(result.payloads.data(), buffer + ArithLabel::byte_length,
                   bit_arith_label, word.payloads.data(), word.width(),
                   default_pk, default_sk);
  } else if (mode == EVAL) {
    eval_mul_vec(result.payloads.data(), buffer + ArithLabel::byte_length,
                 bit_arith_label, word.payloads.data(), word.width(),
                 default_pk);
  }
}

void bit_arith_word_mul(const Bit &bit, const ArithWord &word,
                        ArithWord &result) {
  const size_t buffer_size =
      ArithLabel::byte_length + GARBLE_MULT_SIZE * word.width();
  uint8_t *buffer = new uint8_t[buffer_size];
  Mode mode = bit.get_mode();
  if (mode == EVAL) {
    // read gc from the owner gadget
    bit.get_owner()->get_gc().read_data(buffer, buffer_size);
  }
  bit_arith_word_mul(bit, word, result, buffer);
  if (mode == GARBLE) {
    // write gc to the owner gadget
    bit.get_owner()->get_gc().write_data(buffer, buffer_size);
  }
  delete[] buffer;
}

void batch_garble_bit_word_mul(const std::vector<Bit> &bits,
                               const std::vector<ArithWord> &words,
                               std::vector<ArithWord> &results,
                               std::vector<uint8_t *> &buffers) {
  Gadget *owner = bits[0].get_owner();
  size_t word_width = words[0].width();
  std::vector<ArithLabel> bits_arith_labels(bits.size());
  // Prepare mega-batch for iter_1
  std::vector<const fmpz_t *> all_exps;
  std::vector<const PowmPrecomputeTable *> all_tables;
  size_t inner_iteration_1 = word_width + 1;
  size_t inner_iteration_2 = word_width;
  size_t iter_2_begin_index = bits.size() * inner_iteration_1 * 4;
  size_t total_tasks = iter_2_begin_index + bits.size() * inner_iteration_2 * 2;
  fmpz_t *all_powm_results = new fmpz_t[total_tasks];
  for (size_t i = 0; i < total_tasks; i++) {
    fmpz_init(all_powm_results[i]);
  }
  GarbleMulVecCtx *garble_ctx_arr = new GarbleMulVecCtx[bits.size()];
  all_exps.resize(bits.size() * inner_iteration_2 * 2 + iter_2_begin_index);
  all_tables.resize(bits.size() * inner_iteration_2 * 2 + iter_2_begin_index);

  size_t len_randomness = inner_iteration_1 * bits.size() * (LEN_K / 8);
  uint8_t *randomness = new uint8_t[len_randomness];
  RAND_bytes(randomness, len_randomness);
  size_t iter_1_total_ops = bits.size() * inner_iteration_1 * 4;

  for (size_t i = 0; i < bits.size(); ++i) {
    bits_arith_labels[i] = bit_to_arith_label(bits[i], default_sk, buffers[i]);
    results[i].payloads.resize(word_width);
    garble_ctx_arr[i].init_with_randomness(
        bits_arith_labels[i], words[i].payloads.data(), word_width, default_sk,
        randomness + i * inner_iteration_1 * (LEN_K / 8));
    for (size_t j = 0; j < inner_iteration_1; ++j) {
      size_t idx = (i * inner_iteration_1 + j) * 4;
      // garble_ctx_arr[i].push_iter_1_tasks(default_sk, all_exps, all_tables);
      all_exps[idx] = &garble_ctx_arr[i].r[j];
      all_exps[idx + 1] = &garble_ctx_arr[i].r[j];
      all_exps[idx + 2] = &garble_ctx_arr[i].r[j];
      all_exps[idx + 3] = &garble_ctx_arr[i].r[j];
      all_tables[idx] = &default_sk.table_g_p2;
      all_tables[idx + 1] = &default_sk.table_g_q2;
      all_tables[idx + 2] = &default_sk.table_pk_p2;
      all_tables[idx + 3] = &default_sk.table_pk_q2;
    }
    for (size_t j = 0; j < inner_iteration_2; ++j) {
      size_t idx = iter_2_begin_index + (i * inner_iteration_2 + j) * 2;
      garble_ctx_arr[i].set_iter_2_tasks(
          j, words[i].payloads.data(), default_sk, all_exps[idx],
          all_exps[idx + 1], all_tables[idx], all_tables[idx + 1]);
    }
  }
  ParTracker::start();
// Perform mega-batch
#pragma omp parallel num_threads(omp_get_max_threads())
  {
// powm_precomputed_batch(all_powm_results, all_exps.data(), all_tables.data(),
// total_tasks);
#pragma omp for
    for (size_t i = 0; i < total_tasks; i++) {
      powm_precomputed(all_powm_results[i], *all_exps[i], *all_tables[i]);
    }

#pragma omp for collapse(2)
    for (size_t i = 0; i < bits.size(); ++i) {
      for (size_t j = 0; j < inner_iteration_1 + inner_iteration_2; j++) {
        GarbleMulVecCtx &context = garble_ctx_arr[i];
        if (j < inner_iteration_1) {
          size_t res_start_index = (i * inner_iteration_1 + j) * 4;
          context.retrieve_iter_1_results(
              j, buffers[i] + ArithLabel::byte_length, default_pk, default_sk,
              all_powm_results[res_start_index],
              all_powm_results[res_start_index + 1],
              all_powm_results[res_start_index + 2],
              all_powm_results[res_start_index + 3]);
        } else {
          size_t res_start_index =
              iter_1_total_ops +
              (i * inner_iteration_2 + (j - inner_iteration_1)) * 2;
          context.retrieve_iter_2_results(
              j - inner_iteration_1, results[i].payloads.data(),
              words[i].payloads.data(), default_pk, default_sk,
              all_powm_results[res_start_index],
              all_powm_results[res_start_index + 1]);
        }
      }
    }
  }
  ParTracker::stop();

  // write all gc to the owner gadget
  for (size_t i = 0; i < bits.size(); ++i) {
    owner->get_gc().write_data(
        buffers[i],
        ArithLabel::byte_length + GARBLE_MULT_SIZE * words[i].width());
  }

  // clear
  for (size_t i = 0; i < total_tasks; i++) {
    fmpz_clear(all_powm_results[i]);
  }
  delete[] all_powm_results;
  delete[] garble_ctx_arr;
  delete[] randomness;
}

// Optimized batch evaluation that performs all operations in a single large
// batch
void batch_eval_bit_word_mul(const std::vector<Bit> &bits,
                             const std::vector<ArithWord> &words,
                             std::vector<ArithWord> &results,
                             const std::vector<uint8_t *> &buffers) {
  Assert_eq(bits.size(), words.size());
  Assert_eq(bits.size(), buffers.size());
  Assert_eq(bits.size(), results.size());

  size_t num_operations = bits.size();
  size_t word_width = words[0].width();
  for (const auto &word : words) {
    Assert_eq(word.width(), word_width);
    Assert(word.check());
  }

  // Calculate total modular exponentiations needed: 4 operations per bit-word
  // pair
  size_t total_tasks = num_operations * word_width * 4;

  // Allocate arrays for the mega-batch operation
  const fmpz_t **all_bases = new const fmpz_t *[total_tasks];
  const fmpz_t **all_exps = new const fmpz_t *[total_tasks];
  fmpz_t *all_powm_results = new fmpz_t[total_tasks];

  // Storage for bit arithmetic labels and deserialized ciphertexts
  std::vector<ArithLabel> bit_arith_labels(num_operations);

  // Allocate storage for ciphertext bases using raw arrays
  fmpz_t **ct_bases_1 = new fmpz_t *[num_operations];
  fmpz_t **ct_bases_2 = new fmpz_t *[num_operations];
  fmpz_t *base_3_array = new fmpz_t[num_operations];
  fmpz_t *base_4_array = new fmpz_t[num_operations];

  // Initialize all result arrays
  for (size_t i = 0; i < total_tasks; i++) {
    fmpz_init(all_powm_results[i]);
  }

  // Phase 1: Convert bits to arithmetic label and deserialize all ciphertexts
#pragma omp parallel for
  for (size_t op_idx = 0; op_idx < num_operations; op_idx++) {
    // Convert bit to arithmetic label
    bit_arith_labels[op_idx] =
        bit_to_arith_label(bits[op_idx], default_sk, buffers[op_idx]);

    // Initialize storage for this operation
    ct_bases_1[op_idx] = new fmpz_t[word_width];
    ct_bases_2[op_idx] = new fmpz_t[word_width];
    fmpz_init(base_3_array[op_idx]);
    fmpz_init(base_4_array[op_idx]);

    // Deserialize ciphertexts for this bit-word pair
    for (size_t w = 0; w < word_width; w++) {
      fmpz_init(ct_bases_1[op_idx][w]);
      fmpz_init(ct_bases_2[op_idx][w]);
      PaillierCiphertext::deserialize(
          ct_bases_1[op_idx][w], ct_bases_2[op_idx][w],
          buffers[op_idx] + ArithLabel::byte_length +
              (w + 1) * PAILLIER_CIPHER_TEXT_SIZE);
    }
    PaillierCiphertext::deserialize(base_3_array[op_idx], base_4_array[op_idx],
                                    buffers[op_idx] + ArithLabel::byte_length);
  }

  // Phase 2: Populate the mega-batch arrays
  size_t batch_idx = 0;
  for (size_t op_idx = 0; op_idx < num_operations; op_idx++) {
    for (size_t w = 0; w < word_width; w++) {
      // 4 operations per word element (same as in eval_mul_vec_batched)

      // Operation 1: ct_bases_1[op_idx][w] ^ bit_arith_labels[op_idx].auth
      all_bases[batch_idx] = &ct_bases_1[op_idx][w];
      all_exps[batch_idx] = &bit_arith_labels[op_idx].auth;
      batch_idx++;

      // Operation 2: ct_bases_2[op_idx][w] ^ bit_arith_labels[op_idx].raw
      all_bases[batch_idx] = &ct_bases_2[op_idx][w];
      all_exps[batch_idx] = &bit_arith_labels[op_idx].raw;
      batch_idx++;

      // Operation 3: base_3_array[op_idx] ^ words[op_idx].payloads[w].auth
      all_bases[batch_idx] = &base_3_array[op_idx];
      all_exps[batch_idx] = &words[op_idx].payloads[w].auth;
      batch_idx++;

      // Operation 4: base_4_array[op_idx] ^ words[op_idx].payloads[w].raw
      all_bases[batch_idx] = &base_4_array[op_idx];
      all_exps[batch_idx] = &words[op_idx].payloads[w].raw;
      batch_idx++;
    }
  }

  // Phase 3: Single mega-batch modular exponentiation
  ParTracker::start();
  powm_batch(all_powm_results, all_bases, all_exps, default_pk.N2, total_tasks);
  ParTracker::stop();

  // Phase 4: Post-process results
  for (size_t op_idx = 0; op_idx < num_operations; op_idx++) {
    results[op_idx].payloads.resize(word_width);
  }

#pragma omp parallel for
  for (size_t i = 0; i < num_operations * word_width; i++) {
    size_t op_idx = i / word_width;
    size_t w = i % word_width;
    size_t batch_idx = i * 4;

    // Get the 4 results for this word element
    fmpz_t combined_result, ddlog_result, y_combined, final_result;
    fmpz_init(combined_result);
    fmpz_init(ddlog_result);
    fmpz_init(y_combined);
    fmpz_init(final_result);

    // Combine the 4 results
    fmpz_mod_mul(combined_result, all_powm_results[batch_idx],
                 all_powm_results[batch_idx + 1], &default_pk.mod_N2_ctx);
    fmpz_mod_mul(combined_result, combined_result,
                 all_powm_results[batch_idx + 2], &default_pk.mod_N2_ctx);
    fmpz_mod_mul(combined_result, combined_result,
                 all_powm_results[batch_idx + 3], &default_pk.mod_N2_ctx);

    // Compute discrete log
    ddlog(ddlog_result, combined_result, default_pk);

    // Combine word payload and multiply by bit's raw value
    words[op_idx].payloads[w].merge(y_combined);
    fmpz_mul(y_combined, bit_arith_labels[op_idx].raw, y_combined);

    // final_result = y_combined - ddlog_result
    fmpz_sub(final_result, y_combined, ddlog_result);
    fmpz_mod_set_fmpz(final_result, final_result, &default_pk.mod_N_ctx);

    // Split result into output
    results[op_idx].payloads[w].split(final_result);

    // Clean up temporary variables
    fmpz_clear(combined_result);
    fmpz_clear(ddlog_result);
    fmpz_clear(y_combined);
    fmpz_clear(final_result);
  }

  // Cleanup
  for (size_t i = 0; i < total_tasks; i++) {
    fmpz_clear(all_powm_results[i]);
  }
  for (size_t op_idx = 0; op_idx < num_operations; op_idx++) {
    for (size_t w = 0; w < word_width; w++) {
      fmpz_clear(ct_bases_1[op_idx][w]);
      fmpz_clear(ct_bases_2[op_idx][w]);
    }
    delete[] ct_bases_1[op_idx];
    delete[] ct_bases_2[op_idx];
    fmpz_clear(base_3_array[op_idx]);
    fmpz_clear(base_4_array[op_idx]);
  }

  delete[] ct_bases_1;
  delete[] ct_bases_2;
  delete[] base_3_array;
  delete[] base_4_array;
  delete[] all_bases;
  delete[] all_exps;
  delete[] all_powm_results;
}

std::vector<ArithWord> batch_bit_arith_word_mul(
    const std::vector<Bit> &bits, const std::vector<ArithWord> &words) {
  Assert(bits.size() == words.size());
  Assert(!bits.empty());
  Mode mode = bits[0].get_mode();
  Gadget *owner = bits[0].get_owner();
  size_t word_width = words[0].width();
  for (size_t i = 0; i < bits.size(); ++i) {
    Assert_eq(bits[i].get_owner(), owner);
    Assert_eq(words[i].get_owner(), owner);
    Assert_eq(bits[i].get_mode(), mode);
    Assert_eq(words[i].get_mode(), mode);
    Assert_eq(words[i].width(),
              word_width);  // all words should have the same width
  }
  std::vector<ArithWord> results(bits.size(), ArithWord(bits[0].get_owner()));
  std::vector<uint8_t *> buffers(bits.size());
  for (size_t i = 0; i < bits.size(); ++i) {
    buffers[i] = new uint8_t[ArithLabel::byte_length +
                             GARBLE_MULT_SIZE * words[i].width()];
  }

  if (mode == EVAL) {
    // read all gc from the owner gadget
    for (size_t i = 0; i < bits.size(); ++i) {
      Assert_eq(bits[i].get_owner(), owner);
      Assert_eq(words[i].get_owner(), owner);
      Assert_eq(results[i].get_owner(), owner);
      Assert_eq(bits[i].get_mode(), mode);
      Assert_eq(words[i].get_mode(), mode);
      Assert_eq(results[i].get_mode(), mode);
      owner->get_gc().read_data(
          buffers[i],
          ArithLabel::byte_length + GARBLE_MULT_SIZE * words[i].width());
    }

    // Use optimized batched evaluation (no nested parallelism)
    batch_eval_bit_word_mul(bits, words, results, buffers);
  } else if (mode == GARBLE) {
    // Use optimized batched garbling (no nested parallelism)
    batch_garble_bit_word_mul(bits, words, results, buffers);
  }

  for (size_t i = 0; i < bits.size(); ++i) {
    delete[] buffers[i];
  }
  return results;
}

void batch_swap_arith_words(const std::vector<Bit> &bits,
                            std::vector<ArithWord> &words0,
                            std::vector<ArithWord> &words1) {
  Assert(bits.size() == words0.size());
  Assert(bits.size() == words1.size());
  Assert(!bits.empty());
  Mode mode = bits[0].get_mode();
  Gadget *owner = bits[0].get_owner();
  size_t word_width = words0[0].width();
  // compute words_tmp = words0 - words1
  std::vector<ArithWord> words_tmp(bits.size(), ArithWord(bits[0].get_owner()));
  for (size_t i = 0; i < bits.size(); ++i) {
    Assert_eq(bits[i].get_owner(), owner);
    Assert_eq(words0[i].get_owner(), owner);
    Assert_eq(words1[i].get_owner(), owner);
    Assert_eq(bits[i].get_mode(), mode);
    Assert_eq(words0[i].get_mode(), mode);
    Assert_eq(words1[i].get_mode(), mode);
    Assert_eq(words0[i].width(),
              word_width);  // all words should have the same width
    Assert_eq(words1[i].width(), word_width);
    words_tmp[i] = words0[i] - words1[i];
  }
  std::vector<ArithWord> words_delta =
      batch_bit_arith_word_mul(bits, words_tmp);
  for (size_t i = 0; i < bits.size(); ++i) {
    words0[i] -= words_delta[i];
    words1[i] += words_delta[i];
  }
}

// requires bits to be one-hot
void batch_swap_arith_words(const std::vector<Bit> &bits,
                            std::vector<ArithWord> &words0, ArithWord &word1) {
  Assert(bits.size() == words0.size());
  Assert(!bits.empty());
  Mode mode = bits[0].get_mode();
  Gadget *owner = bits[0].get_owner();
  size_t word_width = words0[0].width();
  // compute words_tmp = words0 - words1
  std::vector<ArithWord> words_tmp(bits.size(), ArithWord(bits[0].get_owner()));
  Assert_eq(word1.get_owner(), owner);
  Assert_eq(word1.get_mode(), mode);
  Assert_eq(word1.width(), word_width);
  for (size_t i = 0; i < bits.size(); ++i) {
    Assert_eq(bits[i].get_owner(), owner);
    Assert_eq(words0[i].get_owner(), owner);

    Assert_eq(bits[i].get_mode(), mode);
    Assert_eq(words0[i].get_mode(), mode);
    Assert_eq(words0[i].width(),
              word_width);  // all words should have the same width
    words_tmp[i] = words0[i] - word1;
  }
  // for (size_t i = 0; i < bits.size(); ++i)
  // {
  //     // auto tmp_revealed = words_tmp[i].reveal();
  //     // if (mode == EVAL) {
  //     //     std::cout << "words_tmp[" << i << "] revealed: " <<
  //     tmp_revealed[0] << std::endl;
  //     // }
  //     std::cout << "input " << (mode == GARBLE ? "GARBLE" : "EVAL") << " " <<
  //     words_tmp[i].payloads[0] << std::endl;
  // }
  // bits[0].get_owner()->dbg_check_gc_sync();
  std::vector<ArithWord> words_delta =
      batch_bit_arith_word_mul(bits, words_tmp);
  // print all words_delta
  // for (size_t i = 0; i < bits.size(); ++i)
  // {
  //     // auto delta_revealed = words_delta[i].reveal();
  //     // if (mode == EVAL) {
  //     //     std::cout << "words_delta[" << i << "] revealed: " <<
  //     delta_revealed[0] << std::endl;
  //     // }
  //     std::cout << "output " << (mode == GARBLE ? "GARBLE" : "EVAL") << " "
  //     << words_delta[i].payloads[0] << std::endl;
  // }
  for (size_t i = 0; i < bits.size(); ++i) {
    words0[i] -= words_delta[i];
    word1 += words_delta[i];
  }
}

ArithWord ArithWord::mux(const Bit &bit, const ArithWord &word1,
                         const ArithWord &word2) {
  Assert_eq(bit.get_owner(), word1.get_owner());
  Assert_eq(bit.get_owner(), word2.get_owner());
  Assert_eq(bit.get_mode(), word1.get_mode());
  Assert_eq(bit.get_mode(), word2.get_mode());
  ArithWord result(bit.get_owner());
  ArithWord word_diff = word2 - word1;
  bit_arith_word_mul(bit, word_diff, result);
  result += word1;
  return result;
}

ArithWord ArithWord::prng_from_label(Label seed, size_t width) {
  uint8_t *buffer = new uint8_t[ArithLabel::byte_length * width];
  prng(seed, buffer, ArithLabel::byte_length * width);
  ArithWord result;
  result.payloads.resize(width);
  for (size_t i = 0; i < width; ++i) {
    result.payloads[i].from_bytes(buffer + i * ArithLabel::byte_length);
  }
  delete[] buffer;
  return result;
}

ArithWord ArithWord::input_dbg(Gadget *owner, std::vector<uint64_t> values) {
  ArithWord result(owner);
  result.payloads.resize(values.size());
  for (uint i = 0; i < values.size(); ++i) {
    if (owner->get_mode() == GARBLE) {
      fmpz_set(result.payloads[i].raw, default_random_mask_raw);
      fmpz_set(result.payloads[i].auth, default_random_mask_auth);
    } else if (owner->get_mode() == EVAL) {
      fmpz_t val;
      fmpz_init_set_ui(val, values[i]);
      fmpz_add(result.payloads[i].raw, val, default_random_mask_raw);
      fmpz_mul(result.payloads[i].auth, val, default_sk.K);
      fmpz_add(result.payloads[i].auth, result.payloads[i].auth,
               default_random_mask_auth);
      fmpz_clear(val);
    }
  }
  return result;
}

ArithWord ArithWord::constant(Gadget *owner, size_t width, uint64_t value) {
  // garbler sets auth to default_random_mask_auth - value * sk.K
  // garbler sets raw to default_random_mask_raw - value
  // evaluator sets raw to default_random_mask_raw
  // evaluator sets auth to default_random_mask_auth
  ArithWord result(owner);
  result.payloads.resize(width);
  for (size_t i = 0; i < width; ++i) {
    if (owner->get_mode() == GARBLE) {
      if (value == 0) {
        fmpz_set(result.payloads[i].raw, default_random_mask_raw);
        fmpz_set(result.payloads[i].auth, default_random_mask_auth);
      } else {
        fmpz_t temp;
        fmpz_init(temp);
        fmpz_sub_ui(result.payloads[i].raw, default_random_mask_raw, value);
        fmpz_mul_ui(temp, default_sk.K, value);
        fmpz_sub(result.payloads[i].auth, default_random_mask_auth, temp);
        fmpz_clear(temp);
      }
    } else if (owner->get_mode() == EVAL) {
      fmpz_set(result.payloads[i].raw, default_random_mask_raw);
      fmpz_set(result.payloads[i].auth, default_random_mask_auth);
    } else {
      Assert(false);
    }
  }
  return result;
}

std::vector<ArithLabel> ArithWord::reveal() const {
  if (width() == 0) return {};
  if (get_mode() == GARBLE) {
    // write payloads to gc
    for (const ArithLabel &label : payloads) {
      get_owner()->get_gc().write_arith_label(label);
    }
  } else if (get_mode() == EVAL) {
    std::vector<ArithLabel> result(payloads.size());
    for (size_t i = 0; i < payloads.size(); ++i) {
      ArithLabel garbler_label;
      get_owner()->get_gc().read_arith_label(garbler_label);
      result[i] = payloads[i] - garbler_label;
    }
    return result;
  }
  return {};
}

Word ArithWord::decompose_all(const Mod2Ctx &mod2, uint num_bits,
                              const PaillierPrivKey &sk,
                              const PaillierPubKey &pk) const {
  Word res(get_owner(), width() * num_bits);
  std::vector<ArithLabel> cache(width());
  for (uint j = 0; j < width(); ++j) {
    cache[j] = payloads[j];
  }
  uint64_t buffer_len = ArithLabel::byte_length * num_bits * width();
  uint8_t *buffer = new uint8_t[buffer_len];
  if (get_mode() == EVAL) {
    // read all gc from the owner gadget
    get_owner()->get_gc().read_data(buffer, buffer_len);
  }
  for (uint i = 0; i < num_bits; ++i) {
    if (get_mode() == GARBLE) {
      // Batch compute all bits for this plane
      std::vector<Label> bool_labels(width());
      mod2.garble_mod2_batch(cache.data(), width(), bool_labels.data(), sk, pk);

      for (uint j = 0; j < width(); ++j) {
        Label bool_label = bool_labels[j];
        res[j * num_bits + i] = Bit(get_owner(), bool_label, 0, false, false);
        if (i == num_bits - 1) {
          continue;
        }
        ArithLabel arith_bool_label =
            bit_to_arith_label(res[j * num_bits + i], sk,
                               buffer + i * ArithLabel::byte_length * width() +
                                   j * ArithLabel::byte_length);
        cache[j] -= arith_bool_label;
        // divide cache[j] by 2
        fmpz_fdiv_q_2exp(cache[j].raw, cache[j].raw, 1);
        fmpz_fdiv_q_2exp(cache[j].auth, cache[j].auth, 1);
        // set cache[j].raw to at most LEN_RAW_SHARE - i
        fmpz_fdiv_r_2exp(cache[j].raw, cache[j].raw, LEN_RAW_SHARE - i);
        // set cache[j].auth to at most LEN_AUTH_SHARE - i
        fmpz_fdiv_r_2exp(cache[j].auth, cache[j].auth, LEN_AUTH_SHARE - i);
      }
    } else if (get_mode() == EVAL) {
      // Batch compute all bits for this plane
      std::vector<Label> bool_labels(width());
      mod2.eval_mod2_batch(cache.data(), width(), bool_labels.data(), pk);

      for (uint j = 0; j < width(); ++j) {
        Label bool_label = bool_labels[j];
        res[j * num_bits + i] = Bit(get_owner(), bool_label, 0, false, false);
        if (i == num_bits - 1) {
          continue;
        }
        ArithLabel arith_bool_label =
            bit_to_arith_label(res[j * num_bits + i], sk,
                               buffer + i * ArithLabel::byte_length * width() +
                                   j * ArithLabel::byte_length);
        cache[j] -= arith_bool_label;
        // divide cache[j] by 2
        fmpz_fdiv_q_2exp(cache[j].raw, cache[j].raw, 1);
        fmpz_fdiv_q_2exp(cache[j].auth, cache[j].auth, 1);
        // set cache[j].raw to at most LEN_RAW_SHARE - i
        fmpz_fdiv_r_2exp(cache[j].raw, cache[j].raw, LEN_RAW_SHARE - i);
        // set cache[j].auth to at most LEN_AUTH_SHARE - i
        fmpz_fdiv_r_2exp(cache[j].auth, cache[j].auth, LEN_AUTH_SHARE - i);
      }
    }
  }
  if (get_mode() == GARBLE) {
    get_owner()->get_gc().write_data(buffer, buffer_len);
  }
  return res;
}
}  // namespace ZebraGRAM