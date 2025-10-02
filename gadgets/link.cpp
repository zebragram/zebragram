#include "link.hpp"

namespace ZebraGRAM {
SIMDLink::SIMDLink(Gadget* callee) : BaseConditionalLink(callee, SIMD_COND) {}

bool SIMDLink::update_link(const Bit& control, const BigInt& caller_label) {
  Assert(control.get_owner() == caller);
  if (last_caller_time == caller->get_time()) {
    return is_active;
  }
  last_caller_time = caller->get_time();
  inc_time();
#ifdef MEASURE_STACK_COST
  measure_stack_flag = true;
#endif

  Bit control_e = control.reveal();
  Assert_eq(control.get_mode(), mode);
  if (mode == GARBLE) {
    // only record the label
    Assert(caller_label.is_set());
    simd_labels.push_back(caller_label);
    controls.push_back(control_e);
    prefix_sum.emplace_back(prefix_sum.back() + control_e);
  } else if (mode == EVAL) {
    // update the relative label
    is_active = !control_e.get_value();
    if (is_active) {
      callee->inc_time();
      uint64_t tau = get_time();
      relative_label = HashZ(control_e.get_label(), tau);
      inv_relative_label.unset();
      // route through the compaction network
      for (uint i = 0;; ++i) {
        const Word& curr_sum = prefix_sum[tau];
        if (i < curr_sum.width()) {
          const Bit& curr_bit = curr_sum[i];
          Assert(curr_bit.is_pub_e());
          uint64_t index = index_calculator(i, tau);
          if (curr_bit.get_value()) {
            // SIMD Buffer
            // TODO check nonce uniqueness
            // relative_label *= HashZ((!curr_bit).get_label(), 2 * index);

            BigInt join_material;
            uint8_t join_material_bytes[BigInt::byte_length];
            gc.peak_data(join_material_bytes, BigInt::byte_length * index,
                         BigInt::byte_length);
            join_material.dec((!curr_bit).get_label(), join_material_bytes,
                              2 * index);
            // SIMD JOIN
            relative_label *= join_material;
            tau -= (1UL << i);
          } else {
            // only SIMD Buffer
            relative_label *= HashZ(curr_bit.get_label(), 2 * index + 1);
          }
        } else {
          break;
        }
      }
    }
    prefix_sum.emplace_back(prefix_sum.back() + control_e);
  } else if (mode == DEBUG) {
    is_active = !control_e.get_value();
    if (is_active) {
      callee->inc_time();
    }
    prefix_sum.emplace_back(prefix_sum.back() + control_e);
  } else if (mode == MEASURE) {
    prefix_sum[0] = prefix_sum[0] + control_e;
  }
#ifdef MEASURE_STACK_COST
  measure_stack_flag = false;
#endif

  return !control_e.get_value();
}

std::vector<SIMDWord> SIMDLink::translate(const std::vector<SIMDWord>& inputs) {
  for (const SIMDWord& input : inputs) {
    Assert(input.get_mode() == EVAL || input.get_mode() == DEBUG);
    Assert(input.get_owner() == caller || input.get_owner() == callee);
  }
  std::vector<SIMDWord> outputs = inputs;
  if (!inputs.empty()) {
    bool is_caller_to_callee = inputs[0].get_owner() == caller;
    for (SIMDWord& output : outputs) {
      Assert_eq(output.get_owner() == caller, is_caller_to_callee);
      output.set_owner(is_caller_to_callee ? callee : caller);
      if (relative_label.is_set()) {
        if (is_caller_to_callee) {
          for (SIMDWord::SIMDEncoding& encoding : output.encodings) {
            encoding.pow(relative_label);
          }
        } else {
          if (!inv_relative_label.is_set()) {
            inv_relative_label = relative_label.inv().to_montgomery();
          }
          for (SIMDWord::SIMDEncoding& encoding : output.encodings) {
            encoding.pow(inv_relative_label);
          }
          // the garbler cannot infer the pub_es of the callee's output,
          // assuming all return values are unknown to the evaluator
          std::fill(output.pub_es.begin(), output.pub_es.end(), false);
        }
      }
    }
  }
  return outputs;
}

GCPtr SIMDLink::garble(const GCPtr& gc_begin) {
  // std::cout << "garbling SIMDLink for " << T << " timesteps" << std::endl;
  const uint64_t inv_batch_size = std::min(32ul, T);
  std::vector<BigInt*> big_int_to_inv_ptrs(inv_batch_size);
  BaseGadget::garble(gc_begin);
  uint depth = log2ceil(T);
  Assert_eq(get_time() + 1, T);
  if (mode == GARBLE) {
    Assert_eq(controls.size(), T);
    Assert_eq(simd_labels.size(), T);
  }
  init_gc_ptr(gc_begin);
  if (mode == MEASURE) {
#ifdef MEASURE_STACK_COST
    measure_stack_flag = true;
#endif
    gc.skip_big_int(T * depth - (1UL << depth) + 1);
#ifdef MEASURE_STACK_COST
    measure_stack_flag = false;
#endif
    return gc;
  }
  // the language of wires joined with another wire whose language
  // hasn't been computed
  std::vector<std::queue<BigInt>> pending_joins;
  pending_joins.resize(depth);
  uint64_t index = -1;
  uint64_t index2 = -1;
  std::vector<std::vector<BigInt>> hashz_stack(inv_batch_size,
                                               std::vector<BigInt>(depth));
  uint64_t inv_batch_offset = 0;

  uint64_t batch_begin_t = 0;
  for (uint64_t t = 0; t < T; ++t) {
    simd_labels[t] *= HashZ(controls[t].get_label(), t);
    for (uint i = 0; i < depth; ++i) {
      if (i < prefix_sum[t].width()) {
        const BigInt& hashz =
            HashZ(prefix_sum[t][i].get_label(), 2 * index + 1);

        // Notice that the simd labels we use for computing the join
        // material forms a matrixes like this

        // h00            h01            h02          h03
        // h00*h10        h01*h11        h02*h12      h03*h13
        // h00*h10*h20    h01*h11*h21    h02*h12*h22  h03*h13*h23

        // The join materials can be computed once we know the inverse of
        // every entry of the matrix.
        // To compute the inverse of every entry in this matrix, we only
        // need to compute the inverse of the last row. The rest can be computed
        // with multiplication on the last row.

        // we could further optimize this by doing batch inversion
        // (by computing the prefix and suffix products of the last row)
        // however, this requires more memory to store the stacks and makes the
        // code more complex
        if (i == prefix_sum[t].width() - 1) {
          hashz_stack[inv_batch_offset][i] = simd_labels[t];
        } else {
          hashz_stack[inv_batch_offset][i] = hashz;
        }
        simd_labels[t] *= hashz;
        ++index;
      }
      if (t + (1UL << i) < T) {
        pending_joins[i].push(simd_labels[t]);
      }
    }
    ++inv_batch_offset;
    uint64_t actual_inv_batch_size =
        std::min(inv_batch_size, T - batch_begin_t);
    if (inv_batch_offset == actual_inv_batch_size) {
      // todo perform the inversions in batch
      for (uint64_t offset = 0; offset < actual_inv_batch_size; ++offset) {
        uint64_t tt = batch_begin_t + offset;
        BigInt& big_int_to_inv =
            hashz_stack[offset][prefix_sum[tt].width() - 1];
        big_int_to_inv_ptrs[offset] = &big_int_to_inv;
      }
      BigInt::inv_batch(big_int_to_inv_ptrs.data(), actual_inv_batch_size);
      for (uint64_t offset = 0; offset < actual_inv_batch_size; ++offset) {
        uint64_t tt = batch_begin_t + offset;
        Assert_less(prefix_sum[tt].width() - 1, depth);
        for (int i = prefix_sum[tt].width() - 2; i >= 0; --i) {
          hashz_stack[offset][i] *= hashz_stack[offset][i + 1];
        }
        for (uint i = 0; i < depth; ++i) {
          if (i < prefix_sum[tt].width()) {
            if (tt >= (1UL << i)) {
              // we can compute a pending join
              Bit neg_s = !prefix_sum[tt][i];
              // const BigInt& h_neg = HashZ(neg_s.get_label(), 2 * index);
              Assert(!pending_joins[i].empty());
              BigInt& pending_join = pending_joins[i].front();
              pending_join *= hashz_stack[offset][i];
              // }

              uint8_t encrypted_join_material[BigInt::byte_length];
              pending_join.enc(neg_s.get_label(), encrypted_join_material,
                               2 * index2);
              pending_joins[i].pop();
              // whenever the join material is computed, write to the GC
              // sequentially
              gc.write_data(encrypted_join_material, BigInt::byte_length);
            }
            ++index2;
          }
        }
      }
      Assert_eq(index, index2);
      inv_batch_offset = 0;
      batch_begin_t = t + 1;
    }
  }

  // free the memory
  clear_and_release(pending_joins);
  clear_and_release(prefix_sum);
  clear_and_release(controls);
  return gc;
}

Link::Link(Gadget* callee) : BaseConditionalLink(callee, COND) {}

bool Link::update_link(const Bit& control) {
  Assert(control.get_owner() == caller);
  if (last_caller_time == caller->get_time()) {
    return is_active;
  }
  last_caller_time = caller->get_time();
  inc_time();
#ifdef MEASURE_STACK_COST
  measure_stack_flag = true;
#endif
  // get_caller()->dbg_check_gc_sync();
  Bit control_e = control.reveal();
  curr_bit_offset = 0;
  curr_arith_digit_offset = 0;
  Assert_eq(control.get_mode(), mode);
  if (mode == GARBLE) {
    // just record the input
    controls.push_back(control_e);
    port_words.emplace_back();
    prefix_sum.emplace_back(prefix_sum.back() + control_e);
  } else if (mode == EVAL || mode == DEBUG) {
    controls.resize(1);
    controls[0] = control_e;
    is_active = !control_e.get_value();
    if (is_active) {
      callee->inc_time();
      // uint64_t tau = get_time();
    }
    prefix_sum.emplace_back(prefix_sum.back() + control_e);
  } else if (mode == MEASURE) {
    if (port_words.empty()) {
      port_words.emplace_back();
    }
    prefix_sum[0] = prefix_sum[0] + control_e;
  }
#ifdef MEASURE_STACK_COST
  measure_stack_flag = false;
#endif
  return !control_e.get_value();
}

void Link::add_caller_words(FuncInput parent_words) {
  Assert(!port_words.empty());
  FuncOutput& latest_words = port_words.back();
  if (get_mode() == GARBLE) {
    latest_words.insert(latest_words.end(), parent_words.begin(),
                        parent_words.end());
  } else {
    Assert(get_mode() == MEASURE);
    if (last_caller_time == 0) {
      // only insert once
      latest_words.insert(latest_words.end(), parent_words.begin(),
                          parent_words.end());
    }
  }
}

FuncOutput Link::retrieve_callee_words(uint64_t num) {
  Assert(garbled);
  Assert(get_mode() == GARBLE || get_mode() == MEASURE);

  auto word_begin =
      port_words[curr_retrieved_t].begin() + curr_retrieved_word_idx;
  Assert_less(curr_retrieved_t, port_words.size());
  Assert(curr_retrieved_word_idx + num <= port_words[curr_retrieved_t].size());
  FuncInput output = std::vector<Word>(word_begin, word_begin + num);
  curr_retrieved_word_idx += num;
  if (curr_retrieved_word_idx == port_words[curr_retrieved_t].size()) {
    if (get_mode() == GARBLE) {
      ++curr_retrieved_t;
    }
    curr_retrieved_word_idx = 0;
  }
  return output;
}

void Link::clear_callee_words() { clear_and_release(port_words); }

FuncOutput Link::translate(FuncInput input) {
  Assert(is_active);
  if (input.empty()) {
    return {};
  }
  Assert(get_mode() == EVAL || get_mode() == DEBUG);
  FuncOutput output = input;
  Gadget* input_owner = input[0].get_owner();
  Assert(input_owner == caller || input_owner == callee);
  bool is_caller_to_callee = input_owner == caller;
  Gadget* output_owner = is_caller_to_callee ? callee : caller;

  for (Word& word : output) {
    // the garbler cannot predict the routing behavior
    word.set_pub_g(false);
    word.set_owner(output_owner);
    if (!is_caller_to_callee) {
      word.set_pub_e(false);
    }
  }
  if (mode == EVAL) {
    if (word_width_sum == (uint32_t)-1) {
      gc.read_default(word_width_sum);
      gc.read_default(arith_word_width_sum);
    }
    uint32_t combined_join_material_len =
        word_width_sum * sizeof(Label) +
        arith_word_width_sum * ArithLabel::byte_length;
    uint depth = log2ceil(T);
    // update the relative label
    if (is_active) {
      const Bit& control_e = controls[0];
      uint64_t tau = get_time();
      uint64_t bit_offset = curr_bit_offset;
      for (Word& word : output) {
        Label h;
        for (uint bit_idx = 0; bit_idx < word.width(); ++bit_idx) {
          h = Hash(control_e.get_label(), t, bit_offset++);
          word[bit_idx].get_label() ^= h;
        }
        if (word.has_payload()) {
          const ArithWord& diff_h =
              ArithWord::prng_from_label(h, word.get_payload().width());
          if (is_caller_to_callee) {
            word.get_payload() += diff_h;
          } else {
            word.get_payload() -= diff_h;
          }
        }
      }

      // route through the compaction network

      for (uint i = 0; i < depth; ++i) {
        const Word& curr_sum = prefix_sum[tau];
        if (i >= curr_sum.width()) {
          break;
        }
        const Bit& curr_bit = curr_sum[i];
        Assert(curr_bit.is_pub_e());
        uint64_t index = index_calculator(i, tau);
        bit_offset = curr_bit_offset;
        if (curr_bit.get_value()) {
          std::vector<Label> join_material(output.total_bit_width());
          std::vector<ArithLabel> arith_join_material(
              output.total_arith_width());
          uint64_t peak_offset =
              index * combined_join_material_len + curr_gc_offset;
          (gc + peak_offset).peak_default(join_material, 0);
          peak_offset += join_material.size() * sizeof(Label);
          if (arith_word_width_sum > 0) {
            // peak arith labels one by one
            for (uint j = 0; j < arith_join_material.size(); ++j) {
              (gc + peak_offset).peak_arith_label(arith_join_material[j], j);
              // std::cout << "Read arith label at level "<< i << " and gc
              // offset " << (gc + peak_offset).get_offset() + j *
              // ArithLabel::byte_length << " for join material: " <<
              // arith_join_material[j] << std::endl;
            }
          }
          uint64_t arith_bit_inner_offset = 0;
          for (Word& word : output) {
            Label h_neg;
            for (uint bit_idx = 0; bit_idx < word.width();
                 ++bit_idx, ++bit_offset) {
              h_neg = Hash((!curr_bit).get_label(), 2 * index, bit_offset);
              word[bit_idx].get_label() ^=
                  join_material[bit_offset - curr_bit_offset] ^ h_neg;
            }
            if (word.has_payload()) {
              const ArithWord& h_neg_arith_word =
                  ArithWord::prng_from_label(h_neg, word.get_payload().width());
              ArithWord join_arith_word(output_owner);
              join_arith_word.set_payload(
                  arith_join_material.begin() + arith_bit_inner_offset,
                  arith_join_material.begin() + arith_bit_inner_offset +
                      word.get_payload().width());
              if (is_caller_to_callee) {
                word.get_payload() += (join_arith_word + h_neg_arith_word);
              } else {
                word.get_payload() -= (join_arith_word + h_neg_arith_word);
              }
              arith_bit_inner_offset += word.get_payload().width();
            }
          }
          tau -= (1UL << i);
        } else if (tau > 0) {
          // only Buffer
          for (Word& word : output) {
            Label h;
            for (uint bit_idx = 0; bit_idx < word.width(); ++bit_idx) {
              h = Hash(curr_bit.get_label(), 2 * index + 1, bit_offset++);
              word[bit_idx].get_label() ^= h;
            }
            if (word.has_payload()) {
              const ArithWord& diff_h =
                  ArithWord::prng_from_label(h, word.get_payload().width());
              if (is_caller_to_callee) {
                word.get_payload() += diff_h;
              } else {
                word.get_payload() -= diff_h;
              }
            }
          }
        }
      }
      curr_bit_offset =
          (curr_bit_offset + output.total_bit_width()) % word_width_sum;
      curr_gc_offset = (curr_gc_offset + output.total_gc_size()) %
                       combined_join_material_len;
      if (arith_word_width_sum > 0) {
        curr_arith_digit_offset =
            (curr_arith_digit_offset + output.total_arith_width()) %
            arith_word_width_sum;
      }
    }
  }
  return output;
}
#ifdef MEASURE_TSC_STACK
uint64_t tsc_stack_cost_helper(int64_t m, uint w, uint64_t T) {
  if (m <= 0) {
    Assert(T <= 0);
    return 0;
  }
  return (m * 3 + m / 2 * 4 + m * w * 4) * LAMBDA_BYTES +
         tsc_stack_cost_helper(m / 2, w * 2, (T - 3) / 2);
}

uint64_t tsc_stack_cost(uint64_t T, uint w) {
  return tsc_stack_cost_helper(T, w, T / 2);
}
#endif

GCPtr Link::garble(const GCPtr& gc_begin) {
  BaseGadget::garble(gc_begin);
  Assert(T);
  uint depth = log2ceil(T);
  Assert_eq(get_time() + 1, T);
  if (mode == GARBLE) {
    Assert_eq(controls.size(), T);
    Assert_eq(port_words.size(), T);
  }
  init_gc_ptr(gc_begin);
  word_width_sum = port_words[0].total_bit_width();
  arith_word_width_sum = port_words[0].total_arith_width();
  if (mode == MEASURE) {
#ifdef MEASURE_STACK_COST
    measure_stack_flag = true;
#endif
    gc.skip_default<uint32_t>(2);
#ifdef MEASURE_TSC_STACK
    global_tsc_stack_cost += tsc_stack_cost(T, word_width_sum);
#endif
    gc.skip_label((T * depth - (1UL << depth) + 1) * word_width_sum);
    gc.skip_arith_label((T * depth - (1UL << depth) + 1) *
                        arith_word_width_sum);

#ifdef MEASURE_STACK_COST
    measure_stack_flag = false;
#endif
  } else {
    // the evaluator needs to learn this to peak data
    gc.write_default(word_width_sum);
    gc.write_default(arith_word_width_sum);
    // the language of wires joined with another wire whose language
    // hasn't been computed
    std::vector<std::queue<FuncOutput>> pending_joins;
    pending_joins.resize(depth);
    // std::cout << "Garbling Link with T = " << T << ", depth = " << depth <<
    // std::endl;
#ifdef MAX_GADGET_TIME
    GCPtr test_end_gc_ptr = gc;
    bool overwrote = false;
#ifdef TOTAL_TIME
    // use chernoff bound to estimate the maximum number of times a gadget is
    // called in TOTAL_TIME timesteps with high probability const double mu =
    // (double)T * MAX_GADGET_TIME / TOTAL_TIME; const double delta = sqrt(30.0
    // / mu); const uint64_t chernoff_bound = ceil(mu * (1 + delta)); const
    // uint64_t threshold_time = std::min((uint64_t)MAX_GADGET_TIME,
    // chernoff_bound);
    const uint64_t threshold_time = chernoff_upper_bound(
        TOTAL_TIME, MAX_GADGET_TIME, (double)TOTAL_TIME / T, -40);
#else
    const uint64_t threshold_time = MAX_GADGET_TIME;
#endif
#endif
    for (uint64_t t = 0; t < T; ++t) {
      uint64_t bit_offset = 0;
      for (Word& word : port_words[t]) {
        Label h;
        for (uint bit_idx = 0; bit_idx < word.width(); ++bit_idx) {
          h = Hash(controls[t].get_label(), t, bit_offset++);
          word[bit_idx].get_label() ^= h;
        }
        if (word.has_payload()) {
          const ArithWord& diff_h =
              ArithWord::prng_from_label(h, word.get_payload().width());
          word.get_payload() += diff_h;
        }
      }
      Assert_eq(bit_offset, word_width_sum);
      for (uint i = 0; i < depth; ++i) {
        if (i < prefix_sum[t].width()) {
          uint64_t index = index_calculator(i, t);
          if (t >= (1UL << i)) {
            // we can compute a pending join
            Bit neg_s = !prefix_sum[t][i];
            Assert(!pending_joins[i].empty());
            FuncInput pending_join = pending_joins[i].front();
            Assert_eq(pending_join.size(), port_words[t].size());
            bit_offset = 0;
            for (uint64_t word_idx = 0; word_idx < pending_join.size();
                 ++word_idx) {
              const Word& pending_join_word = pending_join[word_idx];
              const Word& port_word = port_words[t][word_idx];
              Assert_eq(pending_join_word.width(), port_word.width());
              Label h_neg;
              for (uint bit_idx = 0; bit_idx < pending_join_word.width();
                   ++bit_idx) {
                h_neg = Hash(neg_s.get_label(), 2 * index, bit_offset);

                Label join_material = pending_join_word[bit_idx].get_label() ^
                                      port_word[bit_idx].get_label() ^ h_neg;
                gc.write_label(join_material);
                ++bit_offset;
              }
              if (port_word.has_payload()) {
                Assert(pending_join_word.has_payload());
                ArithWord h_neg_arith_word = ArithWord::prng_from_label(
                    h_neg, port_word.get_payload().width());
                // TODO: check if should reverse
                const ArithWord& arith_join_material =
                    pending_join_word.get_payload() - port_word.get_payload() -
                    h_neg_arith_word;
                for (uint k = 0; k < arith_join_material.width(); ++k) {
                  // std::cout << "Writing arith label at level " << i << " and
                  // gc offset " << gc.get_offset() << "for join material: " <<
                  // arith_join_material.payloads[k] << std::endl;
                  gc.write_arith_label(arith_join_material.payloads[k]);
                }
              }
            }
            bit_offset = 0;
            for (uint64_t word_idx = 0; word_idx < pending_join.size();
                 ++word_idx) {
              Label h =
                  Hash(prefix_sum[t][i].get_label(), 2 * index + 1, bit_offset);
              Word& port_word = port_words[t][word_idx];
              const Word& pending_join_word = pending_join[word_idx];
              Assert_eq(pending_join_word.width(), port_word.width());
              for (uint bit_idx = 0; bit_idx < pending_join_word.width();
                   ++bit_idx) {
                h = Hash(prefix_sum[t][i].get_label(), 2 * index + 1,
                         bit_offset);
                port_word[bit_idx].get_label() ^= h;
                ++bit_offset;
              }
              if (port_word.has_payload()) {
                Assert_eq(pending_join_word.get_payload().width(),
                          port_word.get_payload().width());
                const ArithWord& diff_h = ArithWord::prng_from_label(
                    h, port_word.get_payload().width());
                port_word.get_payload() += diff_h;
              }
            }

            pending_joins[i].pop();
          }
        }
        if (t + (1UL << i) < T) {
          Assert_eq(port_words[t].size(), port_words[0].size());
          pending_joins[i].push(port_words[t]);
        }
      }
#ifdef MAX_GADGET_TIME
      if (t == threshold_time) {
        test_end_gc_ptr = get_gc();
        overwrote = true;
      }
#endif
#ifdef PRINT_PROGRESS_GRANULARITY
      uint64_t curr_gc_offset = get_gc().get_offset();
#ifdef MAX_GADGET_TIME
      curr_gc_offset += overwritten_bytes;
#endif
      print_curr_gc_offset(curr_gc_offset);
#endif
    }

#ifdef MAX_GADGET_TIME
    if (overwrote) {
      uint64_t local_overwritten_bytes =
          get_gc().get_offset() - test_end_gc_ptr.get_offset();
      overwritten_bytes += local_overwritten_bytes;
      this->gc.rewind_write(local_overwritten_bytes);
    }
#endif
    // std::cout << "Arith material gc offset ends at " << gc.get_offset() <<
    // std::endl; Assert(pending_joins.empty());
    clear_and_release(pending_joins);
    // free the memory
    clear_and_release(prefix_sum);
    clear_and_release(controls);
  }
  for (FuncOutput& words : port_words) {
    for (uint64_t word_idx = 0; word_idx < words.size(); ++word_idx) {
      // the garbler cannot predict the routing behavior
      words[word_idx].set_pub_g(false);
      words[word_idx].set_owner(callee);
    }
  }
  return gc;
}

uint64_t BaseConditionalLink::index_calculator(uint row, uint64_t col) {
  if (!row && !col) {
    return -1;
  }
  uint64_t log_col_group = log2ceil(col);
  uint64_t prev_pow2 = 1UL << (log_col_group - 1);
  uint64_t precede_col_group =
      prev_pow2 * (log_col_group - 2) + 1 + log_col_group * (col - prev_pow2);
  return precede_col_group + (uint64_t)row;
}
}  // namespace ZebraGRAM