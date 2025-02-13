#include "simd_word.hpp"

#include "gadget.hpp"
#include "key_manager.hpp"

namespace PicoGRAM {

constexpr uint SIMDWord::max_pack_width;

SIMDWord::SIMDWord(Gadget* owner, uint width, uint64_t bit_offset,
                   const BigInt& label, const std::vector<bool>& pub_es)
    : DataType(owner), label(label), bit_offset(bit_offset), values(width, 0) {
  if (pub_es.empty()) {
    this->pub_es = std::vector<bool>(width, false);
  } else {
    this->pub_es = pub_es;
  }
}

void SIMDWord::from_word(const Word& word, uint64_t bit_offset) {
  std::vector<SIMDWord> simd_words = from_words({word}, bit_offset);
  std::swap(*this, simd_words[0]);
}

Word SIMDWord::to_word() const { return to_words({*this})[0]; }

std::vector<SIMDWord> SIMDWord::from_words(const std::vector<Word>& words,
                                           uint64_t begin_bit_offset) {
  uint num_words = words.size();
  Assert(num_words > 0);
  std::vector<SIMDWord> simd_words;
  Mode mode = words[0].get_mode();
  Gadget* owner = words[0].get_owner();
  GCPtr& gc = owner->get_gc();
  uint64_t bit_offset = begin_bit_offset;
  simd_words.resize(num_words);
  for (uint word_idx = 0; word_idx < num_words; ++word_idx) {
    const Word& word = words[word_idx];
    Assert_eq(word.get_mode(), mode);
    Assert_eq(word.get_owner(), owner);
    SIMDWord& simd_word = simd_words[word_idx];

    uint width = word.width();
    simd_word.set_owner(owner);
    simd_word.values.resize(width);
    simd_word.pub_es.resize(width);
    for (uint i = 0; i < width; ++i) {
      simd_word.values[i] = word[i].get_value();
      simd_word.pub_es[i] = word[i].is_pub_e();
    }
    switch (mode) {
      case GARBLE:
        Assert(owner);
        simd_word.bit_offset = bit_offset;
        bit_offset += word.width();
        if (!owner->get_simd_label().is_set()) {
          owner->set_simd_label(BigInt::random());
        }
        simd_word.label = owner->get_simd_label();
        Assert(simd_word.label.is_set());

        break;
      case EVAL:
        if (!word.is_skip()) {
          simd_word.label = BigInt::unity();
          simd_word.encodings.resize((width + max_pack_width - 1) /
                                     max_pack_width);
          for (uint offset = 0, i = 0; offset < width;
               offset += max_pack_width, ++i) {
            uint pack_width = std::min(max_pack_width, width - offset);
            Label aggr_encoding = word[offset].get_label();
            for (uint j = 1; j < pack_width; ++j) {
              aggr_encoding =
                  Hash(aggr_encoding) ^ word[offset + j].get_label();
            }
            uint8_t row[ECPoint::byte_length];
            uint64_t row_idx = 0;
            for (uint j = 0; j < pack_width; ++j) {
              uint8_t bit = word[offset + j].get_label().LSB();
              row_idx |= (bit << j);
            }
            for (uint j = 0; j < (1UL << pack_width); ++j) {
              if (j == row_idx) {
                gc.read_data(row, ECPoint::byte_length);
              } else {
                gc.skip_data(ECPoint::byte_length);
              }
            }
            simd_word.encodings[i].point.dec(aggr_encoding, row);
          }
          break;
        }
        // fall through
      case MEASURE:  // or skip
      {
#ifdef MEASURE_STACK_COST
        measure_stack_flag = true;
#endif
        uint num_full_encoding = (width / max_pack_width) << max_pack_width;
        owner->get_gc().skip_ec_point(num_full_encoding);
        uint last_pack_width = width % max_pack_width;
        if (last_pack_width) {
          owner->get_gc().skip_ec_point(1UL << last_pack_width);
        }
#ifdef MEASURE_STACK_COST
        measure_stack_flag = false;
#endif
      } break;
      default:
        break;
    }
  }
  if (mode == GARBLE) {
    uint64_t num_tasks = 0;
    for (const SIMDWord& word : simd_words) {
      uint width = word.width();
      uint full_pack_count = width / max_pack_width;
      uint last_pack_width = width % max_pack_width;
      num_tasks += full_pack_count << max_pack_width;
      if (last_pack_width) {
        num_tasks += 1u << last_pack_width;
      }
    }
    std::vector<ECPoint> points(num_tasks);
    for (ECPoint& point : points) {
      point.initialize_temp_point();
    }
    std::vector<Label> merged_encodings(num_tasks);
    std::vector<uint64_t> buffer_idx_to_write(num_tasks);
    std::vector<BigInt> key_delta;

    key_delta.reserve(num_tasks);

    uint64_t curr_buffer_idx = 0;
    uint64_t task_id = 0;

    for (uint word_idx = 0; word_idx < num_words; ++word_idx) {
      const Word& word = words[word_idx];
      SIMDWord& simd_word = simd_words[word_idx];
      uint width = simd_word.width();
      for (uint i = 0; i < width; i += max_pack_width) {
        uint pack_width = std::min(max_pack_width, width - i);
        for (uint j = 0; j < (1u << pack_width); ++j) {
          uint64_t row_idx = 0;
          for (uint k = 0; k < pack_width; ++k) {
            const Label& bit_label = word[i + k].get_label();
            uint8_t bit_val = (j >> k) & 1;
            Label bit_encoding =
                (bit_val) ? bit_label ^ key_manager.get_Delta() : bit_label;
            if (k == 0) {
              merged_encodings[task_id] = bit_encoding;
            } else {
              merged_encodings[task_id] =
                  Hash(merged_encodings[task_id]) ^ bit_encoding;
            }
            uint8_t row_idx_bit = bit_val ^ word[i + k].get_label().LSB();
            row_idx |= (row_idx_bit << k);
          }
          buffer_idx_to_write[task_id] = curr_buffer_idx + row_idx;
          ++task_id;
        }
        curr_buffer_idx += 1UL << pack_width;
        const std::vector<BigInt>& gamma =
            key_manager.get_gamma(simd_word.bit_offset + i, 1UL << pack_width);

        for (uint j = 0; j < (1UL << pack_width); ++j) {
          key_delta.emplace_back(simd_word.label * gamma[j]);
          // from_montgomery is not thread-safe
          key_delta.back().from_montgomery();
        }
      }
    }
    Assert_eq(num_tasks, key_delta.size());
    ParTracker::start();
    if (!num_workers_g) {
#if NUM_THREADS > 1
#pragma omp parallel for num_threads(NUM_THREADS) \
    schedule(static) if (num_tasks > 4)
#endif
      for (uint64_t i = 0; i < num_tasks; ++i) {
        from_word_worker_g(key_delta[i], points[i]);
      }
    } else {
      uint64_t num_threads = num_workers_g + 1;  // include main thread
      uint64_t tasks_per_thread = (num_tasks + num_workers_g) / num_threads;
      for (uint i = 0; i < num_workers_g; ++i) {
        uint64_t begin = std::min(num_tasks, i * tasks_per_thread);
        uint64_t end = std::min(num_tasks, (i + 1) * tasks_per_thread);
        uint64_t num_job = end - begin;
        if (num_job == 0) {
          continue;
        }
        simd_word_worker_g_threads[i].from_word(
            num_job, key_delta.data() + begin, points.data() + begin);
      }
      uint64_t begin = std::min(num_tasks, num_workers_g * tasks_per_thread);
      uint64_t end =
          std::min(num_tasks, (num_workers_g + 1) * tasks_per_thread);
      for (uint64_t i = begin; i < end; ++i) {
        from_word_worker_g(key_delta[i], points[i]);
      }
      for (uint i = 0; i < num_workers_g; ++i) {
        simd_word_worker_g_threads[i].wait();
      }
    }
    ParTracker::stop();
    std::vector<uint64_t> buffer_idx_to_write_rev(num_tasks);
    for (uint64_t i = 0; i < num_tasks; ++i) {
      buffer_idx_to_write_rev[buffer_idx_to_write[i]] = i;
    }
    for (uint64_t i = 0; i < num_tasks; ++i) {
      uint8_t ec_point_ciphertext[ECPoint::byte_length];
      uint64_t from_idx = buffer_idx_to_write_rev[i];
      points[from_idx].enc(merged_encodings[from_idx], ec_point_ciphertext);
      gc.write_data(ec_point_ciphertext, ECPoint::byte_length);
    }
  }
  return simd_words;
}

std::vector<Word> SIMDWord::to_words(const std::vector<SIMDWord>& simd_words) {
  uint num_words = simd_words.size();
  Assert(num_words > 0);
  Mode mode = simd_words[0].get_mode();
  Gadget* owner = simd_words[0].get_owner();
  GCPtr& gc = owner->get_gc();
  std::vector<Word> words;
  words.reserve(num_words);
  if (mode == GARBLE) {
    std::vector<BigInt> exponents;
    static constexpr uint max_pack_width = SIMDWord::max_pack_width;
    struct MetaInfo {
      uint64_t buffer_idx_to_write;
      uint j;
    };

    uint64_t num_tasks = 0;
    uint64_t num_packs = 0;
    for (const SIMDWord& simd_word : simd_words) {
      uint full_pack_count =
          simd_word.width() / max_pack_width - simd_word.shift / max_pack_width;
      uint last_pack_width = simd_word.width() % max_pack_width;
      num_tasks += full_pack_count << max_pack_width;
      num_packs += full_pack_count;
      if (last_pack_width) {
        num_tasks += 1u << last_pack_width;
        ++num_packs;
      }
    }
    std::vector<MetaInfo> meta_info;
    std::vector<uint> pack_widths;
    exponents.reserve(num_tasks);
    meta_info.reserve(num_tasks);
    pack_widths.reserve(num_tasks);
    uint64_t curr_buffer_idx = 0;
    Assert(max_pack_width <= 3);
    std::vector<uint64_t> rand_perms(num_packs);
    // generate all randomness with one call
    secure_random((uint8_t*)rand_perms.data(), num_packs * sizeof(uint64_t));

    uint64_t pack_idx = 0;
    for (const SIMDWord& simd_word : simd_words) {
      uint i_begin = simd_word.shift / max_pack_width * max_pack_width;
      uint begin_offset = simd_word.shift % max_pack_width;
      for (uint i = i_begin; i < simd_word.width(); i += max_pack_width) {
        uint pack_width = std::min(max_pack_width, simd_word.width() - i);
        Assert(begin_offset <= pack_width);
        if (simd_word.bit_offsets_to_aggr.empty()) {
          // regular SIMD word
          const std::vector<BigInt>& gamma = key_manager.get_gamma(
              i + simd_word.bit_offset, 1UL << pack_width);
          for (uint j = 0; j < (1UL << pack_width); ++j) {
            exponents.emplace_back(
                (simd_word.label * gamma[j]).from_montgomery());
          }
        } else {
          // aggregated SIMD word
          uint offset = i + simd_word.bit_offsets_to_aggr[0];

          std::vector<BigInt> gamma =
              key_manager.get_gamma(offset, 1UL << pack_width);
          Assert_less((1UL << pack_width) - 1, gamma.size());
          BigInt extra_gamma = BigInt::zero().to_montgomery();
          for (uint j = 1; j < simd_word.bit_offsets_to_aggr.size(); ++j) {
            uint new_offset = i + simd_word.bit_offsets_to_aggr[j];
            extra_gamma += key_manager.get_gamma(new_offset, 1UL)[0];
          }
          for (uint j = 0; j < (1UL << pack_width); ++j) {
            gamma[j] += extra_gamma;
            exponents.emplace_back(
                (simd_word.label * gamma[j]).from_montgomery());
          }
        }
        uint64_t perm[1UL << max_pack_width];
        for (uint j = 0; j < (1UL << pack_width); ++j) {
          perm[j] = j;
        }
        permute_small(perm, perm + (1UL << pack_width), rand_perms[pack_idx++]);
        for (uint j = 0; j < (1UL << pack_width); ++j) {
          // from_montgomery is not thread-safe
          meta_info.push_back({curr_buffer_idx + perm[j], j});
          pack_widths.push_back(pack_width);
        }
        curr_buffer_idx += 1UL << pack_width;
      }
    }

    Assert_eq(num_tasks, exponents.size());
    Assert_eq(num_tasks, meta_info.size());
    std::vector<std::array<Label, max_pack_width>> hashes(num_tasks);
    std::vector<std::array<Label, max_pack_width>> hashes_pre_shuffle(
        num_tasks);
    std::vector<MAC> macs(num_tasks);
    std::vector<MAC> macs_pre_shuffle(num_tasks);
    std::vector<ECPoint> points(num_tasks);
    for (ECPoint& point : points) {
      point.initialize_temp_point();
    }
    ParTracker::start();
    if (!num_workers_g) {
#if NUM_THREADS > 1
#pragma omp parallel for num_threads(NUM_THREADS) \
    schedule(static) if (num_tasks > 4)
#endif
      for (uint64_t i = 0; i < num_tasks; ++i) {
        to_word_worker_g(hashes_pre_shuffle[i], pack_widths[i], exponents[i],
                         points[i], macs_pre_shuffle[i]);
      }
    } else {
      uint64_t num_threads = num_workers_g + 1;
      uint64_t tasks_per_thread = (num_tasks + num_workers_g) / num_threads;
      for (uint i = 0; i < num_workers_g; ++i) {
        uint64_t begin = std::min(num_tasks, i * tasks_per_thread);
        uint64_t end = std::min(num_tasks, (i + 1) * tasks_per_thread);
        uint64_t num_job = end - begin;
        if (num_job == 0) {
          continue;
        }
        simd_word_worker_g_threads[i].to_word(
            num_job, hashes_pre_shuffle.data() + begin,
            pack_widths.data() + begin, exponents.data() + begin,
            points.data() + begin, macs_pre_shuffle.data() + begin);
      }
      uint64_t begin = std::min(num_tasks, num_workers_g * tasks_per_thread);
      uint64_t end =
          std::min(num_tasks, (num_workers_g + 1) * tasks_per_thread);
      for (uint64_t i = begin; i < end; ++i) {
        to_word_worker_g(hashes_pre_shuffle[i], pack_widths[i], exponents[i],
                         points[i], macs_pre_shuffle[i]);
      }
      for (uint i = 0; i < num_workers_g; ++i) {
        simd_word_worker_g_threads[i].wait();
      }
    }
    ParTracker::stop();
    for (uint64_t i = 0; i < num_tasks; ++i) {
      uint64_t to_write_idx = meta_info[i].buffer_idx_to_write;
      uint pack_width = pack_widths[i];
      uint j = meta_info[i].j;
      for (uint k = 0; k < pack_width; ++k) {
        uint8_t bit = (j >> k) & 1;
        hashes[to_write_idx][k] = hashes_pre_shuffle[i][k];
        if (bit) {
          hashes[to_write_idx][k] ^= key_manager.get_Delta();
        }
      }
      macs[to_write_idx] = macs_pre_shuffle[i];
    }
    curr_buffer_idx = 0;
    for (const SIMDWord& simd_word : simd_words) {
      words.emplace_back(owner, simd_word.width() - simd_word.shift);
      uint i_begin = simd_word.shift / max_pack_width * max_pack_width;
      uint begin_offset = simd_word.shift % max_pack_width;
      for (uint i = i_begin; i < simd_word.width();
           i += simd_word.max_pack_width) {
        uint pack_width = std::min(max_pack_width, simd_word.width() - i);
        uint k_begin = (i == i_begin) ? begin_offset : 0;
        for (uint k = k_begin; k < pack_width; ++k) {
          words.back()[i + k - simd_word.shift] =
              Bit(owner, hashes[curr_buffer_idx][k], simd_word.values[i + k],
                  false, simd_word.pub_es[i + k]);
        }
        // although the gc has a buffer, it's still faster to first cache the
        // contents here and do a single write
        uint8_t write_buffer[(sizeof(MAC) + max_pack_width * sizeof(Label))
                             << max_pack_width];
        uint8_t* write_buffer_ptr = write_buffer;
        uint64_t macs_bytes = sizeof(MAC) * ((1UL << pack_width) - 1);
        memcpy(write_buffer_ptr, &macs[curr_buffer_idx + 1], macs_bytes);
        write_buffer_ptr += macs_bytes;
        for (uint j = 1; j < (1UL << pack_width); ++j) {
          for (uint k = k_begin; k < pack_width; ++k) {
            hashes[curr_buffer_idx + j][k] ^= hashes[curr_buffer_idx][k];
          }
          uint64_t labels_bytes = (pack_width - k_begin) * sizeof(Label);
          memcpy(write_buffer_ptr, hashes[curr_buffer_idx + j].data() + k_begin,
                 labels_bytes);
          write_buffer_ptr += labels_bytes;
        }
        gc.write_data(write_buffer, write_buffer_ptr - write_buffer);
        curr_buffer_idx += 1UL << pack_width;
      }
    }

  } else if (mode == EVAL) {
    std::vector<SIMDEncoding> encodings_to_compute;
    std::vector<uint> pack_widths;
    uint64_t num_tasks = 0;
    for (const SIMDWord& simd_word : simd_words) {
      if (simd_word.encodings.empty()) {
        continue;
      }
      num_tasks += (simd_word.width() + max_pack_width - 1) / max_pack_width -
                   simd_word.shift / max_pack_width;
    }
    encodings_to_compute.reserve(num_tasks);
    pack_widths.reserve(num_tasks);
    for (const SIMDWord& simd_word : simd_words) {
      if (simd_word.encodings.empty()) {
        continue;
      }
      for (uint encoding_idx = simd_word.shift / max_pack_width;
           encoding_idx * max_pack_width < simd_word.width(); ++encoding_idx) {
        uint pack_width = std::min(
            max_pack_width, simd_word.width() - encoding_idx * max_pack_width);
        const auto& encoding = simd_word.encodings[encoding_idx];
        encodings_to_compute.emplace_back(encoding);
        // from_montgomery is not thread-safe
        encodings_to_compute.back().from_montgomery();
        pack_widths.push_back(pack_width);
      }
    }
    Assert_eq(num_tasks, encodings_to_compute.size());
    Assert_eq(num_tasks, pack_widths.size());
    std::vector<std::array<Label, max_pack_width>> hashes(num_tasks);
    std::vector<MAC> macs(num_tasks);
    ParTracker::start();
    if (!num_workers_e) {
#if NUM_THREADS > 1
#pragma omp parallel for num_threads(NUM_THREADS) \
    schedule(static) if (num_tasks > 4)
#endif
      for (uint64_t i = 0; i < num_tasks; ++i) {
        to_word_worker_e(encodings_to_compute[i], pack_widths[i], hashes[i],
                         macs[i]);
      }
    } else {
      uint64_t num_threads = num_workers_e + 1;
      uint64_t tasks_per_thread = (num_tasks + num_workers_e) / num_threads;
      for (uint i = 0; i < num_workers_e; ++i) {
        uint64_t begin = std::min(num_tasks, i * tasks_per_thread);
        uint64_t end = std::min(num_tasks, (i + 1) * tasks_per_thread);
        simd_word_worker_e_threads[i].to_word(
            end - begin, encodings_to_compute.data() + begin,
            pack_widths.data() + begin, hashes.data() + begin,
            macs.data() + begin);
      }
      // main thread execute the rest
      uint64_t begin = std::min(num_tasks, num_workers_e * tasks_per_thread);
      uint64_t end = num_tasks;
      for (uint64_t i = 0; i < end - begin; ++i) {
        to_word_worker_e(encodings_to_compute[begin + i],
                         pack_widths[begin + i], hashes[begin + i],
                         macs[begin + i]);
      }

      for (uint i = 0; i < num_workers_e; ++i) {
        simd_word_worker_e_threads[i].wait();
      }
    }
    ParTracker::stop();
    uint64_t curr_buffer_idx = 0;

    for (const SIMDWord& simd_word : simd_words) {
      words.emplace_back(owner, simd_word.width() - simd_word.shift);
      uint i_begin = simd_word.shift / max_pack_width * max_pack_width;
      uint begin_offset = simd_word.shift % max_pack_width;
      for (uint i = i_begin; i < simd_word.width(); i += max_pack_width) {
        uint pack_width = std::min(max_pack_width, simd_word.width() - i);
        uint k_begin = (i == i_begin) ? begin_offset : 0;
        if (simd_word.encodings.empty()) {
          gc.skip_mac((1UL << pack_width) - 1);
          gc.skip_label(((1UL << pack_width) - 1) * (pack_width - k_begin));
          for (uint k = k_begin; k < pack_width; ++k) {
            words.back()[i + k - simd_word.shift] =
                Bit(owner, Label(), simd_word.values[i + k], false,
                    simd_word.pub_es[i + k], true);
          }
          continue;
        }
        uint64_t matched_row_idx = 0;
        MAC table_mac;
        for (uint64_t j = 1; j < (1UL << pack_width); ++j) {
          gc.read_mac(table_mac);
          if (table_mac == macs[curr_buffer_idx]) {
            matched_row_idx = j;
          }
        }
        Label rows[max_pack_width];
        for (uint64_t j = 1; j < (1UL << pack_width); ++j) {
          for (uint k = k_begin; k < pack_width; ++k) {
            if (j == matched_row_idx) {
              gc.read_label(rows[k]);
            } else {
              gc.skip_label();
            }
          }
          for (uint k = k_begin; k < pack_width; ++k) {
            words.back()[i + k - simd_word.shift] =
                Bit(owner, rows[k] ^ hashes[curr_buffer_idx][k],
                    simd_word.values[i + k], false, simd_word.pub_es[i + k]);
          }
        }

        ++curr_buffer_idx;
      }
    }
    Assert_eq(curr_buffer_idx, num_tasks);
  } else if (mode == MEASURE) {
#ifdef MEASURE_STACK_COST
    measure_stack_flag = true;
#endif
    for (const SIMDWord& simd_word : simd_words) {
      words.emplace_back(owner, simd_word.width() - simd_word.shift);
      uint i_begin = simd_word.shift / max_pack_width * max_pack_width;
      uint begin_offset = simd_word.shift % max_pack_width;
      for (uint i = i_begin; i < simd_word.width(); i += max_pack_width) {
        uint pack_width = std::min(max_pack_width, simd_word.width() - i);
        uint k_begin = (i == i_begin) ? begin_offset : 0;

        gc.skip_mac((1UL << pack_width) - 1);
        gc.skip_label(((1UL << pack_width) - 1) * (pack_width - k_begin));
        for (uint k = k_begin; k < pack_width; ++k) {
          words.back()[i + k - simd_word.shift] =
              Bit(owner, Label(), simd_word.values[i + k], false,
                  simd_word.pub_es[i + k], true);
        }
      }
    }
#ifdef MEASURE_STACK_COST
    measure_stack_flag = false;
#endif
  } else if (mode == DEBUG) {
    for (const SIMDWord& simd_word : simd_words) {
      words.emplace_back(owner, simd_word.width() - simd_word.shift);
      for (uint i = 0; i < words.back().width(); ++i) {
        words.back()[i] =
            Bit(owner, Label(), simd_word.values[i + simd_word.shift], false,
                simd_word.pub_es[i + simd_word.shift]);
      }
    }
  } else {
    Assert(false);
  }
  return words;
}

SIMDWord SIMDWord::slice(uint begin, uint end) const {
  SIMDWord result = *this;
  uint internal_end = end + result.shift;
  if (internal_end > result.width()) {
    internal_end = result.width();
  }
  result.values.resize(internal_end);
  result.pub_es.resize(internal_end);
  result.encodings.resize(internal_end);
  result >>= begin;
  return result;
}

// std::unordered_set<std::string> used_key_delta;

void SIMDWord::from_word_worker_g(const BigInt& key_delta, ECPoint& output) {
  // uint8_t key_delta_bytes[ECPoint::byte_length];
  // key_delta.to_bytes(key_delta_bytes);
  // std::string key_delta_str((char*)key_delta_bytes, ECPoint::byte_length);
  // if (used_key_delta.find(key_delta_str) != used_key_delta.end()) {
  //   std::cerr << "key_delta already used: " << key_delta << std::endl;
  //   exit(1);
  // }
  // used_key_delta.insert(key_delta_str);
  ECPoint::generator_pow(key_delta, output);
}

void SIMDWord::to_word_worker_g(std::array<Label, max_pack_width>& hashes,
                                uint pack_width, const BigInt& exponent,
                                ECPoint& tmp_point, MAC& mac) {
  // uint8_t key_delta_bytes[ECPoint::byte_length];
  // exponent.to_bytes(key_delta_bytes);
  // std::string key_delta_str((char*)key_delta_bytes, ECPoint::byte_length);
  // if (used_key_delta.find(key_delta_str) != used_key_delta.end()) {
  //   std::cerr << "exponent already used: " << exponent << std::endl;
  //   exit(1);
  // }
  // used_key_delta.insert(key_delta_str);
  ECPoint::generator_pow(exponent, tmp_point);
  hashes[0] = Hash(tmp_point);
  for (uint k = 1; k < pack_width; ++k) {
    hashes[k] = Hash(hashes[k - 1]);  // prng
  }
  mac = Mac(hashes[pack_width - 1]);
}

void SIMDWord::to_word_worker_e(
    SIMDEncoding& encoding_to_compute, uint pack_width,
    std::array<Label, max_pack_width>& output_hashes, MAC& output_mac) {
  output_hashes[0] = Hash(encoding_to_compute.simplify());
  for (uint k = 1; k < pack_width; ++k) {
    output_hashes[k] = Hash(output_hashes[k - 1]);
  }
  output_mac = Mac(output_hashes[pack_width - 1]);
}

uint SIMDWord::num_workers_g = 0;
uint SIMDWord::num_workers_e = 0;
SIMDWord::SIMDWordWorkerGThread* SIMDWord::simd_word_worker_g_threads = nullptr;
SIMDWord::SIMDWordWorkerEThread* SIMDWord::simd_word_worker_e_threads = nullptr;

}  // namespace PicoGRAM