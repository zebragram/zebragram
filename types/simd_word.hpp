#pragma once
#include <atomic>
#include <thread>
// #include <unordered_set>

#include "ec.hpp"
#include "lock.hpp"
#include "word.hpp"

namespace ZebraGRAM {
/**
 * @brief The garbling of cable in ZebraGRAM, i.e., a sequence of bits
 * garbled with a single local language.
 *
 */
struct SIMDWord : DataType {
  // The local language of the SIMD word. While the evaluator does not know the
  // local language, it can use this field to store the pending exponents of the
  // encodings and defer the exponentiation.
  BigInt label;

  uint shift = 0;

  // used by garbler to record the starting offset of the SIMD word in the
  // current time step. The offset is used to retrieve gamma for each bit.
  uint bit_offset;

  // used by the garbler to store bit offsets of merged SIMD words
  std::vector<uint> bit_offsets_to_aggr;

  // pack two bits into one ECPoint
  // change to 1 will further decrease the communication cost but
  // increases evaluator's computation
  static constexpr uint max_pack_width = 2;

  /**
   * @brief Used by the evaluator to record the encoding of each pack of bits in
   * a SIMD word. While in principle the encoding can be expressed as a single
   * ECPoint, it is much more efficient to defer EC operation mult-adds by first
   * operating on the exponents. The SIMDEncoding buffers the pending
   * exponentiation and aggregation of ECPoints.
   *
   */
  struct SIMDEncoding {
    // For non-aggregated bits
    ECPoint point;
    BigInt exp;
    // For aggregated bits
    std::vector<ECPoint> points_to_aggr;
    std::vector<BigInt> exps_to_aggr;
    // Whether the encoding is simplified to a single ECPoint
    bool simplified = false;

    SIMDEncoding() { simplified = true; }

    SIMDEncoding(const SIMDEncoding& other) = default;

    SIMDEncoding(SIMDEncoding&& other) = default;

    /**
     * @brief A scalar multiplication on the encoding but performed lazily
     *
     * @param exp
     */
    void pow(const BigInt& exp) {
      Assert(exp.is_set());
      if (exps_to_aggr.empty()) {
        // non-aggregated bits
        if (this->exp.is_set()) {
          this->exp *= exp;
        } else {
          this->exp = exp;
        }
      } else {
        // aggregated bits
        for (BigInt& e : exps_to_aggr) {
          if (e.is_set()) {
            e *= exp;
          } else {
            e = exp;
          }
        }
      }
      simplified = false;
    }

    /**
     * @brief Convert the exponents to standard form
     *
     */
    void from_montgomery() {
      if (simplified) {
        return;
      }
      if (exps_to_aggr.empty()) {
        if (exp.is_set()) {
          exp.from_montgomery();
        }
      } else {
        for (BigInt& e : exps_to_aggr) {
          if (e.is_set()) {
            e.from_montgomery();
          }
        }
      }
    }

    /**
     * @brief Simplify the encoding to a single ECPoint
     *
     * @return const ECPoint& the simplified ECPoint
     */
    const ECPoint& simplify() {
      if (simplified) {
        return point;
      }
      if (points_to_aggr.empty()) {
        // non-aggregated bits
        Assert(exps_to_aggr.empty());
        point.pow(exp, point);
        exp.unset();
      } else {
        // aggregated bits
        Assert_eq(points_to_aggr.size(), exps_to_aggr.size());
        ECPoint::aggr_pow(exps_to_aggr, points_to_aggr, point);
        points_to_aggr.clear();
        exps_to_aggr.clear();
      }
      simplified = true;
      return point;
    }

    /**
     * @brief Add an encoding to the aggregation list
     *
     * @param encoding
     */
    SIMDEncoding& operator+=(const SIMDEncoding& encoding) {
      Assert(!exp.is_set());  // the ec_point and exp of this should be unset
      // do not support nested aggregation
      Assert(encoding.points_to_aggr.empty());
      points_to_aggr.push_back(encoding.point);
      exps_to_aggr.push_back(encoding.exp);
      simplified = false;
      return *this;
    }
  };

  // used by the evaluator to store the encodings of the SIMD word, the first
  // div second BigInt in the pair is the pending exponent for each bit
  std::vector<SIMDEncoding> encodings;

  // the values of the bits
  std::vector<bool> values;

  // whether the bits are public to the evaluator
  std::vector<bool> pub_es;

 private:
  static void from_word_worker_g(const BigInt& key_delta, ECPoint& output);

  static void to_word_worker_g(std::array<Label, max_pack_width>& hashes,
                               uint pack_width, const BigInt& exponent,
                               ECPoint& tmp_point, MAC& mac);

  static void to_word_worker_e(SIMDEncoding& encoding_to_compute,
                               uint pack_width,
                               std::array<Label, max_pack_width>& output_hashes,
                               MAC& output_mac);

 public:
  // the data structure is not supposed to be thread-safe
  struct SIMDWordWorkerGThread : BusyThread {
    bool is_to_word;
    uint job_count = 0;

    uint* num_hashes_to_merge;

    std::array<Label, SIMDWord::max_pack_width>* hashes = NULL;

    BigInt* exponent = NULL;
    ECPoint* output = NULL;
    MAC* mac = NULL;

    void exec_step() override {
      for (uint i = 0; i < job_count; ++i) {
        if (is_to_word) {
          SIMDWord::to_word_worker_g(hashes[i], num_hashes_to_merge[i],
                                     exponent[i], output[i], mac[i]);
        } else {
          SIMDWord::from_word_worker_g(exponent[i], output[i]);
        }
      }
    }

    void to_word(uint job_count, std::array<Label, max_pack_width>* hashes,
                 uint* pack_width, BigInt* exponent, ECPoint* tmp_point,
                 MAC* mac) {
      Assert(!stopped.load(std::memory_order_acquire));
      Assert(!has_job.load(std::memory_order_acquire));
      this->job_count = job_count;
      this->is_to_word = true;
      this->hashes = hashes;
      this->num_hashes_to_merge = pack_width;
      this->exponent = exponent;
      this->output = tmp_point;
      this->mac = mac;

      has_job.store(true, std::memory_order_release);
    }

    void from_word(uint job_count, BigInt* key_delta, ECPoint* output) {
      Assert(!stopped.load(std::memory_order_acquire));
      Assert(!has_job.load(std::memory_order_acquire));
      this->job_count = job_count;
      this->is_to_word = false;
      this->exponent = key_delta;
      this->output = output;

      has_job.store(true, std::memory_order_release);
    }
  };

  // the data structure is not supposed to be thread-safe
  struct SIMDWordWorkerEThread : BusyThread {
    uint job_count = 0;
    SIMDEncoding* encoding_to_compute = nullptr;
    uint* pack_width = nullptr;
    std::array<Label, max_pack_width>* output_hashes = nullptr;
    MAC* output_mac = nullptr;

    void exec_step() override {
      for (uint i = 0; i < job_count; ++i) {
        SIMDWord::to_word_worker_e(encoding_to_compute[i], pack_width[i],
                                   output_hashes[i], output_mac[i]);
      }
    }

    void to_word(uint job_count, SIMDEncoding* encoding_to_compute,
                 uint* pack_width,
                 std::array<Label, max_pack_width>* output_hashes,
                 MAC* output_mac) {
      Assert(!stopped.load(std::memory_order_acquire));
      Assert(!has_job.load(std::memory_order_acquire));
      this->job_count = job_count;
      this->encoding_to_compute = encoding_to_compute;
      this->pack_width = pack_width;
      this->output_hashes = output_hashes;
      this->output_mac = output_mac;

      has_job.store(true, std::memory_order_release);
    }
  };

  static uint num_workers_g;
  static uint num_workers_e;
  static SIMDWordWorkerGThread* simd_word_worker_g_threads;
  static SIMDWordWorkerEThread* simd_word_worker_e_threads;

  // notice that even if the garbler knows the values of a SIMD word, the values
  // cannot not be preserved after the SIMD word is translated through a SIMD
  // link. Therefore, we do not include the pub_g field.

  SIMDWord() = default;

  /**
   * @brief Mock a SIMDWord word
   *
   * @param owner owner gadget
   * @param width the width of the SIMD word
   * @param bit_offset the bit offset of the SIMD word
   * @param label the local language of the SIMD word
   * @param pub_es (optional) the public flags for the evaluator
   */
  SIMDWord(Gadget* owner, uint width, uint64_t bit_offset, const BigInt& label,
           const std::vector<bool>& pub_es = std::vector<bool>());

  SIMDWord(const Word& word, uint64_t bit_offset) : DataType(word.get_owner()) {
    from_word(word, bit_offset);
  }

  static SIMDWord new_aggr_word(const SIMDWord& simd_word) {
    SIMDWord result(simd_word.get_owner(), 0, 0, BigInt());
    result.label = simd_word.label;
    result.shift = simd_word.shift;
    result.pub_es = simd_word.pub_es;
    result.values = simd_word.values;
    Mode mode = simd_word.get_mode();
    if (mode == GARBLE) {
      result.bit_offsets_to_aggr.push_back(simd_word.bit_offset);
    } else if (mode == EVAL) {
      uint num_encodings = simd_word.encodings.size();
      result.encodings.resize(num_encodings);
      for (uint i = 0; i < num_encodings; ++i) {
        result.encodings[i] += simd_word.encodings[i];
      }
    }

    return result;
  }

  void aggregate_with(const SIMDWord& other) {
    Mode mode = get_mode();
    uint64_t num_encodings = encodings.size();
    if ((mode == EVAL && num_encodings) || mode == GARBLE) {
      Assert_eq(label, other.label);
    }
    Assert_eq(mode, other.get_mode());
    Assert_eq(get_owner(), other.get_owner());

    Assert_eq(shift, other.shift);
    Assert_eq(width(), other.width());
    for (uint i = 0; i < width(); ++i) {
      pub_es[i] = pub_es[i] && other.pub_es[i];
      values[i] = values[i] || other.values[i];
    }
    Assert_eq(num_encodings, other.encodings.size());

    if (mode == GARBLE) {
      bit_offsets_to_aggr.push_back(other.bit_offset);
    } else if (mode == EVAL) {
      for (uint64_t i = 0; i < num_encodings; ++i) {
        encodings[i] += other.encodings[i];
      }
    }
  }

  static std::vector<SIMDWord> aggr_simd_words(
      const std::vector<SIMDWord>& simd_words, uint out_word_count) {
    Assert(out_word_count);
    uint aggr_group_size = simd_words.size() / out_word_count;
    Assert(aggr_group_size * out_word_count == simd_words.size());
    std::vector<SIMDWord> result_words;
    result_words.reserve(out_word_count);
    for (uint i = 0; i < out_word_count; ++i) {
      result_words.emplace_back(new_aggr_word(simd_words[i]));
    }
    for (uint i = 1; i < aggr_group_size; ++i) {
      for (uint j = 0; j < out_word_count; ++j) {
        result_words[j].aggregate_with(simd_words[i * out_word_count + j]);
      }
    }
    return result_words;
  }

  static void start_workers_g(uint num_workers) {
    if (!num_workers) {
      return;
    }
    num_workers_g = num_workers;
    simd_word_worker_g_threads = new SIMDWordWorkerGThread[num_workers];
    // spawn a thread for each worker with omp nowait

    for (uint i = 0; i < num_workers; ++i) {
      simd_word_worker_g_threads[i].start();
    }
  }

  static void stop_workers_g() {
    if (!simd_word_worker_g_threads) {
      return;
    }
    for (uint i = 0; i < num_workers_g; ++i) {
      simd_word_worker_g_threads[i].stop();
    }
    num_workers_g = 0;
    delete[] simd_word_worker_g_threads;
  }

  static void start_workers_e(uint num_workers) {
    if (!num_workers) {
      return;
    }
    num_workers_e = num_workers;
    simd_word_worker_e_threads = new SIMDWordWorkerEThread[num_workers];
    // spawn a thread for each worker with omp nowait

    for (uint i = 0; i < num_workers; ++i) {
      simd_word_worker_e_threads[i].start();
    }
  }

  static void stop_workers_e() {
    if (!simd_word_worker_e_threads) {
      return;
    }
    for (uint i = 0; i < num_workers_e; ++i) {
      simd_word_worker_e_threads[i].stop();
    }
    num_workers_e = 0;
    delete[] simd_word_worker_e_threads;
  }

  /**
   * @brief Construct a SIMD word from a Word. The bit offset is set
   * according to the bit_offset of the owner.
   *
   * @param word the Word to be converted
   */
  void from_word(const Word& word, uint64_t bit_offset);

  static std::vector<SIMDWord> from_words(const std::vector<Word>& words,
                                          uint64_t begin_bit_offset);

  static std::vector<Word> to_words(const std::vector<SIMDWord>& simd_words);

  /**
   * @brief Convert the SIMD word to a Word
   *
   * @return Word
   */
  Word to_word() const;

  // /**
  //  * @brief Get the bit at a given index, or a constant bit 0 if the index is
  //  * out of range
  //  *
  //  * @param i the index of the bit
  //  * @return Bit
  //  */
  // Bit operator[](uint i) const;

  uint width() const { return values.size(); }

  void operator>>=(uint shift) {
    // logical right shift
    this->shift += shift;
    Assert(this->shift <= values.size());
  }

  SIMDWord slice(uint begin, uint end) const;

  /**
   * @brief Print the encodings of the SIMD word. This function is for debugging
   * purposes.
   *
   */
  void print_encodings() const {
    if (get_mode() == GARBLE) {
      for (uint i = 0; i < width(); ++i) {
        const std::vector<BigInt>& gamma =
            key_manager.get_gamma(bit_offset + i, 2);
        std::cout << ECPoint::generator_pow(label * gamma[0]) << " / "
                  << ECPoint::generator_pow(label * gamma[1]) << std::endl;
      }

    } else if (get_mode() == EVAL) {
      for (const auto& encoding : encodings) {
        auto encoding_copy = encoding;
        encoding_copy.simplify();
        std::cout << encoding_copy.point << std::endl;
      }
    }
  }
};
}  // namespace ZebraGRAM