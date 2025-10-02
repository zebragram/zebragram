#pragma once
#include <algorithm>
#include <vector>

#include "arith_label.hpp"
#include "arith_mul.hpp"
#include "bit.hpp"
#include "data_type.hpp"
#include "label.hpp"
#include "mod2.hpp"

namespace ZebraGRAM {

struct Word;
/**
 * @brief A data type that represents a word, which is a sequence of bits. A
 * word < 64 bits can also represent an unsigned integer.
 *
 */
struct ArithWord : DataType {
  std::vector<ArithLabel> payloads;

  ArithWord() : DataType() {}
  ArithWord(Gadget *owner) : DataType(owner) {}

  static ArithWord input_dbg(Gadget *owner, std::vector<uint64_t> values);
  static ArithWord constant(Gadget *owner, size_t width, uint64_t value);

  std::vector<ArithLabel> reveal() const;

  template <typename Iterator>
  void set_payload(Iterator begin, Iterator end) {
    // initialize payload vector using the range constructor
    payloads = std::vector<ArithLabel>(begin, end);
  }

  /**
   * @brief Translate a word src to a word dst. In garble mode, both the src and
   * dst label need to be set. In eval mode, the src bit is converted to the dst
   * bit. A usecase is to break down a circuit into multiple sub-circuits to
   * garble, and use join between sub-circuits to connect them.
   * Incurs lambda bits of communication per bit.
   *
   * @param src the word to be translated from
   * @param dst the word to be translated to
   */
  static void join(const ArithWord &src, ArithWord &dst);

  /**
   * @brief The width (number of binary bits) of the word
   *
   * @return uint
   */
  uint width() const { return payloads.size(); }

  ArithWord operator+(const ArithWord &other) const {
    ArithWord result(get_owner() ? get_owner() : other.get_owner());
    Assert(!get_owner() || !other.get_owner() ||
           get_owner() == other.get_owner());

    Assert(width() == other.width() || payloads.empty() ||
           other.payloads.empty());
    if (!payloads.empty() && !other.payloads.empty()) {
      result.payloads.resize(width());
      for (size_t i = 0; i < width(); ++i) {
        result.payloads[i] = payloads[i] + other.payloads[i];
      }
    } else if (!payloads.empty()) {
      result.payloads = payloads;
    } else if (!other.payloads.empty()) {
      result.payloads = other.payloads;
    }
    return result;
  }

  // only for payloads
  ArithWord operator-(const ArithWord &other) const {
    ArithWord result(get_owner() ? get_owner() : other.get_owner());
    Assert(!get_owner() || !other.get_owner() ||
           get_owner() == other.get_owner());

    Assert(width() == other.width() || payloads.empty() ||
           other.payloads.empty());
    if (!payloads.empty() && !other.payloads.empty()) {
      result.payloads.resize(width());
      for (size_t i = 0; i < width(); ++i) {
        result.payloads[i] = payloads[i] - other.payloads[i];
      }
    } else if (!payloads.empty()) {
      result.payloads = payloads;
    } else if (!other.payloads.empty()) {
      result.payloads.resize(other.width());
      for (size_t i = 0; i < other.width(); ++i) {
        result.payloads[i] = -other.payloads[i];
      }
    }
    return result;
  }

  ArithWord operator*(uint64_t scalar) const {
    ArithWord result(get_owner());
    result.payloads.resize(width());
    for (size_t i = 0; i < width(); ++i) {
      fmpz_mul_ui(result.payloads[i].auth, payloads[i].auth, scalar);
      fmpz_fdiv_r_2exp(result.payloads[i].auth, result.payloads[i].auth,
                       LEN_AUTH_SHARE);
      fmpz_mul_ui(result.payloads[i].raw, payloads[i].raw, scalar);
      fmpz_fdiv_r_2exp(result.payloads[i].raw, result.payloads[i].raw,
                       LEN_RAW_SHARE);
    }
    return result;
  }

  ArithWord &operator*=(uint64_t scalar) {
    *this = *this * scalar;
    return *this;
  }
  ArithWord &operator+=(const ArithWord &other) {
    *this = *this + other;
    return *this;
  }
  ArithWord &operator-=(const ArithWord &other) {
    *this = *this - other;
    return *this;
  }

  static ArithWord mux(const Bit &bit, const ArithWord &word1,
                       const ArithWord &word2);

  static ArithWord prng_from_label(Label seed, size_t width);

  static ArithWord rand_label_arith_word(Gadget *owner, size_t width) {
    ArithWord word(owner);
    word.payloads.resize(width);
    for (size_t i = 0; i < width; ++i) {
      secure_random_fmpz(word.payloads[i].auth, LEN_AUTH_SHARE);
      secure_random_fmpz(word.payloads[i].raw, LEN_RAW_SHARE);
    }
    return word;
  }

  Word decompose_all(const Mod2Ctx &mod2, uint num_bits,
                     const PaillierPrivKey &sk, const PaillierPubKey &pk) const;

  bool check() const {
    for (const auto &label : payloads) {
      if (!label.check()) {
        return false;
      }
    }
    return true;
  }
};

void bit_arith_word_mul(const Bit &bit, const ArithWord &word,
                        ArithWord &result);

// todo: unify the style
std::vector<ArithWord> batch_bit_arith_word_mul(
    const std::vector<Bit> &bits, const std::vector<ArithWord> &words);

// requires bits to be one-hot
void batch_swap_arith_words(const std::vector<Bit> &bits,
                            std::vector<ArithWord> &words0, ArithWord &word1);

void batch_swap_arith_words(const std::vector<Bit> &bits,
                            std::vector<ArithWord> &words0,
                            std::vector<ArithWord> &words1);
}  // namespace ZebraGRAM
