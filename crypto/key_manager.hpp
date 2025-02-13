#pragma once

#include "ec.hpp"
#include "label.hpp"

namespace PicoGRAM {
/**
 * @brief A global key manager for the garbler
 *
 */
struct KeyManager {
 private:
  Label Delta;  // Δ in free-XOR, i.e., the XOR of two random labels of the same
                // wire
  BigInt delta;  // δ in PicoGRAM, i.e., the discrete logarithm of the two
                 // encodings of the same sub-wire of a cable.
  Label rand_encoding;  // the encoding of all random input wires, intended to
                        // be shared with the evaluator (since each gate has a
                        // nonce, it's fine to reuse this encoding) Set to zero
                        // by default
  Label neg_rand_encoding;  // the negation of rand_encoding
  // Γ values in PicoGRAM for each cable bit offset
  std::vector<std::vector<BigInt>> Gamma;

 public:
  KeyManager() : Delta(Label::random_odd()), delta(BigInt::random()) {
    neg_rand_encoding = rand_encoding ^ Delta;
    delta.to_montgomery();
  }

  /**
   * @brief Get the γ value for a given bit offset
   *
   * @param offset the bit offset
   * @return const BigInt&
   */
  std::vector<BigInt> get_gamma(uint offset, uint64_t base = 1) {
    Assert_less(offset, (1UL << 24));
    if (offset >= Gamma.size()) {
      Gamma.reserve(offset + 1);
      while (offset >= Gamma.size()) {
        Gamma.emplace_back();
        std::vector<BigInt>& back = Gamma.back();
        back.reserve(base);
        back.push_back(BigInt::random());
      }
    }
    std::vector<BigInt>& result = Gamma[offset];
    while (result.size() < base) {
      Gamma[offset].push_back(result.back() + delta);
    }
    for (BigInt& bn : result) {
      bn.to_montgomery();
    }
    return Gamma[offset];
  }

  const Label& get_Delta() const { return Delta; }

  const Label& get_rand_encoding() const { return rand_encoding; }

  const Label& get_neg_rand_encoding() const { return neg_rand_encoding; }

  void set_Delta(const Label& Delta) {
    Assert(Delta.is_odd());
    this->Delta = Delta;
    neg_rand_encoding = rand_encoding ^ Delta;
  }
};

extern KeyManager key_manager;
}  // namespace PicoGRAM