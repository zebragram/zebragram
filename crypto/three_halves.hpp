#pragma once
#include "global.hpp"
#include "label.hpp"

namespace ZebraGRAM {
namespace ThreeHalves {
/**
 * @brief The garble function
 * nonce should increment by at least 3 for each call
 */
Label Garble(Label A0, Label B0, uint64_t nonce, uint8_t material[3 * 8 + 1]);

Label Eval(Label A, Label B, uint64_t nonce, const uint8_t material[3 * 8 + 1]);
}  // namespace ThreeHalves

}  // namespace ZebraGRAM