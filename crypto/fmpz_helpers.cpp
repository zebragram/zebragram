#include "fmpz_helpers.hpp"
#include <sstream>
#include <stdexcept>
#include <cstdlib>  // For malloc/free

namespace flint_ipcl {

BigNumber fmpzToBigNumber(const fmpz_t fmpz_val) {
#ifdef IPCL_AVAILABLE
    // Handle zero case
    if (fmpz_is_zero(fmpz_val)) {
        return BigNumber(0);
    }
    
    // Handle negative numbers (note: BigNumber typically handles positive values)
    bool is_negative = fmpz_sgn(fmpz_val) < 0;
    fmpz_t abs_val;
    fmpz_init(abs_val);
    fmpz_abs(abs_val, fmpz_val);
    
    // Use FLINT's direct conversion to array of limbs
    slong size = fmpz_size(abs_val);
    if (size <= 0) {
        fmpz_clear(abs_val);
        return BigNumber(0);
    }
    
    // Get FLINT limbs directly
    ulong* limbs = (ulong*)malloc(size * sizeof(ulong));
    if (!limbs) {
        fmpz_clear(abs_val);
        throw std::runtime_error("Memory allocation failed");
    }
    
    fmpz_get_ui_array(limbs, size, abs_val);
    
    // Convert FLINT limbs (64-bit) to 32-bit words for IPCL
    std::vector<Ipp32u> words;
    words.reserve(size * 2);
    
    for (slong i = 0; i < size; ++i) {
        ulong limb = limbs[i];
        words.push_back(static_cast<Ipp32u>(limb & 0xFFFFFFFFUL));
        words.push_back(static_cast<Ipp32u>(limb >> 32));
    }
    
    // Remove trailing zeros from words
    while (!words.empty() && words.back() == 0) {
        words.pop_back();
    }
    
    if (words.empty()) {
        free(limbs);
        fmpz_clear(abs_val);
        return BigNumber(0);
    }
    
    // Create BigNumber using the constructor
    IppsBigNumSGN sgn = is_negative ? IppsBigNumNEG : IppsBigNumPOS;
    BigNumber result(words.data(), words.size(), sgn);
    
    free(limbs);
    fmpz_clear(abs_val);
    return result;
#else
    // Fallback to string conversion when IPCL is not available
    char* str = fmpz_get_str(nullptr, 16, fmpz_val);
    if (!str) {
        throw std::runtime_error("Failed to convert fmpz_t to string");
    }
    
    std::string hex_str(str);
    bool is_negative = false;
    if (hex_str[0] == '-') {
        is_negative = true;
        hex_str = hex_str.substr(1);
    }
    
    if (hex_str.substr(0, 2) != "0x") {
        hex_str = "0x" + hex_str;
    }
    
    BigNumber result(hex_str.c_str());
    flint_free(str);
    return result;
#endif
}

void bigNumberToFmpz(const BigNumber& bn, fmpz_t fmpz_val) {
#ifdef IPCL_AVAILABLE
    // Use the efficient num2vec method to get 32-bit words
    std::vector<Ipp32u> words;
    bn.num2vec(words);
    
    if (words.empty()) {
        fmpz_zero(fmpz_val);
        return;
    }
    
    // Convert 32-bit words to 64-bit FLINT limbs
    size_t limb_count = (words.size() + 1) / 2; // Each limb holds 2 words
    ulong* limbs = (ulong*)malloc(limb_count * sizeof(ulong));
    if (!limbs) {
        throw std::runtime_error("Memory allocation failed");
    }
    
    // Pack 32-bit words into 64-bit limbs
    for (size_t i = 0; i < limb_count; ++i) {
        ulong limb = 0;
        
        // Low 32 bits
        if (i * 2 < words.size()) {
            limb |= static_cast<ulong>(words[i * 2]);
        }
        
        // High 32 bits
        if (i * 2 + 1 < words.size()) {
            limb |= static_cast<ulong>(words[i * 2 + 1]) << 32;
        }
        
        limbs[i] = limb;
    }
    
    // Convert limbs to FLINT fmpz
    fmpz_set_ui_array(fmpz_val, limbs, (slong)limb_count);
    free(limbs);
#else
    // Fallback to string conversion when IPCL is not available
    std::string hex_str;
    bn.num2hex(hex_str);
    
    if (hex_str.substr(0, 2) == "0x") {
        hex_str = hex_str.substr(2);
    }
    
    int result = fmpz_set_str(fmpz_val, hex_str.c_str(), 16);
    if (result != 0) {
        throw std::runtime_error("Failed to convert hex string to fmpz_t");
    }
#endif
}

std::vector<BigNumber> fmpzArrayToBigNumberVector(
    const fmpz_t* fmpz_array, size_t size) {
    std::vector<BigNumber> result;
    result.resize(size);
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < size; ++i) {
        result[i] = fmpzToBigNumber(fmpz_array[i]);
    }
    
    return result;
}

std::vector<BigNumber> fmpzVectorToBigNumberVector(
    const std::vector<fmpz_t>& fmpz_vec) {
    size_t size = fmpz_vec.size();
    return fmpzArrayToBigNumberVector(fmpz_vec.data(), size);
}



void bigNumberVectorToFmpzVector(
    const std::vector<BigNumber>& bn_vec,
    std::vector<fmpz_t>& fmpz_vec) {
    if (fmpz_vec.size() != bn_vec.size()) {
        throw std::invalid_argument("Vector sizes must match");
    }
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < bn_vec.size(); ++i) {
        bigNumberToFmpz(bn_vec[i], fmpz_vec[i]);
    }
}

} // namespace flint_ipcl