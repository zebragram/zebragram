#pragma once
#include <flint/fmpz.h>
#include <omp.h>

#include <vector>
#include <string>
#include <flint/fmpz.h>

// Forward declaration for BigNumber when IPCL is not available
#ifdef IPCL_AVAILABLE
#include <ipcl/ipcl.hpp>
using Ipp32u = unsigned int;
#else
// Dummy BigNumber class for compilation when IPCL is not available
class BigNumber {
public:
    BigNumber() = default;
    BigNumber(int val) {}
    BigNumber(const char* hex_str) {}
    
    // Dummy methods to satisfy interface
    bool create(const void* pData, int length, int sgn = 0) { return false; }
    void num2vec(std::vector<unsigned int>& v) const {}
    void num2hex(std::string& hex_str) const { hex_str = "0x0"; }
    std::string toHexString() const { return "0x0"; }
};
using Ipp32u = unsigned int;
enum { IppsBigNumPOS = 0, IppsBigNumNEG = 1 };
using IppsBigNumSGN = int;
#endif

namespace flint_ipcl {

/**
 * Convert FLINT fmpz_t to IPCL BigNumber
 * @param[in] fmpz_val FLINT integer
 * @return IPCL BigNumber
 */
BigNumber fmpzToBigNumber(const fmpz_t fmpz_val);

/**
 * Convert IPCL BigNumber to FLINT fmpz_t
 * @param[in] bn IPCL BigNumber
 * @param[out] fmpz_val FLINT integer (must be initialized)
 */
void bigNumberToFmpz(const BigNumber& bn, fmpz_t fmpz_val);

/**
 * Convert vector of FLINT fmpz_t to vector of IPCL BigNumber
 * @param[in] fmpz_vec vector of FLINT integers
 * @return vector of IPCL BigNumbers
 */
std::vector<BigNumber> fmpzVectorToBigNumberVector(
    const std::vector<fmpz_t>& fmpz_vec);

/**
 * Convert vector of FLINT fmpz_t to vector of IPCL BigNumber (C-style array)
 * @param[in] fmpz_array C-style array of FLINT integers
 * @param[in] size size of the array
 * @return vector of IPCL BigNumbers
 */
std::vector<BigNumber> fmpzArrayToBigNumberVector(
    const fmpz_t* fmpz_array, size_t size);

/**
 * Convert vector of IPCL BigNumber to vector of FLINT fmpz_t
 * @param[in] bn_vec vector of IPCL BigNumbers
 * @param[out] fmpz_vec vector of FLINT integers (must be pre-initialized)
 */
void bigNumberVectorToFmpzVector(
    const std::vector<BigNumber>& bn_vec,
    std::vector<fmpz_t>& fmpz_vec);

} // namespace flint_ipcl