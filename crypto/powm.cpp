#include "powm.hpp"
#ifdef IPCL_AVAILABLE
#include "fmpz_helpers.hpp"
#endif

void powm_cpu(fmpz_t result, const fmpz_t base, const fmpz_t exp,
              const fmpz_t mod) {
#ifdef IPCL_AVAILABLE
  try {
    // Convert FLINT integers to IPCL BigNumbers
    BigNumber bn_base = flint_ipcl::fmpzToBigNumber(base);
    BigNumber bn_exp = flint_ipcl::fmpzToBigNumber(exp);
    BigNumber bn_mod = flint_ipcl::fmpzToBigNumber(mod);

    // Perform modular exponentiation using IPCL
    BigNumber bn_result = ipcl::ippModExp(bn_base, bn_exp, bn_mod);

    // Convert result back to FLINT
    flint_ipcl::bigNumberToFmpz(bn_result, result);
  } catch (...) {
    // Fallback to FLINT if IPCL fails
    fmpz_powm(result, base, exp, mod);
  }
#else
  fmpz_powm(result, base, exp, mod);
#endif
}

void powm_cpu_batch(fmpz_t *results, const fmpz_t *bases, const fmpz_t exp,
                    const fmpz_t mod, int batch_size) {
#ifdef IPCL_AVAILABLE
  try {
    // Convert arrays to IPCL BigNumber vectors
    std::vector<BigNumber> bn_bases =
        flint_ipcl::fmpzArrayToBigNumberVector(bases, batch_size);
    BigNumber bn_exp = flint_ipcl::fmpzToBigNumber(exp);
    BigNumber bn_mod = flint_ipcl::fmpzToBigNumber(mod);

    // Create vector of same exponent and modulus for batch operation
    std::vector<BigNumber> bn_exps(batch_size, bn_exp);
    std::vector<BigNumber> bn_mods(batch_size, bn_mod);

    // Perform batch modular exponentiation using IPCL
    std::vector<BigNumber> bn_results =
        ipcl::ippModExp(bn_bases, bn_exps, bn_mods);

// Convert results back to FLINT
#pragma omp parallel for schedule(static)
    for (int i = 0; i < batch_size; i++) {
      flint_ipcl::bigNumberToFmpz(bn_results[i], results[i]);
    }
  } catch (...) {
    // Fallback to FLINT if IPCL fails
#pragma omp parallel for schedule(static)
    for (int i = 0; i < batch_size; i++) {
      fmpz_powm(results[i], bases[i], exp, mod);
    }
  }
#else
#pragma omp parallel for schedule(static)
  for (int i = 0; i < batch_size; i++) {
    fmpz_powm(results[i], bases[i], exp, mod);
  }
#endif
}

void powm_cpu_batch(fmpz_t *results, const fmpz_t base, const fmpz_t *exps,
                    const fmpz_t mod, int batch_size) {
#ifdef IPCL_AVAILABLE
  try {
    // Convert arrays to IPCL BigNumber vectors
    BigNumber bn_base = flint_ipcl::fmpzToBigNumber(base);
    std::vector<BigNumber> bn_exps =
        flint_ipcl::fmpzArrayToBigNumberVector(exps, batch_size);
    BigNumber bn_mod = flint_ipcl::fmpzToBigNumber(mod);

    // Create vector of same base and modulus for batch operation
    std::vector<BigNumber> bn_bases(batch_size, bn_base);
    std::vector<BigNumber> bn_mods(batch_size, bn_mod);

    // Perform batch modular exponentiation using IPCL
    std::vector<BigNumber> bn_results =
        ipcl::ippModExp(bn_bases, bn_exps, bn_mods);
// Convert results back to FLINT
#pragma omp parallel for schedule(static)
    for (int i = 0; i < batch_size; i++) {
      flint_ipcl::bigNumberToFmpz(bn_results[i], results[i]);
    }
  } catch (...) {
    // Fallback to FLINT if IPCL fails
#pragma omp parallel for schedule(static)
    for (int i = 0; i < batch_size; i++) {
      fmpz_powm(results[i], base, exps[i], mod);
    }
  }
#else
#pragma omp parallel for schedule(static)
  for (int i = 0; i < batch_size; i++) {
    fmpz_powm(results[i], base, exps[i], mod);
  }
#endif
}

void powm_cpu_batch(fmpz_t *results, const fmpz_t *base, const fmpz_t *exps,
                    const fmpz_t mod, int batch_size) {
#ifdef IPCL_AVAILABLE
  try {
    // Convert arrays to IPCL BigNumber vectors
    std::vector<BigNumber> bn_bases =
        flint_ipcl::fmpzArrayToBigNumberVector(base, batch_size);
    std::vector<BigNumber> bn_exps =
        flint_ipcl::fmpzArrayToBigNumberVector(exps, batch_size);
    BigNumber bn_mod = flint_ipcl::fmpzToBigNumber(mod);

    // Create vector of same modulus for batch operation
    std::vector<BigNumber> bn_mods(batch_size, bn_mod);

    // Perform batch modular exponentiation using IPCL
    std::vector<BigNumber> bn_results =
        ipcl::ippModExp(bn_bases, bn_exps, bn_mods);

// Convert results back to FLINT
#pragma omp parallel for schedule(static)
    for (int i = 0; i < batch_size; i++) {
      flint_ipcl::bigNumberToFmpz(bn_results[i], results[i]);
    }
  } catch (...) {
    // Fallback to FLINT if IPCL fails
#pragma omp parallel for schedule(static)
    for (int i = 0; i < batch_size; i++) {
      fmpz_powm(results[i], base[i], exps[i], mod);
    }
  }
#else
#pragma omp parallel for schedule(static)
  for (int i = 0; i < batch_size; i++) {
    fmpz_powm(results[i], base[i], exps[i], mod);
  }
#endif
}

void powm_cpu_batch(fmpz_t *results, const fmpz_t **base, const fmpz_t **exps,
                    const fmpz_t mod, int batch_size) {
#ifdef IPCL_AVAILABLE
  try {
    // Convert pointer arrays to IPCL BigNumber vectors
    std::vector<BigNumber> bn_bases;
    std::vector<BigNumber> bn_exps;
    bn_bases.reserve(batch_size);
    bn_exps.reserve(batch_size);

    for (int i = 0; i < batch_size; i++) {
      bn_bases.push_back(flint_ipcl::fmpzToBigNumber(*(base[i])));
      bn_exps.push_back(flint_ipcl::fmpzToBigNumber(*(exps[i])));
    }

    BigNumber bn_mod = flint_ipcl::fmpzToBigNumber(mod);

    // Create vector of same modulus for batch operation
    std::vector<BigNumber> bn_mods(batch_size, bn_mod);

    // Perform batch modular exponentiation using IPCL
    std::vector<BigNumber> bn_results =
        ipcl::ippModExp(bn_bases, bn_exps, bn_mods);

// Convert results back to FLINT
#pragma omp parallel for schedule(static)
    for (int i = 0; i < batch_size; i++) {
      flint_ipcl::bigNumberToFmpz(bn_results[i], results[i]);
    }
  } catch (...) {
    // Fallback to FLINT if IPCL fails
#pragma omp parallel for schedule(static)
    for (int i = 0; i < batch_size; i++) {
      fmpz_powm(results[i], *(base[i]), *(exps[i]), mod);
    }
  }
#else
#pragma omp parallel for schedule(static)
  for (int i = 0; i < batch_size; i++) {
    fmpz_powm(results[i], *(base[i]), *(exps[i]), mod);
  }
#endif
}

void powm_precomputed_cpu(fmpz_t result, const fmpz_t exp,
                          const PowmPrecomputeTable &table) {
  table.powm_cpu(result, exp);
}

void powm_precomputed_cpu_batch(fmpz_t *results, const fmpz_t *exps,
                                const PowmPrecomputeTable &table,
                                int batch_size) {
  if (batch_size >= 4) {
#pragma omp parallel for schedule(static)
    for (int i = 0; i < batch_size; i++) {
      table.powm_cpu(results[i], exps[i]);
    }
  } else {
    for (int i = 0; i < batch_size; i++) {
      table.powm_cpu(results[i], exps[i]);
    }
  }
}

void powm_precomputed_cpu_batch(fmpz_t *results, const fmpz_t *exps,
                                const PowmPrecomputeTable **tables,
                                int batch_size) {
  if (batch_size >= 4) {
#pragma omp parallel for
    for (int i = 0; i < batch_size; i++) {
      tables[i]->powm_cpu(results[i], exps[i]);
    }
  } else {
    for (int i = 0; i < batch_size; i++) {
      tables[i]->powm_cpu(results[i], exps[i]);
    }
  }
}

void powm_precomputed_cpu_batch(fmpz_t *results, const fmpz_t **exps,
                                const PowmPrecomputeTable **tables,
                                int batch_size) {
  if (batch_size >= 4) {
#pragma omp parallel for
    for (int i = 0; i < batch_size; i++) {
      tables[i]->powm_cpu(results[i], *(exps[i]));
    }
  } else {
    for (int i = 0; i < batch_size; i++) {
      tables[i]->powm_cpu(results[i], *(exps[i]));
    }
  }
}

// Generate a cryptographically secure random fmpz with exactly nbits using
// OpenSSL, filling FLINT limbs directly (no endian conversions needed)
void secure_random_fmpz(fmpz_t r, size_t nbits) {
  if (nbits == 0) {
    fmpz_zero(r);
    return;
  }
  const size_t limb_bits = (size_t)FLINT_BITS;
  size_t nlimbs = (nbits + limb_bits - 1) / limb_bits;
  ulong *limbs = (ulong *)malloc(nlimbs * sizeof(ulong));
  if (!limbs) {
    fprintf(stderr, "secure_random_fmpz: malloc failed\n");
    exit(1);
  }
  // Fill limbs with CSPRNG bytes
  if (RAND_bytes((unsigned char *)limbs, (int)(nlimbs * sizeof(ulong))) != 1) {
    fprintf(stderr, "secure_random_fmpz: RAND_bytes failed\n");
    free(limbs);
    exit(1);
  }
  // Mask high bits in the most significant limb and set the top bit to ensure
  // exact nbits
  size_t top_bits = nbits % limb_bits;
  if (top_bits == 0) {
    // nbits is a multiple of limb_bits; ensure highest bit of top limb is set
    limbs[nlimbs - 1] |= (ulong)1UL << (limb_bits - 1);
  } else {
    if (top_bits < limb_bits) {
      if (top_bits == sizeof(ulong) * 8) {
        // nothing to mask
      } else {
        limbs[nlimbs - 1] &= ((ulong)1UL << top_bits) - 1UL;
      }
    }
    limbs[nlimbs - 1] |= (ulong)1UL << (top_bits - 1);
  }
  fmpz_set_ui_array(r, limbs, (slong)nlimbs);
  free(limbs);
}