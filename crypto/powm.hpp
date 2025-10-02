#pragma once
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mod.h>
#include <inttypes.h>
#include <math.h>
#include <omp.h>
#include <stdlib.h>
// openssl
#include <openssl/rand.h>

#include <functional>
// #undef IPCL_AVAILABLE
#ifdef IPCL_AVAILABLE
#include <ipcl/ipcl.hpp>
#include <ipcl/mod_exp.hpp>
#endif

struct PowmPrecomputeTable {
  fmpz_t *table = nullptr;  // precompute table in Montgomery form
  fmpz_mod_ctx_t mont_ctx;  // Montgomery context for modular arithmetic
  size_t table_size;
  int log2b = 0;  // log2(radix)
  int b = 0;      // radix
  int L = 0;      // number of digits

  PowmPrecomputeTable() : table(nullptr), table_size(0) {}

  void cleanup() {
    if (table) {
      for (size_t i = 0; i < table_size; i++) {
        fmpz_clear(table[i]);
      }
      delete[] table;
      table = nullptr;
      table_size = 0;
    }
    if (table)  // Only clear if we initialized the context
      fmpz_mod_ctx_clear(mont_ctx);
  }

  ~PowmPrecomputeTable() { cleanup(); }

  void init_with_window_size(const fmpz_t base, const fmpz_t mod,
                             uint64_t exp_bits, int window_size) {
    cleanup();  // clear previous table if any
    // set mod and initialize Montgomery context
    fmpz_mod_ctx_init(mont_ctx, mod);

    // check log2b <= 16
    if (window_size < 1 || window_size > 20) {
      fprintf(stderr, "window_size must be in [1, 20]\n");
      exit(1);
    }
    this->log2b = window_size;
    b = 1 << log2b;

    L = (exp_bits + log2b - 1) / log2b;
    if (L < 1) L = 1;

    // Safety check on total entries
    size_t total_entries = (size_t)(b - 1) * (size_t)L + 1;
    if (total_entries > (size_t)1e8) {
      fprintf(stderr,
              "Refusing absurdly large table (%lu entries). Reduce "
              "precompute_mb.\n",
              (size_t)total_entries);
      exit(1);
    }
    table_size = total_entries;
    // allocate table of fmpz_t entries (index j*b + i), init all
    table = new fmpz_t[total_entries];
    if (!table) {
      fprintf(stderr, "malloc table failed\n");
      exit(1);
    }
    for (size_t idx = 0; idx < total_entries; ++idx) fmpz_init(table[idx]);

    // Convert base to Montgomery form
    fmpz_t base_mont, pow_tmp;
    fmpz_init(base_mont);
    fmpz_init(pow_tmp);
    fmpz_mod_set_fmpz(base_mont, base, mont_ctx);

    // Precompute table in Montgomery form:
    // table[j*b + i - 1] = g^{ i * b^j } mod mod (in Montgomery form)
    fmpz_mod_set_fmpz(table[0], base, mont_ctx);  // g^(b^0) in Montgomery form

    for (int j = 0; j < L; ++j) {
      unsigned long idx_base = (unsigned long)j * (unsigned long)(b - 1);

      // i = 2..b
      for (int i = 1; i < b; ++i) {
        // table[idx_base + i] = table[idx_base + i-1] * base_mont mod mod
        // (Montgomery multiplication)
        fmpz_mod_mul(table[idx_base + i], table[idx_base + i - 1],
                     table[idx_base], mont_ctx);
      }
    }
    fmpz_clear(base_mont);
    fmpz_clear(pow_tmp);
  }

  void init(const fmpz_t base, const fmpz_t mod, uint64_t exp_bits,
            size_t precompute_mb) {
    // Estimate bytes per entry (approx)
    size_t bytes_per_entry = (size_t)((fmpz_bits(mod) + 7) / 8) + 64;
    size_t allowed_bytes = precompute_mb << 20;
    if (allowed_bytes < bytes_per_entry) {
      fprintf(stderr,
              "Requested precompute size too small for even one entry (need >= "
              "%zu bytes)\n",
              bytes_per_entry);
      exit(1);
    }

    // Choose radix b in [2 .. BMAX]
    const int LOGBMAX = 20;

    for (log2b = 1; log2b <= LOGBMAX; log2b++) {
      b = 1 << log2b;
      size_t Ld = (exp_bits + log2b - 1) / log2b;
      if (Ld <= 0) Ld = 1;
      size_t entries = (size_t)(b - 1) * Ld + 1;
      long double needed = (long double)entries * (long double)bytes_per_entry;
      if (needed > (long double)allowed_bytes) {
        log2b--;
        b = 1 << log2b;
        break;
      }
    }
    if (log2b < 1) log2b = 1;

    init_with_window_size(base, mod, exp_bits, log2b);
  }

  void powm_cpu(fmpz_t result, const fmpz_t exp) const {
    if (!table) {
      fprintf(stderr, "powm called before init\n");
      fmpz_set_ui(result, 0);
      return;
    }

    if (fmpz_is_zero(exp)) {
      fmpz_set_ui(result, 1);
      return;
    }

    // Initialize result to 1 in normal form
    fmpz_one(result);

    // Use FLINT's built-in bit extraction for better performance
    slong exp_bits = fmpz_bits(exp);
    int iterations = (exp_bits + log2b - 1) / log2b;

    // check iterations not exceed L
    if (iterations > L) {
      fprintf(stderr,
              "Error: exponent has more bits (%ld) than precomputed (%d * %d = "
              "%d bits). Ignoring high bits.\n",
              exp_bits, L, log2b, L * log2b);
      // perform the regular powm
      exit(1);
    }

    for (int i = 0; i < iterations; i++) {
      // Extract next log2b bits using FLINT's bit operations
      int bits_this_iter = log2b;
      if ((i + 1) * log2b > exp_bits) bits_this_iter = exp_bits - i * log2b;

      int digit = 0;
      for (int j = 0; j < bits_this_iter; j++) {
        if (fmpz_tstbit(exp, i * log2b + j)) digit |= (1 << j);
      }

      if (digit != 0) {
        // Multiply result by precomputed table entry
        // Since table[idx] is in Montgomery form and mont_result is in normal
        // form, fmpz_mod_mul will handle the conversion automatically
        unsigned long idx = (unsigned long)i * (unsigned long)(b - 1) +
                            (unsigned long)(digit - 1);
        fmpz_mod_mul(result, result, table[idx], mont_ctx);
      }
    }

    // Result is already in normal form
  }
};

void powm_cpu(fmpz_t result, const fmpz_t base, const fmpz_t exp,
              const fmpz_t mod);

void powm_cpu_batch(fmpz_t *results, const fmpz_t *bases, const fmpz_t exp,
                    const fmpz_t mod, int batch_size);
void powm_cpu_batch(fmpz_t *results, const fmpz_t base, const fmpz_t *exps,
                    const fmpz_t mod, int batch_size);

void powm_cpu_batch(fmpz_t *results, const fmpz_t *base, const fmpz_t *exps,
                    const fmpz_t mod, int batch_size);
void powm_cpu_batch(fmpz_t *results, const fmpz_t **base, const fmpz_t **exps,
                    const fmpz_t mod, int batch_size);

void powm_cpu_batch(fmpz_t *results, const fmpz_t **base, const fmpz_t **exps,
                    const fmpz_t mod, int batch_size);

void powm_precomputed_cpu(fmpz_t result, const fmpz_t exp,
                          const PowmPrecomputeTable &table);

void powm_precomputed_cpu_batch(fmpz_t *results, const fmpz_t *exps,
                                const PowmPrecomputeTable &table,
                                int batch_size);

void powm_precomputed_cpu_batch(fmpz_t *results, const fmpz_t *exps,
                                const PowmPrecomputeTable **tables,
                                int batch_size);

void powm_precomputed_cpu_batch(fmpz_t *results, const fmpz_t **exps,
                                const PowmPrecomputeTable **tables,
                                int batch_size);
// Generate a cryptographically secure random fmpz with exactly nbits using
// OpenSSL, filling FLINT limbs directly (no endian conversions needed)
void secure_random_fmpz(fmpz_t r, size_t nbits);