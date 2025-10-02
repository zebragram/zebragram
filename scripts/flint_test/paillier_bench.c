#include <openssl/bn.h>
#include <openssl/rand.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <inttypes.h>
#include <omp.h>

// Generate a random BIGNUM with the specified number of bits
BIGNUM *rand_bn_bits(int bits) {
    BIGNUM *bn = BN_new();
    if (!bn) {
        fprintf(stderr, "BN_new failed\n");
        exit(1);
    }
    if (!BN_rand(bn, bits, 1, 0)) {  // top=1 ensures full bit length
        fprintf(stderr, "BN_rand failed\n");
        exit(1);
    }
    return bn;
}

void benchmark_eval(int n_bits, int exp_bits) {
    BN_CTX *ctx = BN_CTX_new();
    if (!ctx) {
        fprintf(stderr, "BN_CTX_new failed\n");
        exit(1);
    }

    // Generate random N and compute N^2
    BIGNUM *N = rand_bn_bits(n_bits);
    BIGNUM *N2 = BN_new();
    BN_mul(N2, N, N, ctx);

    // Random base g in [0, N^2)
    BIGNUM *g = rand_bn_bits(n_bits * 2 - 1);
    BN_mod(g, g, N2, ctx);

    // Random exponent
    BIGNUM *exp = rand_bn_bits(exp_bits);

    // Allocate result
    BIGNUM *res = BN_new();

    // Benchmark loop
    int iterations = 100;
    double start = omp_get_wtime();
    for (int i = 0; i < iterations; i++) {
        BN_mod_exp(res, g, exp, N2, ctx);
    }
    double end = omp_get_wtime();

    double total_time = end - start;
    printf("\nOpenSSL benchmark results:\n");
    printf("Performed %d modular exponentiations mod N^2 (N=%d bits, exp=%d bits)\n",
           iterations, n_bits, exp_bits);
    printf("Total time: %.4f seconds\n", total_time);
    printf("Average time per operation: %.4f seconds\n", total_time / iterations);

    // Cleanup
    BN_free(N); BN_free(N2); BN_free(g); BN_free(exp); BN_free(res);
    BN_CTX_free(ctx);
}

#include <flint/fmpz.h>
#include <flint/flint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// Generate a random fmpz with specified bits
void rand_fmpz_bits(fmpz_t x, int bits, flint_rand_t state) {
    fmpz_randbits(x, state, bits);
    // Ensure top bit is set to get full bit-length
    fmpz_setbit(x, bits - 1);
}

void benchmark_flint_eval(int n_bits, int exp_bits, int batch_size, int num_threads) {
    flint_rand_t state;
    flint_randinit(state);

    fmpz_t N, N2, g, exp;
    fmpz_init(N); fmpz_init(N2); fmpz_init(g); fmpz_init(exp);
    // Random N and compute N^2
    rand_fmpz_bits(N, n_bits, state);
    fmpz_mul(N2, N, N);

    // Random base g in [0, N^2)
    rand_fmpz_bits(g, n_bits * 2 - 1, state);
    fmpz_mod(g, g, N2);

    // Random exponent
    rand_fmpz_bits(exp, exp_bits, state);

    int iterations = 100;
    int adjusted_iterations = iterations;
    double start = omp_get_wtime();
    if (num_threads == 1) {
        fmpz_t res;
        fmpz_init(res);
        for (int i = 0; i < iterations * batch_size; i++) {
            fmpz_powm(res, g, exp, N2);
        }
        fmpz_clear(res);
    } else {
        fmpz_t res[num_threads];
        for (int t = 0; t < num_threads; t++) {
            fmpz_init(res[t]);
        }
        adjusted_iterations = (iterations * num_threads + batch_size - 1) / batch_size;
        for (int i = 0; i < adjusted_iterations; i++) {
            #pragma omp parallel for num_threads(num_threads) schedule(static)
            for (int j = 0; j < batch_size; j++) {
                fmpz_powm(res[omp_get_thread_num()], g, exp, N2);
            }
        }
        for (int t = 0; t < num_threads; t++) {
            fmpz_clear(res[t]);
        }
    }
    double end = omp_get_wtime();
    printf("\nFLINT benchmark results:\n");
    double total_time = end - start;
    printf("Performed %d iterations of modular exponentiations mod N^2 (N=%d bits, exp=%d bits)\n",
           adjusted_iterations, n_bits, exp_bits);
    printf("Total time: %.4f seconds\n", total_time);
    printf("Average time per operation: %.4f seconds\n", total_time / adjusted_iterations);

    fmpz_clear(N); fmpz_clear(N2); fmpz_clear(g); fmpz_clear(exp); 
    flint_randclear(state);
}

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <string.h>

#include <openssl/bn.h>
#include <openssl/rand.h>

static double now_sec(void) {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + t.tv_nsec / 1e9;
}


/*
  benchmark_garbler:
    base_bits    - number of bits for N (N ~ base_bits)
    exp_bits     - exponent size in bits
    precompute_mb- requested precompute memory size in megabytes (approx)
*/
void benchmark_garbler(int base_bits, int exp_bits, double precompute_mb) {
    BN_CTX *ctx = BN_CTX_new();
    if (!ctx) { fprintf(stderr, "BN_CTX_new failed\n"); return; }

    // 1) Build modulus N and N^2 (use random N with top bit set)
    BIGNUM *N = BN_new();
    BN_rand(N, base_bits, BN_RAND_TOP_ONE, BN_RAND_BOTTOM_ODD);
    BIGNUM *N2 = BN_new();
    if (!BN_mul(N2, N, N, ctx)) { fprintf(stderr, "BN_mul(N,N) failed\n"); goto cleanup0; }

    // 2) random base g in [1, N2-1]
    BIGNUM *g = BN_new();
    if (!BN_rand_range(g, N2)) { fprintf(stderr, "BN_rand_range failed\n"); goto cleanup0; }
    if (BN_is_zero(g)) BN_one(g);

    // 3) random exponent of exp_bits
    BIGNUM *exp = BN_new();
    BN_rand(exp, exp_bits, BN_RAND_TOP_ANY, BN_RAND_BOTTOM_ANY);

    // 4) Prepare montgomery context for N2
    BN_MONT_CTX *mont = BN_MONT_CTX_new();
    if (!mont || !BN_MONT_CTX_set(mont, N2, ctx)) {
        fprintf(stderr, "Montgomery context setup failed\n"); goto cleanup0;
    }

    // Convert g to montgomery form
    BIGNUM *g_mont = BN_new();
    if (!BN_to_montgomery(g_mont, g, mont, ctx)) { fprintf(stderr, "BN_to_montgomery failed\n"); goto cleanup0; }

    // Estimate size per precomputed entry (approx): use BN_num_bytes(N2) and add small overhead
    size_t bytes_per_entry = (size_t)BN_num_bytes(N2) + 64; // +64 bytes as overhead estimate
    size_t allowed_bytes = (size_t)(precompute_mb * 1024 * 1024ULL);
    if (allowed_bytes < bytes_per_entry) {
        fprintf(stderr, "Requested precompute size too small for even one entry (need >= %zu bytes)\n", bytes_per_entry);
        goto cleanup0;
    }

    // Choose radix b automatically: find largest b in [2 .. BMAX] such that
    // (b-1) * L(b) * bytes_per_entry <= allowed_bytes
    // where L(b) = ceil(exp_bits / log2(b))
    const int BMAX = 1 << 16; // cap radix to 65536 (safe: BN_div uses BIGNUM divisor if needed)
    int best_b = 2;
    for (int b = 2; b <= BMAX; ++b) {
        double log2b = log2((double)b);
        double Ld = ceil((double)exp_bits / log2b);
        if (Ld <= 0) Ld = 1;
        uint64_t entries = (uint64_t)(b - 1) * (uint64_t)Ld;
        // check overflow and memory
        long double needed = (long double)entries * (long double)bytes_per_entry;
        if (needed <= (long double)allowed_bytes) best_b = b;
        else break; // as b increases, entries usually increase strongly -> break early
    }
    int b = best_b;
    double log2b = log2((double)b);
    int L = (int)ceil((double)exp_bits / log2b);
    if (L < 1) L = 1;

    printf("Selected radix b = %d, digits L = %d (approx precompute memory %lf MB, bytes/entry=%zu)\n",
           b, L, precompute_mb, bytes_per_entry);

    // Safety check: compute total entries and check allocation size
    uint64_t total_entries = (uint64_t)(b) * (uint64_t)L; // we will index j*b + i for i in [0..b-1]
    if (total_entries > (uint64_t)1e8) {
        fprintf(stderr, "Refusing absurdly large table (%llu entries). Reduce precompute_mb.\n", (unsigned long long)total_entries);
        goto cleanup0;
    }

    // allocate pointer table of size L * b (we'll use index j*b + i), initialize to NULL
    BIGNUM **table = calloc((size_t)total_entries, sizeof(BIGNUM*));
    if (!table) { fprintf(stderr, "calloc table failed\n"); goto cleanup0; }

    // Precompute table:
    // table[j*b + i] holds g^{ i * b^j } (in MONT form), for i = 1..b-1, j = 0..L-1
    // We'll compute base_bj = g^{b^j} progressively:
    BIGNUM *b_bn = BN_new();
    BN_set_word(b_bn, (unsigned long)b);

    double t_pre_start = now_sec();
    BIGNUM *prev_base = NULL; // prev_base = g^{b^{j-1}} in mont form
    for (int j = 0; j < L; ++j) {
        BIGNUM *base_bj = BN_new();
        if (!base_bj) { fprintf(stderr, "BN_new failed\n"); goto precompute_fail; }

        if (j == 0) {
            // g^{b^0} = g
            if (!BN_copy(base_bj, g_mont)) { fprintf(stderr, "BN_copy failed\n"); BN_free(base_bj); goto precompute_fail; }
        } else {
            // base_bj = prev_base ^ b  (prev_base already in mont domain)
            // note: BN_mod_exp_mont_consttime expects base in normal domain; prev_base is mont form.
            // So convert prev_base back to normal, exponentiate, then convert to mont.
            BIGNUM *prev_normal = BN_new();
            if (!BN_from_montgomery(prev_normal, prev_base, mont, ctx)) { fprintf(stderr, "from_mont failed\n"); BN_free(prev_normal); BN_free(base_bj); goto precompute_fail; }

            if (!BN_mod_exp_mont_consttime(base_bj, prev_normal, b_bn, N2, ctx, mont)) {
                fprintf(stderr, "BN_mod_exp_mont_consttime(base^b) failed\n"); BN_free(prev_normal); BN_free(base_bj); goto precompute_fail;
            }
            BN_free(prev_normal);
        }

        // Now base_bj is in mont form and equals g^{b^j} (mont)
        // Fill i = 1..b-1: table[j*b + i] = base_bj^i = cumulative multiply
        unsigned long idx_base = (unsigned long)j * (unsigned long)b;
        // i=0 unused (we'll leave it NULL)
        // i=1
        table[idx_base + 1] = BN_new();
        if (!BN_copy(table[idx_base + 1], base_bj)) { fprintf(stderr, "BN_copy failed\n"); goto precompute_fail; }

        for (int i = 2; i <= b - 1; ++i) {
            table[idx_base + i] = BN_new();
            if (!table[idx_base + i]) { fprintf(stderr, "BN_new failed table fill\n"); goto precompute_fail; }
            if (!BN_mod_mul_montgomery(table[idx_base + i],
                                       table[idx_base + i - 1],
                                       base_bj, mont, ctx)) {
                fprintf(stderr, "BN_mod_mul_montgomery failed at j=%d i=%d\n", j, i);
                goto precompute_fail;
            }
        }

        // keep base_bj as prev_base for next iteration
        if (prev_base) BN_free(prev_base);
        prev_base = base_bj; // will be freed later
    }

    double t_pre_end = now_sec();
    printf("Precomputation finished in %.6f s. Table entries stored: %llu\n",
           t_pre_end - t_pre_start, (unsigned long long)total_entries);

    // Decompose exponent into base-b digits once (least-significant-first)
    unsigned int *digits = calloc((size_t)L, sizeof(unsigned int));
    if (!digits) { fprintf(stderr, "digits calloc failed\n"); goto precompute_fail; }

    {
        BIGNUM *e_copy = BN_dup(exp);
        BIGNUM *quo = BN_new();
        BIGNUM *rem = BN_new();
        for (int j = 0; j < L; ++j) {
            if (BN_is_zero(e_copy)) {
                digits[j] = 0;
                // keep e_copy = 0 for remaining digits
            } else {
                if (!BN_div(quo, rem, e_copy, b_bn, ctx)) {
                    fprintf(stderr, "BN_div failed in decomposition\n"); BN_free(quo); BN_free(rem); BN_free(e_copy); goto precompute_fail;
                }
                unsigned long dv = BN_get_word(rem); // rem < b so fits in unsigned long
                digits[j] = (unsigned int)dv;
                BN_free(e_copy);
                e_copy = BN_dup(quo);
            }
        }
        BN_free(quo); BN_free(rem); BN_free(e_copy);
    }

    // Prepare accumulator one in mont form
    BIGNUM *one = BN_new();
    BN_one(one);
    BIGNUM *one_mont = BN_new();
    BN_to_montgomery(one_mont, one, mont, ctx);

    // Benchmark: multiply table lookups according to digits (mont multiplies)
    const int iterations = 100; // reasonable default
    BIGNUM *acc = BN_new();
    double t0 = now_sec();
    for (int it = 0; it < iterations; ++it) {
        if (!BN_copy(acc, one_mont)) { fprintf(stderr, "BN_copy acc failed\n"); goto bench_fail; }
        for (int j = 0; j < L; ++j) {
            unsigned int d = digits[j];
            if (d == 0) continue;
            BIGNUM *tbl = table[(size_t)j * b + d];
            if (!tbl) { fprintf(stderr, "Missing table[%d,%u]\n", j, d); goto bench_fail; }
            if (!BN_mod_mul_montgomery(acc, acc, tbl, mont, ctx)) { fprintf(stderr, "mont mul failed runtime\n"); goto bench_fail; }
        }
        // (Optionally convert result out of mont form)
        // BIGNUM *res_normal = BN_new();
        // BN_from_montgomery(res_normal, acc, mont, ctx);
        // BN_free(res_normal);
    }
    double t1 = now_sec();
    double avg_precomp_method = (t1 - t0) / iterations;

    // Baseline: OpenSSL BN_mod_exp_mont_consttime
    BIGNUM *res_baseline = BN_new();
    double t2 = now_sec();
    for (int it = 0; it < iterations; ++it) {
        if (!BN_mod_exp_mont_consttime(res_baseline, g, exp, N2, ctx, mont)) { fprintf(stderr, "baseline exp failed\n"); goto bench_fail; }
    }
    double t3 = now_sec();
    double avg_baseline = (t3 - t2) / iterations;

    printf("\nBenchmark results (averaged over %d iterations):\n", iterations);
    printf("  Precomputed-table method: %.8f s/op\n", avg_precomp_method);
    printf("  OpenSSL BN_mod_exp_mont_consttime baseline: %.8f s/op\n", avg_baseline);

    // Cleanup and free table
bench_fail:
    free(digits);
    // free table elements
    for (uint64_t idx = 0; idx < total_entries; ++idx) {
        if (table[idx]) BN_free(table[idx]);
    }
    free(table);
    if (prev_base) BN_free(prev_base);
    BN_free(b_bn);
    BN_free(one); BN_free(one_mont);
    BN_free(acc); BN_free(res_baseline);
    BN_MONT_CTX_free(mont);

cleanup0:
    BN_free(g_mont);
    BN_free(g);
    BN_free(exp);
    BN_free(N2);
    BN_free(N);
    BN_CTX_free(ctx);
    return;

precompute_fail:
    fprintf(stderr, "Precompute failed; cleaning up\n");
    // free partial table
    if (table) {
        for (uint64_t idx = 0; idx < total_entries; ++idx) if (table[idx]) BN_free(table[idx]);
        free(table);
    }
    if (prev_base) BN_free(prev_base);
    BN_free(b_bn);
    BN_free(one); BN_free(one_mont);
    BN_MONT_CTX_free(mont);
    BN_free(g_mont);
    BN_free(g);
    BN_free(exp);
    BN_free(N2);
    BN_free(N);
    BN_CTX_free(ctx);
}

void benchmark_flint_garbler(int base_bits, int exp_bits, double precompute_mb) {
    flint_rand_t state;
    flint_randinit(state);

    // 1) Build modulus N and N^2
    fmpz_t N, N2, g, exp;
    fmpz_init(N); fmpz_init(N2); fmpz_init(g); fmpz_init(exp);

    rand_fmpz_bits(N, base_bits, state);
    fmpz_mul(N2, N, N);            // N2 = N * N

    // 2) random base g in [1, N2-1]
    fmpz_randm(g, state, N2);
    if (fmpz_is_zero(g)) fmpz_one(g);

    // 3) random exponent
    rand_fmpz_bits(exp, exp_bits, state);

    // Estimate bytes per entry (approx)
    size_t bytes_per_entry = (size_t)((fmpz_bits(N2) + 7) / 8) + 64;
    size_t allowed_bytes = (size_t)(precompute_mb * 1024 * 1024.0);
    if (allowed_bytes < bytes_per_entry) {
        fprintf(stderr, "Requested precompute size too small for even one entry (need >= %zu bytes)\n",
                bytes_per_entry);
        goto cleanup_all;
    }

    // Choose radix b in [2 .. BMAX]
    const int BMAX = 1 << 16;
    int best_b = 2;
    for (int b = 2; b <= BMAX; ++b) {
        double log2b = log2((double)b);
        double Ld = ceil((double)exp_bits / log2b);
        if (Ld <= 0) Ld = 1;
        unsigned long long entries = (unsigned long long)(b - 1) * (unsigned long long)Ld;
        long double needed = (long double)entries * (long double)bytes_per_entry;
        if (needed <= (long double)allowed_bytes) best_b = b;
        else break;
    }
    int b = best_b;
    double log2b = log2((double)b);
    int L = (int)ceil((double)exp_bits / log2b);
    if (L < 1) L = 1;

    // approximate memory (used entries = (b-1)*L)
    long double approx_bytes = (long double)( (unsigned long long)(b - 1) * (unsigned long long)L )
                                * (long double)bytes_per_entry;
    double approx_mb = (double)(approx_bytes / (1024.0L * 1024.0L));

    printf("Selected radix b = %d, digits L = %d (approx precompute memory %.6f MB, bytes/entry=%zu)\n",
           b, L, approx_mb, bytes_per_entry);

    // Safety check on total entries
    unsigned long long total_entries = (unsigned long long)b * (unsigned long long)L;
    if (total_entries > (unsigned long long)1e8) {
        fprintf(stderr, "Refusing absurdly large table (%llu entries). Reduce precompute_mb.\n",
                (unsigned long long)total_entries);
        goto cleanup_all;
    }

    // allocate table of fmpz_t entries (index j*b + i), init all
    fmpz_t *table = malloc(sizeof(fmpz_t) * (size_t)total_entries);
    if (!table) { fprintf(stderr, "malloc table failed\n"); goto cleanup_all; }
    for (unsigned long long idx = 0; idx < total_entries; ++idx) fmpz_init(table[idx]);

    // Precompute:
    // table[j*b + i] = g^{ i * b^j } mod N2, for i=1..b-1, j=0..L-1
    double t_pre_start = now_sec();

    fmpz_t base_bj, pow_tmp, b_f;
    fmpz_init(base_bj);
    fmpz_init(pow_tmp);
    fmpz_init(b_f);
    fmpz_set_ui(b_f, (unsigned long)b);

    fmpz_set(base_bj, g); // g^(b^0)

    for (int j = 0; j < L; ++j) {
        unsigned long idx_base = (unsigned long)j * (unsigned long)b;

        // i = 1
        fmpz_set(table[idx_base + 1], base_bj);
        // i = 2..b-1
        for (int i = 2; i <= b - 1; ++i) {
            // table[idx_base + i] = table[idx_base + i-1] * base_bj mod N2
            fmpz_mul(pow_tmp, table[idx_base + i - 1], base_bj);
            fmpz_mod(table[idx_base + i], pow_tmp, N2);
        }

        // compute base_bj := base_bj^b mod N2 for next j (avoid huge intermediates)
        if (j < L - 1) {
            fmpz_powm(pow_tmp, base_bj, b_f, N2);
            fmpz_set(base_bj, pow_tmp);
        }
    }

    fmpz_clear(base_bj);
    fmpz_clear(pow_tmp);

    double t_pre_end = now_sec();
    printf("Precomputation finished in %.6f s. Table entries stored (approx) = %llu\n",
           t_pre_end - t_pre_start, (unsigned long long)((unsigned long long)(b - 1) * (unsigned long long)L));

    // Decompose exponent into base-b digits (LSB-first)
    unsigned int *digits = calloc((size_t)L, sizeof(unsigned int));
    if (!digits) { fprintf(stderr, "digits calloc failed\n"); goto precompute_fail; }

    fmpz_t tmp_q, tmp_r, tmp_e;
    fmpz_init(tmp_q); fmpz_init(tmp_r); fmpz_init(tmp_e);
    fmpz_set(tmp_e, exp);

    for (int j = 0; j < L; ++j) {
        if (fmpz_is_zero(tmp_e)) {
            digits[j] = 0;
        } else {
            // tmp_q = tmp_e / b; tmp_r = tmp_e mod b
            fmpz_fdiv_qr(tmp_q, tmp_r, tmp_e, b_f);
            digits[j] = (unsigned int) fmpz_get_ui(tmp_r); // rem < b so fits
            fmpz_set(tmp_e, tmp_q);
        }
    }
    fmpz_clear(tmp_q); fmpz_clear(tmp_r); fmpz_clear(tmp_e);

    // Prepare accumulator and run benchmark loops
    const int iterations = 100;
    fmpz_t acc;
    fmpz_init(acc);

    double t0 = now_sec();
    for (int it = 0; it < iterations; ++it) {
        fmpz_one(acc);
        for (int j = 0; j < L; ++j) {
            unsigned int d = digits[j];
            if (d == 0) continue;
            // multiply acc = acc * table[j*b + d] mod N2
            fmpz_mul(acc, acc, table[(size_t)j * b + d]);
            fmpz_mod(acc, acc, N2);
        }
    }
    double t1 = now_sec();
    double avg_precomp_method = (t1 - t0) / iterations;

    // Baseline: fmpz_powm (modular exponentiation)
    fmpz_t res_baseline;
    fmpz_init(res_baseline);
    double t2 = now_sec();
    for (int it = 0; it < iterations; ++it) {
        fmpz_powm(res_baseline, g, exp, N2);
    }
    double t3 = now_sec();
    double avg_baseline = (t3 - t2) / iterations;

    // Single-run correctness check: compute one precomp result and compare with fmpz_powm
    fmpz_t acc_check;
    fmpz_init(acc_check);
    fmpz_one(acc_check);
    for (int j = 0; j < L; ++j) {
        unsigned int d = digits[j];
        if (d == 0) continue;
        fmpz_mul(acc_check, acc_check, table[(size_t)j * b + d]);
        fmpz_mod(acc_check, acc_check, N2);
    }
    // compute baseline once
    fmpz_t base_once;
    fmpz_init(base_once);
    fmpz_powm(base_once, g, exp, N2);

    if (fmpz_cmp(acc_check, base_once) == 0) {
        printf("Correctness: precomputed result matches fmpz_powm.\n");
    } else {
        printf("ERROR: mismatch between precomputed result and fmpz_powm!\n");
    }

    printf("\nBenchmark results (averaged over %d iterations):\n", iterations);
    printf("  Precomputed-table method: %.8f s/op\n", avg_precomp_method);
    printf("  FLINT fmpz_powm baseline: %.8f s/op\n", avg_baseline);

    // cleanup single-run temps
    fmpz_clear(acc_check);
    fmpz_clear(base_once);

precompute_fail:
    free(digits);
    for (unsigned long long idx = 0; idx < total_entries; ++idx) {
        fmpz_clear(table[idx]);
    }
    free(table);
    fmpz_clear(b_f);
    fmpz_clear(acc);
    fmpz_clear(res_baseline);

cleanup_all:
    fmpz_clear(N); fmpz_clear(N2); fmpz_clear(g); fmpz_clear(exp);
    flint_randclear(state);
}

int main() {
    // benchmark_garbler(1536, 256, 100);
    printf("\n ==== OpenSSL garbler per half exponentiation === \n");
    benchmark_garbler(3108, 256, 128); // g^r and (g^k)^r mod p^2 * q^2 using CRT

    printf("\n ==== FLINT garbler per half exponentiation === \n");
    benchmark_flint_garbler(3108, 256, 128);

    // evaluator time
    // benchmark_eval(3072, (3072 - 256) / 2 - 60 + 256);
    // use 3108 bit so that 4kB block ops can be efficiently parallelized with 48 threads
    printf("\n ==== OpenSSL evaluator per exponentiation === \n");
    benchmark_eval(3108, (3108 - 256) / 2 + 256);
    // garbler time

    printf("\n ==== FLINT evaluator per exponentiation === \n");
    benchmark_flint_eval(3108, (3108 - 256) / 2 + 256, 1, 1);

    int num_cores = omp_get_max_threads();
    int batch_size = 192;
    printf("\n ==== FLINT evaluator batch %d exponentiation with %d threads === \n", batch_size, num_cores);
    benchmark_flint_eval(3108, (3108 - 256) / 2 + 256, batch_size, num_cores);
    return 0;
}
