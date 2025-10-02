#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#define MAX_GADGET_TIME 100
#define TOTAL_TIME 16384

constexpr double LN2 = 0.693147180559945309417232121458176568; // ln(2)

// epsilon_exp is the base-2 exponent: e.g. -60 means failure prob = 2^-60
double chernoff_upper_bound(uint64_t N, uint64_t M, uint64_t K,
                                   int epsilon_exp) {
    if (K == 0 || N == 0 || M == 0) return 0.0;

    // Number of "good" items (< N/K)
    uint64_t S = N / K; // floor
    double mu = static_cast<double>(M) * static_cast<double>(S) / static_cast<double>(N);

    // log(1/epsilon) = -epsilon_exp * ln(2)
    double log_term = -static_cast<double>(epsilon_exp) * LN2;

    // Chernoff bound: mu + sqrt(3 * mu * log(1/epsilon))
    double bound = ceil(mu + std::sqrt(3.0 * mu * log_term));

    // Can't exceed min(M, S)
    return std::min(bound, static_cast<double>(std::min(M, S)));
}

int main() {
    for (uint64_t T = 1; T <= TOTAL_TIME; T *= 2) {
        std::cout << "T: " << T << std::endl;
        const double mu = (double)T * MAX_GADGET_TIME / TOTAL_TIME;
        const double delta = sqrt(30.0 / mu);
        const uint64_t chernoff_bound = ceil(mu * (1 + delta));
        const uint64_t threshold_time = std::min((uint64_t)MAX_GADGET_TIME, chernoff_bound);
        std::cout << "threshold_time: " << threshold_time << std::endl;
        const uint64_t threshold_time_2 = chernoff_upper_bound(TOTAL_TIME, MAX_GADGET_TIME, (double)TOTAL_TIME / T, -20);
        std::cout << "threshold_time_2: " << threshold_time_2 << std::endl;
    }
    return 0;
}
