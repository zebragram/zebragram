#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include "integer_permutation.hpp"  // Include the relevant libsnark files
#include "test_util.hpp"
#include "waksman.hpp"

std::vector<size_t> generate_random_permutation(size_t size) {
  std::vector<size_t> permutation(size);
  for (size_t i = 0; i < size; ++i) {
    permutation[i] = i;
  }

  // Shuffle the permutation
  PicoGRAM::secure_permute(permutation.begin(), permutation.end());

  return permutation;
}

TEST(PermuteVector, Waksman) {
  for (size_t size = 2; size <= 128; ++size) {
    // Create a vector of integers from 0 to size - 1
    std::vector<int> original_vector(size);
    for (size_t i = 0; i < size; ++i) {
      original_vector[i] = i;
    }

    // Generate a random permutation
    std::vector<size_t> random_permutation = generate_random_permutation(size);

    // Apply the AS-Waksman permutation
    std::vector<int> permuted_vector = PicoGRAM::waksman_permute_vector(
        original_vector, random_permutation, [](bool cond, int &a, int &b) {
          if (cond) std::swap(a, b);
        });

    // Verify the result
    for (size_t i = 0; i < size; ++i) {
      // The element at random_permutation[i] in the original vector
      // should now be at index i in the permuted vector
      assert(permuted_vector[random_permutation[i]] == original_vector[i]);
    }
  }
}
