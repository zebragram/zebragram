#include "rand.hpp"
namespace ZebraGRAM {
void secure_random_128bit(unsigned char* random_bytes) {
  // Generate 128 bits (16 bytes) of secure random data
  if (RAND_bytes(random_bytes, 16) != 1) {
    // If the function fails, throw an exception
    throw std::runtime_error("Error generating secure random number");
  }
}

void secure_random(uint8_t* random_bytes, size_t len) {
  if (RAND_bytes(random_bytes, len) != 1) {
    // If the function fails, throw an exception
    throw std::runtime_error("Error generating secure random number");
  }
}

uint64_t secure_random_uint64() {
  uint64_t random_uint64;
  if (RAND_bytes((unsigned char*)&random_uint64, 8) != 1) {
    // If the function fails, throw an exception
    throw std::runtime_error("Error generating secure random number");
  }
  return random_uint64;
}

uint8_t secure_random_byte() {
  uint8_t random_byte;
  if (RAND_bytes(&random_byte, 1) != 1) {
    // If the function fails, throw an exception
    throw std::runtime_error("Error generating secure random number");
  }
  return random_byte;
}

RandBitPool rand_bit_pool;

}  // namespace ZebraGRAM