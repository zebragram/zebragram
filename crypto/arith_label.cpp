#include "arith_label.hpp"
namespace ZebraGRAM {
void secret_share(ArithLabel &share1, ArithLabel &share2, const fmpz_t secret,
                  size_t max_bit, const PaillierPrivKey &sk) {
  if (max_bit > LEN_ARITH_DIGIT) {
    fprintf(stderr, "Error: max_bit exceeds LEN_ARITH_DIGIT in secret_share\n");
    exit(1);
  }
  const uint auth_share_len = max_bit + (LEN_K + LEN_STAT_PARAM);
  const uint raw_share_len = max_bit + LEN_STAT_PARAM;
  // sample random auth1 and raw1
  secure_random_fmpz(share1.auth, auth_share_len);
  secure_random_fmpz(share1.raw, raw_share_len);

  // compute share2.raw = secret + share1.raw
  fmpz_add(share2.raw, secret, share1.raw);
  // compute share2.auth = secret * K + share1.auth
  fmpz_mul(share2.auth, secret, sk.K);
  fmpz_add(share2.auth, share2.auth, share1.auth);
}

void secret_share(ArithLabel &share1, ArithLabel &share2, const fmpz_t secret,
                  const PaillierPrivKey &sk) {
  secret_share(share1, share2, secret, LEN_ARITH_DIGIT, sk);
}

void open_share(fmpz_t result, const ArithLabel &share1,
                const ArithLabel &share2) {
  // compute result = share2.raw - share1.raw
  fmpz_sub(result, share2.raw, share1.raw);
}
}  // namespace ZebraGRAM