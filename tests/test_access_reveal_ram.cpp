// #include <algorithm>

// #include "access_reveal_ram.hpp"
// #include "test_util.hpp"
// using namespace PicoGRAM;
// struct Main : Gadget {
//   uint addr_width;
//   uint word_width;
//   AccessRevealRAM sam;
//   std::vector<uint64_t> addrs;
//   std::vector<uint64_t> ref;

//   DEFINE_FUNC(main, {}, [&](FuncInput) {
//     uint64_t address_val = addrs[get_time()];
//     uint64_t op_val = rand() % 2;
//     uint64_t data_val = rand() % (1UL << word_width);
//     Word address = Word::input_dbg(this, addr_width, address_val);
//     Word op = Word::input_dbg(this, 1, op_val);
//     Word data = Word::input_dbg(this, word_width, data_val);
//     FuncOutput result = sam.access({address, op, data});
//     uint64_t read_val = result[0].to_int();

//     if (get_mode() == EVAL || get_mode() == DEBUG) {
//       Assert_eq(read_val, ref[address_val]);
//       if (op_val == 0) {
//         ref[address_val] = data_val;
//       }
//     }
//     return FuncOutput();
//   });

//   Main(uint64_t T, Mode mode, uint addr_width, uint word_width,
//        uint log_way = 1)
//       : Gadget(mode, T),
//         addr_width(addr_width),
//         word_width(word_width),
//         sam(this, NONE, T, addr_width, word_width, log_way) {
//     set_name("Main");
//     for (uint64_t i = 0; i < T; ++i) {
//       addrs.push_back(i % (1UL << addr_width));
//     }
//     // shuffle the addresses
//     std::random_shuffle(addrs.begin(), addrs.end());
//     ref.resize(1UL << addr_width, 0);
//   }
// };

// TEST(AccessRevealRAM, AccessTest) {
//   uint64_t T_min = 256;
//   uint64_t T_max = 256;
//   uint addr_width = 5;
//   uint word_width = 5;
//   uint log_way = 2;
//   test_gadget<Main>(POW2, {T_min, T_max}, 1024 * 1024 * 32, addr_width,
//                     word_width, log_way);
// }

// TEST(AccessRevealRAM, AccessDbg) {
//   uint64_t T_min = 256;
//   uint64_t T_max = 1024;
//   uint addr_width = 7;
//   uint word_width = 6;
//   uint log_way = 2;
//   test_gadget_dbg<Main>(POW2, {T_min, T_max}, addr_width, word_width,
//   log_way);
// }