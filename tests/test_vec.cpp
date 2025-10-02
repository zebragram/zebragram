#include "gadget.hpp"
#include "test_util.hpp"
#include "vec.hpp"
using namespace ZebraGRAM;
struct TestVecMain : Gadget {
  Vec vec;
  std::vector<uint64_t> ref;
  uint data_width;
  std::function<FuncOutput(FuncInput)> main_func = [&](FuncInput) {
    uint64_t index = rand() % vec.size();
    Word index_word = Word::input_dbg(this, log2ceil(vec.size()), index);
    if (get_time() % 2) {
      uint64_t new_data = rand() % (1UL << data_width);
      const Word& old_data =
          vec.write(index_word, Word::input_dbg(this, data_width, new_data));
      uint64_t old_data_val = old_data.to_int();
      if (get_mode() == EVAL || get_mode() == DEBUG) {
        EXPECT_EQ(old_data_val, ref[index]);
        ref[index] = new_data;
      }
    } else {
      uint64_t data = vec.read(index_word).to_int();
      if (get_mode() == EVAL || get_mode() == DEBUG) {
        EXPECT_EQ(data, ref[index]);
      }
    }
    return FuncOutput();
  };
  DEFINE_FUNC(main, {}, main_func);
  TestVecMain(uint64_t T, Mode mode, uint size, uint data_width)
      : Gadget(mode, T),
        vec(this, size),
        ref(size, 0),
        data_width(data_width) {}
};

TEST(Vec, RandTest) {
  test_gadget<TestVecMain>(ANY, {100, 200}, 1024 * 1024 * 128, rand() % 20 + 1,
                           rand() % 6 + 1);
}

TEST(Vec, RandDbgTest) {
  test_gadget_dbg<TestVecMain>(ANY, {100, 200}, rand() % 30 + 1,
                               rand() % 10 + 1);
}