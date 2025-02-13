// #include "func.hpp"
// #include "gadget.hpp"
// #include "gadget_group.hpp"
// #include "test_util.hpp"
// using namespace PicoGRAM;
// struct Child : Gadget {
//   DEFINE_FUNC(negate, std::vector<uint>({3, 2}), [&](FuncInput inputs) {
//     FuncOutput outputs;
//     for (const Word& word : inputs) {
//       const Word& neg_word = ~word;
//       outputs.emplace_back(neg_word);
//     }
//     return outputs;
//   });
//   DEFINE_FUNC(plus, {4}, [&](FuncInput inputs) {
//     FuncOutput outputs;
//     Word sum = inputs[0] + inputs[1];
//     outputs.emplace_back(sum);
//     return outputs;
//   });
//   Child(Gadget* main, LinkType link_type, uint64_t T)
//       : Gadget(main, link_type, T) {
//     set_name("Child");
//   }
// };

// struct Parent : Gadget {
//   Child child;
//   uint64_t real_count = 0;
//   uint64_t T;
//   // SIMDFunc main;
//   DEFINE_FUNC(main, {}, [&](FuncInput) {
//     uint64_t remain_real_count = T / 2 - real_count;
//     uint8_t is_real = (rand() % (T - get_time())) < remain_real_count;
//     real_count += is_real;
//     Bit is_real_bit = Bit::input_dbg(this, !is_real);
//     const uint data0_width = 3;
//     const uint data1_width = 2;
//     uint64_t data0 = 0;
//     uint64_t data1 = 0;
//     uint64_t neg_data0 = (~data0) & ((1UL << data0_width) - 1);
//     uint64_t neg_data1 = (~data1) & ((1UL << data1_width) - 1);
//     uint64_t ref_sum = neg_data0 + neg_data1;

//     Word query0 = Word::input_dbg(this, data0_width, data0);
//     Word query1 = Word::input_dbg(this, data1_width, data1);
//     FuncOutput out = child.negate(is_real_bit, {query0, query1});
//     FuncOutput sum = child.plus(is_real_bit, out);
//     uint64_t sum_value = sum[0].to_int();
//     Assert(!out[0].is_pub_e());
//     Assert(!out[1].is_pub_e());
//     if (get_mode() == EVAL && is_real) {
//       EXPECT_EQ(sum_value, ref_sum);
//     }
//     out = child.negate(is_real_bit, {~query0, query1});
//     sum = child.plus(is_real_bit, {out[1], out[0]});
//     sum_value = sum[0].to_int();
//     if (get_mode() == EVAL && is_real) {
//       EXPECT_EQ(sum_value, data0 + neg_data1);
//     }
//     return FuncOutput();
//   });
//   Parent(uint64_t T, Mode mode, LinkType link_type = SIMD_COND)
//       : Gadget(mode, T), child(this, link_type, T / 2), T(T) {
//     set_name("Parent");
//   }
// };

// // struct ParentSameT : Gadget {
// //   Child child;
// //   // SIMDFunc main;
// //   DEFINE_FUNC(main, {}, [&](FuncInput) {
// //     Bit is_real_bit = Bit::input_dbg(this, 0);
// //     const uint data0_width = 3;
// //     const uint data1_width = 2;
// //     uint64_t data0 = rand() % (1UL << data0_width);
// //     uint64_t data1 = rand() % (1UL << data1_width);
// //     uint64_t neg_data0 = (~data0) & ((1UL << data0_width) - 1);
// //     uint64_t neg_data1 = (~data1) & ((1UL << data1_width) - 1);
// //     uint64_t ref_sum = neg_data0 + neg_data1;

// //     Word query0 = Word::input_dbg(this, data0_width, data0);
// //     Word query1 = Word::input_dbg(this, data1_width, data1);
// //     FuncOutput out = child.negate(is_real_bit, {query0, query1});
// //     FuncOutput sum = child.plus(is_real_bit, out);
// //     uint64_t sum_value = sum[0].to_int();
// //     if (get_mode() == EVAL) {
// //       EXPECT_EQ(sum_value, ref_sum);
// //     }
// //     out = child.negate(is_real_bit, {~query0, query1});
// //     sum = child.plus(is_real_bit, {out[1], out[0]});
// //     sum_value = sum[0].to_int();
// //     if (get_mode() == EVAL) {
// //       EXPECT_EQ(sum_value, data0 + neg_data1);
// //     }
// //     return FuncOutput();
// //   });
// //   ParentSameT(uint64_t T, Mode mode) : Gadget(T, mode), child(T, this) {
// //     set_name("Parent");
// //   }
// // };

// struct DuoChildParent : Gadget {
//   static constexpr uint64_t num_children = 2;
//   Group<Child> childs;
//   uint64_t left_count = 0;
//   uint64_t T;
//   std::function<FuncOutput(FuncInput)> func = [&](FuncInput) {
//     uint64_t remain_real_count = T / 2 - left_count;
//     uint64_t selector_value = (rand() % (T - get_time())) <
//     remain_real_count; left_count += selector_value; Word selector =
//         Word::input_dbg(this, log2ceil(num_children), selector_value);
//     const uint data0_width = 3;
//     const uint data1_width = 2;
//     uint64_t data0 = rand() % (1UL << data0_width);
//     uint64_t data1 = rand() % (1UL << data1_width);
//     uint64_t neg_data0 = (~data0) & ((1UL << data0_width) - 1);
//     uint64_t neg_data1 = (~data1) & ((1UL << data1_width) - 1);
//     uint64_t ref_sum = neg_data0 + neg_data1;

//     Word query0 = Word::input_dbg(this, data0_width, data0);
//     Word query1 = Word::input_dbg(this, data1_width, data1);
//     FuncOutput out = childs[selector].call("negate", {query0, query1});
//     uint64_t out0 = out[0].to_int();
//     uint64_t out1 = out[1].to_int();
//     if (get_mode() == EVAL) {
//       EXPECT_EQ(out0, neg_data0);
//       EXPECT_EQ(out1, neg_data1);
//     }
//     FuncOutput sum = childs[selector].call("plus", out);
//     uint64_t sum_value = sum[0].to_int();
//     if (get_mode() == EVAL) {
//       EXPECT_EQ(sum_value, ref_sum);
//     }
//     return FuncOutput();
//   };
//   // SIMDFunc main;
//   DEFINE_FUNC(main, {}, func);
//   DuoChildParent(uint64_t T, Mode mode)
//       : Gadget(mode, T),
//         childs(this, std::vector<uint64_t>(num_children, T / num_children)),
//         T(T) {
//     set_name("Duo Child Parent");
//   }
// };

// TEST(TestFunc, CallSIMD) { test_gadget<Parent>(ANY, {2, 100}); }

// TEST(TestFunc, CallNonSIMD) {
//   test_gadget<Parent>(ANY, {2, 100}, 16UL << 20, COND);
// }

// TEST(TestFunc, MeasureNonSIMD) { measure_gadget_gc<Parent>(2, COND); }

// // TEST(TestFunc, CallSameT) { test_gadget<ParentSameT>(ANY, {2, 100}); }

// TEST(TestFunc, CallGroup) { test_gadget<DuoChildParent>(EVEN, {2, 50}); }

// struct ChildStateful : Gadget {
//   Word sum = Word::constant(this, 8, (uint64_t)0);
//   const uint data_width = 8;
//   DEFINE_FUNC(inc, std::vector<uint>({data_width}), [&](FuncInput inputs) {
//     FuncOutput outputs;
//     sum = (sum + inputs[0]).slice(data_width);
//     outputs.emplace_back(sum);
//     return outputs;
//   });
//   ChildStateful(Gadget* caller, LinkType link_type, uint64_t T)
//       : Gadget(caller, link_type, T) {
//     set_name("Child");
//   }
// };

// struct ParentStateful : Gadget {
//   const uint data_width = 8;
//   static constexpr uint64_t num_children = 2;
//   Group<ChildStateful> childs;
//   uint64_t left_count = 0;
//   uint64_t ref_sums[num_children] = {0};
//   uint64_t T;
//   std::function<FuncOutput(FuncInput)> func = [&](FuncInput) {
//     uint64_t remain_real_count = T / 2 - left_count;
//     uint64_t selector_value = (rand() % (T - get_time())) <
//     remain_real_count; left_count += selector_value;
//     // std::cout << "Selector value: " << selector_value << std::endl;
//     Word selector =
//         Word::input_dbg(this, log2ceil(num_children), selector_value);

//     uint64_t data = rand() % (1UL << data_width);

//     Word query = Word::input_dbg(this, data_width, data);
//     FuncOutput out = childs[selector].call("inc", {query});
//     ref_sums[selector_value] += data;
//     ref_sums[selector_value] &= (1UL << data_width) - 1;
//     uint64_t out_sum = out[0].to_int();
//     if (get_mode() == EVAL) {
//       EXPECT_EQ(out_sum, ref_sums[selector_value]);
//     }
//     return FuncOutput();
//   };
//   DEFINE_FUNC(main, {}, func);
//   ParentStateful(uint64_t T, Mode mode)
//       : Gadget(mode, T),
//         childs(this, std::vector<uint64_t>(num_children, T / num_children)),
//         T(T) {
//     set_name("Parent");
//   }
// };

// TEST(TestFunc, CallStateful) { test_gadget<ParentStateful>(EVEN, {2, 50}); }

// struct Level0 : Gadget {
//   Word sum = Word::constant(this, 8, (uint64_t)0);
//   const uint data_width = 8;
//   DEFINE_FUNC(inc, std::vector<uint>({data_width}), [&](FuncInput inputs) {
//     FuncOutput outputs;
//     sum = (sum + inputs[0]).slice(data_width);
//     outputs.emplace_back(sum);
//     return outputs;
//   });
//   Level0(Gadget* caller, LinkType link_type, uint64_t T)
//       : Gadget(caller, link_type, T) {
//     set_name("Level 0");
//   }
// };

// struct Level1 : Gadget {
//   const uint data_width = 8;
//   static constexpr uint64_t num_children = 2;
//   Group<Level0> childs;
//   uint64_t left_count = 0;
//   uint64_t ref_sums[num_children] = {0};
//   uint64_t T;
//   std::function<FuncOutput(FuncInput)> func = [&](FuncInput inputs) {
//     uint64_t remain_real_count = T / 2 - left_count;

//     uint64_t selector_value = (rand() % (T - get_time())) <
//     remain_real_count; left_count += selector_value; Word selector =
//         Word::input_dbg(this, log2ceil(num_children), selector_value);

//     Word query = inputs[0];
//     FuncOutput out = childs[selector].call("inc", {query});
//     out.push_back(selector);
//     return out;
//   };
//   DEFINE_FUNC(inc, std::vector<uint>({data_width, 1}), func);
//   Level1(Gadget* caller, LinkType link_type, uint64_t T)
//       : Gadget(caller, link_type, T),
//         childs(this, std::vector<uint64_t>(num_children, T / num_children)),
//         T(T) {
//     set_name("Level 1");
//   }
// };

// struct Level2 : Gadget {
//   const uint data_width = 8;
//   static constexpr uint64_t num_children = 2;
//   Group<Level1> childs;
//   uint64_t left_count = 0;
//   uint64_t ref_sums[4] = {0};
//   uint64_t T;
//   std::function<FuncOutput(FuncInput)> func = [&](FuncInput) {
//     uint64_t remain_real_count = T / 2 - left_count;
//     uint64_t selector_value = (rand() % (T - get_time())) <
//     remain_real_count; left_count += selector_value; Word selector =
//         Word::input_dbg(this, log2ceil(num_children), selector_value);

//     uint64_t data = rand() % (1UL << data_width);

//     Word query = Word::input_dbg(this, data_width, data);
//     FuncOutput out = childs[selector].call("inc", {query});

//     Word child_selector = out[1];
//     uint64_t child_selector_value = child_selector.to_int();
//     uint64_t index = selector_value * 2 + child_selector_value;
//     ref_sums[index] += data;
//     ref_sums[index] &= (1UL << data_width) - 1;
//     uint64_t out_sum = out[0].to_int();
//     if (get_mode() == EVAL || get_mode() == DEBUG) {
//       EXPECT_EQ(out_sum, ref_sums[index]);
//     }
//     return FuncOutput();
//   };
//   DEFINE_FUNC(main, {}, func);
//   Level2(uint64_t T, Mode mode)
//       : Gadget(mode, T),
//         childs(this, std::vector<uint64_t>(num_children, T / num_children)),
//         T(T) {
//     set_name("Level 2");
//   }
// };

// TEST(TestFunc, CallTwoLevels) { test_gadget<Level2>(POW2, {4, 32}); }
// TEST(TestFunc, CallTwoLevelsDbg) { test_gadget_dbg<Level2>(POW2, {64, 128});
// }

// struct Level2MultiCalls : Gadget {
//   const uint data_width = 8;
//   Level1 child;
//   uint64_t ref_sums[2] = {0};
//   uint64_t access_time;
//   uint64_t T;
//   std::function<FuncOutput(FuncInput)> func = [&](FuncInput) {
//     for (uint64_t i = 0; i < access_time / T; ++i) {
//       uint64_t data = rand() % (1UL << data_width);

//       Word query = Word::input_dbg(this, data_width, data);
//       child.inc_time();
//       FuncOutput out = child.inc({query});

//       Word child_selector = out[1];
//       uint64_t child_selector_value = child_selector.to_int();
//       uint64_t index = child_selector_value;
//       ref_sums[index] += data;
//       ref_sums[index] &= (1UL << data_width) - 1;
//       uint64_t out_sum = out[0].to_int();
//       if (get_mode() == EVAL || get_mode() == DEBUG) {
//         EXPECT_EQ(out_sum, ref_sums[index]);
//       }
//     }
//     return FuncOutput();
//   };
//   DEFINE_FUNC(main, {}, func);
//   Level2MultiCalls(uint64_t T, Mode mode, uint64_t access_time)
//       : Gadget(mode, T),
//         child(this, DIRECT, access_time),
//         access_time(access_time),
//         T(T) {
//     set_name("Level 2");
//   }
// };

// TEST(TestFunc, CallTwoLevelsMultiCalls) {
//   test_gadget<Level2MultiCalls>(POW2, {2, 2}, 1024 * 1024 * 128, 2);
// }

// TEST(TestFunc, CallTwoLevelsMultiCallsDbg) {
//   test_gadget_dbg<Level2MultiCalls>(POW2, {2, 2}, 2);
// }
