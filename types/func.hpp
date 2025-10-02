#pragma once
#include <functional>
#include <vector>

#include "arith_word.hpp"
#include "simd_word.hpp"

namespace ZebraGRAM {
using SIMDFuncInput = const std::vector<SIMDWord>&;
using SIMDFuncOutput = std::vector<SIMDWord>;

// using FuncInput = const std::vector<Word>&;
template <typename Word_Type>
struct FuncParamVec : public std::vector<Word_Type> {
  // ArithWord payload;

  FuncParamVec() = default;

  FuncParamVec(std::initializer_list<Word_Type> init)
      : std::vector<Word_Type>(init) {}

  // constructor from std::vector<Word>
  FuncParamVec(const std::vector<Word_Type>& vec)
      : std::vector<Word_Type>(vec) {}

  // FuncParamVec(const std::vector<Word_Type>& vec, const ArithWord& payload) :
  // std::vector<Word_Type>(vec), payload(payload) {}

  // get the total bit width of the input words
  uint total_bit_width() const {
    return std::accumulate(
        this->begin(), this->end(), 0UL,
        [](uint acc, const Word_Type& w) { return acc + w.width(); });
  }

  uint total_arith_width() const {
    return std::accumulate(
        this->begin(), this->end(), 0UL, [](uint acc, const Word_Type& w) {
          return acc + (w.has_payload() ? w.get_payload().width() : 0);
        });
  }

  uint64_t total_gc_size() const {
    return total_bit_width() * sizeof(Label) +
           total_arith_width() * ArithLabel::byte_length;
  }
};

using FuncInput = const FuncParamVec<Word>&;

using FuncOutput = FuncParamVec<Word>;

#define DEFINE_SIMDFUNC(name, out_widths_param, func) \
  SIMDFunc name = SIMDFunc(self, func, out_widths_param, STR(name));

#define DEFINE_FUNC(name, out_widths_param, func) \
  Func name = Func(self, func, out_widths_param, STR(name));

#define DEFINE_FUNC_ARITH(name, out_widths_param, out_arith_widths_param, \
                          func)                                           \
  Func name =                                                             \
      Func(self, func, out_widths_param, out_arith_widths_param, STR(name));

struct Link;
struct SIMDLink;

/**
 * @brief The base class for function types
 *
 */
struct FuncType : DataType {
 protected:
  // the output SIMDWord widths of the SIMD function
  std::vector<uint> out_widths;
  std::vector<uint> out_arith_widths;

  // the name of the function, must be unique in the same gadget
  std::string name;

 public:
  FuncType(Gadget* owner, const std::vector<uint>& out_widths,
           const std::string& name)
      : DataType(owner), out_widths(out_widths), name(name) {
    Assert(owner);
  }

  FuncType(Gadget* owner, const std::vector<uint>& out_widths,
           const std::vector<uint>& out_arith_widths, const std::string& name)
      : DataType(owner),
        out_widths(out_widths),
        out_arith_widths(out_arith_widths),
        name(name) {
    Assert(owner);
  }

  // delete copy constructor
  FuncType(const FuncType&) = delete;
  // delete assignment operator
  FuncType& operator=(const FuncType&) = delete;

  // move constructor
  FuncType(FuncType&&) = default;
  // move assignment operator
  FuncType& operator=(FuncType&&) = default;

  /**
   * @brief Garble the function for one time step
   *
   */
  virtual void garble() = 0;

  /**
   * @brief Get the name of the function
   *
   * @return const std::string&
   */
  const std::string& get_name() const { return name; }

  /**
   * @brief Get the caller gadget
   *
   * @return Gadget*
   */
  Gadget* get_caller() const;

  /**
   * @brief Get the type of link between the function's callee and
   * caller
   *
   * @return LinkType
   */
  LinkType get_link_type() const;

  const std::vector<uint>& get_out_widths() const { return out_widths; }

  const std::vector<uint>& get_out_arith_widths() const {
    return out_arith_widths;
  }

  uint sum_out_width() const {
    return std::accumulate(out_widths.begin(), out_widths.end(), 0UL);
  }
};

/**
 * @brief A wrapper for functions that take SIMD words as inputs and outputs.
 * The SIMDFunc may be called multiple times during eachtimestep with
 * different input widths, but the output width should always be fixed.
 *
 */
struct SIMDFunc : FuncType {
 private:
  // the width of the input SIMD words for each call in atimestep
  std::vector<std::vector<uint>> in_widths;
  // the public flags of the input SIMD words for each call in atimestep
  std::vector<std::vector<std::vector<bool>>> in_pub_es;

  // the begin bit offsets of the input SIMD words for each call in atimestep
  std::vector<std::vector<uint>> in_bit_offsets;

  std::vector<std::vector<uint>> in_shifts;

  // the latesttimestep of the caller, used in operator() to reset
  // called_count
  uint64_t caller_last_time = -1;

  // the latest timestep of the callee, used in garble() to reset called_count
  uint64_t callee_last_time = -1;

  // the number of times the function is called during the current timestep
  uint64_t called_count = 0;

  // the internal function to be called
  std::function<SIMDFuncOutput(SIMDFuncInput)> func;

  /**
   * @brief Get the simd link for calling the function
   *
   * @return SIMDLink*
   */
  SIMDLink* get_simd_link() const;

 public:
  /**
   * @brief Construct a new SIMDFunc
   *
   * @param owner the owner gadget
   * @param func the internal function to be called
   * @param out_widths the output SIMDWord widths of the SIMD function
   * @param name the name of the function
   */
  SIMDFunc(Gadget* owner = NULL,
           const std::function<SIMDFuncOutput(SIMDFuncInput)>& func = NULL,
           const std::vector<uint>& out_widths = {},
           const std::string& name = "");

  // delete copy constructor
  SIMDFunc(const SIMDFunc&) = delete;
  // delete assignment operator
  SIMDFunc& operator=(const SIMDFunc&) = delete;

  // move constructor
  SIMDFunc(SIMDFunc&&) = default;
  // move assignment operator
  SIMDFunc& operator=(SIMDFunc&&) = default;

  /**
   * @brief A conditional function call.
   * For the evaluator, if the call is real, translate the inputs to
   * callee's timezone, make the function call, and translate the outputs back.
   * Otherwise, just update the SIMD link to prepare for the next call.
   * For the garbler, only update the SIMD link and defer the real function
   * call.
   *
   * @param control the control bit. If the control bit is 0, the call is real,
   * otherwise it is fake.
   * @param inputs the SIMD inputs
   * @return std::vector<SIMDWord> the SIMD outputs, the output bit offsets must
   * follow immediately after the input bit offsets. If the call is fake, the
   * outputs will be not contain any encoding. The pub_e flags in the output
   * words should always be false.
   */
  std::vector<SIMDWord> operator()(const Bit& control, SIMDFuncInput inputs);

  /**
   * @brief Alternative interface of operator()
   *
   * @param control the control bit. If the control bit is 0, the call is real,
   * otherwise it is fake.
   * @param inputs
   * @param outputs intended to be modified to output the execution result
   */
  void exec(const Bit& control, SIMDFuncInput inputs, SIMDFuncOutput& outputs);

  /**
   * @brief Garble the function for one time step. Mock all function calls
   * performed by this->func, and add these functions to the garble list of
   * the callee gadget so that they can be garbled later.
   *
   */
  void garble() override;
};

/**
 * @brief  A wrapper for functions that take words as inputs and outputs.
 * Similar to the SIMDFunc, the same Func may be called multiple times during
 * eachtimestep with different input widths, but the output width should
 * always be fixed.
 *
 */
struct Func : FuncType {
 private:
  // the original function, used when translation is not needed (main function
  // or when child has the same T as parent).
  std::function<FuncOutput(FuncInput)> original_func;

  // When the link is SIMD, the function delegates all the operations to a
  // SIMDFunc. We first convert the Words into SIMDWords under the caller's
  // timezone, then call simd_func, which translates the SIMDWord into the
  // callee's timezone. Finally, the internal function of the simd_func converts
  // the SIMDWord back to a vector of output Words.
  SIMDFunc* simd_func = NULL;

  // when the link is non-simd cond, used by the garbler to record the sizes of
  // the input words (the input words are cached in the link)
  std::vector<uint> input_sizes;

  // when the link is direct, used by the garbler to record the input/output
  // words after mocking the function call
  std::vector<std::vector<Word>> in_words;
  std::vector<std::vector<Word>> out_words;

  // only store the meta data for measurement to save memory
  std::vector<std::vector<WordMetaData>> in_words_meta_data;
  std::vector<std::vector<WordMetaData>> out_words_meta_data;

  // for the garbler to keep track of the current input/output words offset when
  // garble the function
  uint64_t garble_counter = 0;

 public:
  Func(Gadget* owner = NULL,
       const std::function<FuncOutput(FuncInput)>& func = NULL,
       const std::vector<uint>& out_widths = {}, const std::string& name = "");

  Func(Gadget* owner = NULL,
       const std::function<FuncOutput(FuncInput)>& func = NULL,
       const std::vector<uint>& out_widths = {},
       const std::vector<uint>& out_arith_widths = {},
       const std::string& name = "");

  /**
   * @brief Function call only for the main function, which has no input/output.
   *
   */
  void operator()();

  /**
   * @brief A conditional function call.
   * For the evaluator, if the call is real, translate the inputs to
   * callee's timezone, make the function call, and translate the outputs back.
   * Otherwise, just update the link to prepare for the next call.
   * For the garbler, only update the link and defer the real function
   * call.
   *
   * @param control the control bit. If the control bit is 0, the call is real,
   * otherwise it is fake.
   * @param inputs the input words
   * @return FuncOutput the output words. If the call is fake, the
   * outputs will return a vector of constant zeros. The pub_e flags in the
   * output words should always be false.
   */
  FuncOutput operator()(const Bit& control, FuncInput inputs);

  /**
   * @brief Alternative interface of operator(). Allows the output languages to
   * be predefined for the direct link and non-simd conditional link.
   *
   * @param control the control bit. If the control bit is 0, the call is real,
   * otherwise it is fake.
   * @param inputs the input words
   * @param outputs when the mode is GARBLE and the link is direct or non-simd
   * conditional, outputs should contain the predefined languages for the output
   * words. This way, we can let the output of different children in a call
   * group share the same language.
   *
   */
  void exec(const Bit& control, FuncInput inputs, FuncOutput& outputs);

  /**
   * @brief A wrapper of conditional function call where the condition is
   * always true.
   *
   * @param inputs
   * @return FuncOutput
   */
  FuncOutput operator()(FuncInput inputs);

  /**
   * @brief Garble the function for one timestep. Mock all function calls
   * performed by this->func, and add these functions to the garble list of the
   * callee gadget so that they can be garbled later.
   *
   */
  void garble() override;

  // copy constructor
  Func(const Func& other) = delete;

  // move constructor
  Func(Func&& other) = default;

  ~Func() {
    if (simd_func) {
      delete simd_func;
    }
  }
};
}  // namespace ZebraGRAM