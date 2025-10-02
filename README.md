# gramimpl
An implementation of ZebraGRAM - a garbled RAM with communication cost of O(w * log N + lambda * log^3 N). This project is build upon the implementation of PicoGRAM from Crypto' 25. On a high level, we extend PicoGRAM with the arithmetic circuit garbling techniques, achieving an order of magnitude better communication cost when the data block is large.

## Install Dependencies
First update the googletest sub repository
```
git submodule update --init --recursive
```
You can use our docker image to install the rest of the dependencies
```
docker build -t zebragram-builder .
```
On Linux
```
docker run -it --rm -v ${pwd}/zebragram -w /zebragram zebragram-builder
```
On Windows
```
docker run -it --rm -v ${PWD}:/zebragram -w /zebragram zebragram-builder
```

## Build in Debug Mode
```
rm -r build
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug
ninja -C build
```

## Build in Release Mode
```
rm -r build
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
ninja -C build
```

## Test Arithmetic GRAM
```
./build/test_basic --gtest_filter=*RandPayloadMedium*
```

## Bench 256MB Database with 4kB blocks (position map not included)
(Warning: OOM error may occur for machine with insufficient memory)
```
./build/test_basic --gtest_filter=*Payload256M4k*
```
## Bench 256MB Database with 1kB blocks (position map not included)
(Warning: OOM error may occur for machine with insufficient memory)
```
./build/test_basic --gtest_filter=*Payload256M1k*
```

## Bench Bit Decomposition
```
./build/test_basic --gtest_filter=*DecomposeAll*
```

## Evaluation of Garbled Boolean/Arithmetic Circuit
Given a boolean/arithmetic circuit, a garbler can garble each gate in the circuit and send the garbled circuit (GC) to the evaluator. The evaluator evaluates every gate of the circuit in the same order as the garbler.

In this implementation, we define basic data types such as `Bit` and `Word` and overload common operators in `types/` to facilitate the garbling and evaluation of boolean circuits. When operating on these data types, the garbler writes to the GC sequentially, and later when the evaluator performs the same operation, it will read from the GC also sequentially.

## GRAM and Dynamic Language Translation Problem
The key difficulty in constructing a garbled RAM is that the evaluation order cannot be determined statically. Specifically, we want a function $f_a$ to call another function $f_b$ conditionally based on the runtime value of a control wire. If the condition is unmet, no gate in $f_b$ should be evaluated.
ZebraGRAM extends the idea of Tri-state Circuit by Heath et al. and constructed an improved garbled stack to address this challenge.

## Abstractions
### Basic Data Types
A Bit represents the garbling of a Boolean wire in the circuit. We implement XOR (`^`) and AND (`&`) using the techniques of Free XOR and half gates. 

A Word consists of multiple bits and can be interpreted as an integer. We define bitwise operations on words as well as arithmetic operations such as addition and comparison. For each Word, we additionally let it hold a payload of type `ArithWord`, which stores the labels of arithmetic wires with `fmpz_t` type from the FLINT library.

A SIMDWord represents the garbling of a cable in ZebraGRAM. The garbler's share is a single finite field element (`BigInt`), and the evaluator's share is an array of group elements (`ECPoint`). A SIMDWord can be translated from/to a Word.

### Gadget (`gadgets/gadget.hpp`)
A gadget is a boolean sub-circuit whose evaluation order is deterministic. As an example, each node in circuit oram tree always runs the same set of read and evict functions at every timestep. Every data object holds a pointer to a gadget that owns it. Each gadget owns a contigueous segment in the GC of the whole circuit, and can be garbled and evaluated like a normal boolean circuit. A gadget may hold states, just like a C++ struct. 

To initialize a gadget, one needs to specify `T`, the maximum number of timesteps it can be repeatedly executed. We also require gadgets to be mutually connected into a tree structure, representing the call graph. When initializing a non-root gadget, one also needs to specify its parent gadget.

### Function(`types/func.hpp`)
We allow a gadget to have multiple `Func` or `SIMDFunc`. 

A `Func` takes in a `vector` of `Word`s and returns a `vector` of `Word`s. When declaring a function, one needs to specify the number and width of each return word.

As an example, the following negate function outputs the negation of each input word, and we declare the function returns two words of width 3 and 2 respectively. The plus function adds two inputs, and outputs a word of width 4.
```cpp
DEFINE_FUNC(negate, std::vector<uint>({3, 2}), [&](FuncInput inputs) {
    FuncOutput outputs;
    for (const Word& word : inputs) {
      outputs.emplace_back(~word);
    }
    return outputs;
  });

DEFINE_FUNC(plus, {4}, [&](FuncInput inputs) {
    FuncOutput outputs;
    Word sum = inputs[0] + inputs[1];
    outputs.emplace_back(sum);
    return outputs;
});
```
A function may also access the state variables in the gadget.

A `SIMDFunc` is used either in the internal implementation of `Func` and or can be explicitily defined to avoid unnecessary translations between `Word` and `SIMDWord`. The interface is similar to `Func` except that the `Word`s are replaced with `SIMDWord`s in the input and output.

### Link (`gadgets/link.hpp`)
Links are used to connect a gadget with its parent (caller) gadget. A `Link` or `SIMDLink` implements the stack and reversed stack in ZebraGRAM, which takes inspiration from the oblivious compaction algorithm by Goodrich. In particular, a `SIMDLink` utilizes a novel garbling scheme based on the DDH assumption to reduce the communication cost, but may increase computational overhead due to EC operations.

At each timestamp, the parent gadget may call the `update_link` method to conditionally connect itself with the child gadget based on a control bit. Then, the evaluator can perform any number of `translate` operations to convert between a `Word`/`SIMDWord` owned by the parent gadget and a `Word`/`SIMDWord` owned by the child gadget.

A direct link is used to connect two gadgets with synchronous timestamps, i.e. the calling schedule must be fixed.

### Group (`gadgets/gadget_group.hpp`)
When designing a tree-based data structure, it's common that a parent node calls one of its children at each timestep. A `Group` defines multiple children nodes of the same type and allows the function of the parent gadget to call a function of one of the children node. Below, we give an example where the parent has two children each with an incrementer circuit that can be called `T/2` times.
```cpp
struct ParentStateful : Gadget {
  const uint data_width = 8;
  Group<ChildStateful> childs;
  uint64_t left_count = 0; // keep track of how many times the left child has been called
  
  std::function<FuncOutput(FuncInput)> func = [&](FuncInput) {
    uint64_t selector_value = (rand() % (T - get_time())) < (T / 2 - left_count);
    left_count += selector_value;
    Word selector = Word::input_dbg(this, 1, selector_value);
    uint64_t data = rand() % (1UL << data_width);
    Word query = Word::input_dbg(this, data_width, data);
    FuncOutput out = childs[selector].call("inc", {query});
    return FuncOutput();
  };
  DEFINE_FUNC(main, {}, func);
  ParentStateful(uint64_t T, Mode mode)
      : Gadget(T, mode), childs(T / 2, this, 2) {}
};
```
When T is small, non-simd links are instantiated to reduce computational overhead; for large T, simd links are applied as they can significantly reduce communication cost. The threshold can be customized.

### Circuit ORAM (`circuit_oram.hpp` and `rec_oram.hpp`)
Using the primitives described above, we implemented a non-recursive circuit ORAM in `circuit_oram.hpp` and a recursive circuit ORAM in `rec_oram.hpp`. The non-recursive circuit ORAM is faster but requires the client to manage the positions of the blocks in the ORAM tree. The recursive circuit ORAM manages the positions of the blocks for the client using a recursive position map. Both ORAMs allow clients to update a block with a custom function. We also provide a simpler example of non-oblivious RAM in `access_reveal_ram.hpp`, which fully leverages the SIMD acceleration.

## Order of Garbling and Evaluation
Within each gadget, the garbler and the evaluator execute the circuit in the same order. However, in terms of the whole circuit, their execution order are different.
### Garbler's Execution Order
The garbler starts with the `main` function root gadget, and always garble the circuitry of all the timestamps for each gadget before moving on to the next gadget. This design choice ensures good locality when writing to the GC, as well as lower memory consumption. Crucially, the garbler mocks the output of each function call and defer the real garbling process of the function. 

### Evaluator's Execution Order
The evaluator also starts with the main function and only execute a function call if the link to the child gadget is active. Namely, the execution order may depend on the runtime input. To ensure privacy, the algorithm must be oblivious so that the execution order does not reveal the secret inputs.

### EMP-Tool Interface (`ram.h`, `emp/`)
The garbled RAM can be used as an extension of the emp-sh2pc to achieve semi-honest two-party computation. We define a minimal ORAM interface in `emp/emp_interface.hpp`, and a wrapper in `ram.h`. We utilize the IO channel implementation of the emp-tool for communication between the garbler and the evaluator. 

The garbler first garbles the ORAM trees and send the material to the evaluator, and the evaluator caches it (in memory or on disk). Then the parties run the rest of computational tasks in a pipelined fashion just like ordinary emp-tool execution. We need this extra preprocessing phase because the execution order of the evaluator is not deterministic and should be hidden from the garbler.

## Evaluate costs of basic cryptographic operations
```
./build/test_basic --gtest_filter=*TestCrypto*
```

## Scripts for simulating baselines and plotting graphs
```
cd ./scripts
python plot.py
```