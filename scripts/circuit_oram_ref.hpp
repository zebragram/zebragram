#pragma once
#include <sys/types.h>

#include <cstddef>
#include <cstdint>
#include <functional>
#include <vector>

#include "util.hpp"
namespace ZebraGRAM {
/**
 * @brief Data block used when initializing the ORAM
 *
 */
struct PlainTextBlock {
  uint64_t position;  // the position of the block
  uint64_t v_addr;    // the virtual address of the block
  uint64_t data;      // the data of the block

  PlainTextBlock(uint64_t position, uint64_t v_addr, uint64_t data)
      : position(position), v_addr(v_addr), data(data) {}
};

// a vector of blocks in the same position
using BlocksInPos = std::vector<PlainTextBlock>;

using BlocksInPosIter = std::vector<BlocksInPos>::iterator;

using Word = uint64_t;
using Bit = bool;

using FuncOutput = std::vector<Word>;
using FuncInput = std::vector<Word>;

Word mux(const Bit& sel, const Word& a, const Word& b) { return sel ? b : a; }

Bit mux(const Bit& sel, const Bit& a, const Bit& b) { return sel ? b : a; }

Bit get_bit(const Word& word, const uint& idx) { return (word >> idx) & 1; }

void set_bit(Word& word, const uint& idx, const Bit& bit) {
  word = (word & ~(1ULL << idx)) | ((uint64_t)bit << idx);
}

Word trailing_zeros(const Word& word, const uint& width, const uint& pack) {
  uint max_tz = (width + pack - 1) / pack;
  if (!word) return max_tz;
  // use intrinsic function to calculate trailing zeros
  return std::min((uint)__builtin_ctzll(word) / pack, max_tz);
}

// get the bit mask that only keeps the most significant one bit
Word ms_one(const Word& word) {
  if (word == 0) return 0;
  return 1UL << log2(word);
}

Word slice(const Word& word, const uint& start, const uint& end) {
  return (word >> start) & ((1ULL << (end - start)) - 1);
}

Word slice(const Word& word, const uint& end) { return slice(word, 0, end); }

void cond_swap(const Bit& control, Word& a, Word& b) {
  if (control) {
    std::swap(a, b);
  }
}

void cond_swap(const Bit& control, Bit& a, Bit& b) {
  if (control) {
    std::swap(a, b);
  }
}

/**
 * @brief A k-ary circuit ORAM tree (or sub-tree). The data structure
 * accepts either a read or evict request at each time step. To maximize
 * efficiency, eviction is also performed during a read request on the read
 * path. If the element to read is found in a bucket, no other element will
 * be evicted from the same bucket.
 *
 */
struct CircuitORAMRefTree {
  enum Op { READ, EVICT };
  uint v_addr_width;    // the width of the virtual address word
  uint word_width;      // the width of the data word
  uint level_width;     // the width of the total level of the tree
  uint position_width;  // the width of the position word
  uint level;           // the level in the tree (root has level 0)
  Word level_word;      // the level in the tree as a Word
  uint remain_level;    // the remaining level in the tree, i.e., total_level -
                        // level
  uint bkt_size;    // the capacity of the bucket (or stash if this is the root)
  uint log_way;     // the log of the number of ways at the current level
  Word bkt_valids;  // a bit vector indicating whether each slot in the bucket
                    // is valid (i.e. occupied)
  std::vector<Word> bkt_positions;  // the positions of the data in the bucket,
                                    // if level > 0, the repetitive bits in
                                    // position is trimmed
  std::vector<Word> bkt_v_addrs;    // the virtual addresses of the data in the
                                    // bucket
  std::vector<Word> bkt_datas;      // the data in the bucket

  CircuitORAMRefTree** children =
      NULL;  // the children (sub-tree) of the current node

  /**
   * @brief Find the data in the bucket that can be evicted to the deepest
   * level, or if it's a read request, the one that matches v_addr if exists. If
   * multiple addresses can go to the last level and none of them matches the
   * input virtual address, the one with the largest index will be chosen.
   * @param position the path to evict
   * @param is_read whether the operation is a read
   * @param v_addr the virtual address to match. If there is a match, treat it
   * as deepest.
   * @return FuncOutput[0] the max level an address can go to
   *        FuncOutput[1] a one-hot word indicating the location of the deepest
   * address FuncOutput[2] a bit indicating whether the deepest address
   * perfectly matches the input address
   */
  FuncOutput get_local_deepest(const Word& position, const Bit& is_read,
                               const Word& v_addr) {
    Word max_common_suffix_len = 0;
    Word is_deepest_one_hot = 0;
    // whether the deepest address perfectly matches
    // the input address
    Bit has_match = 0;

    for (uint i = 0; i < bkt_size; ++i) {
      Bit is_valid = get_bit(bkt_valids, i);
      const Word& position_xor = bkt_positions[i] ^ position;
      const Word& common_suffix_len =
          trailing_zeros(position_xor, position_width, log_way);
      // if there's a match before, we should not consider it as the deepest
      Bit is_deepest =
          is_valid & (common_suffix_len >= max_common_suffix_len) & !has_match;
      max_common_suffix_len =
          mux(is_deepest, max_common_suffix_len, common_suffix_len);
      // the implementation below is equivalent to !addr_xor but more efficient
      has_match |= is_valid & is_read & (!(bkt_v_addrs[i] ^ v_addr));
      set_bit(is_deepest_one_hot, i, is_deepest);
    }
    // remove 1 in one-hot encoding except the last one
    is_deepest_one_hot = ms_one(is_deepest_one_hot);
    const Word& deepest_level = level_word + max_common_suffix_len;
    Word has_match_word = has_match;
    return FuncOutput({deepest_level, is_deepest_one_hot, has_match_word});
  }

  // below are temporary variables used at each time step
  Word op;        // the operation to perform in the current time step
  Word position;  // the path to read / evict in the current time step
  Word v_addr;  // the virtual address to read from in the current time step, or
                // 0 if it's an evict operation
  Word local_deepest_one_hot;  // marks which location in bkt can be evicted
                               // deepest
  Bit local_has_match;  // whether there's a match for read in the current
                        // bucket
  Bit should_send;      // whether the current level should send an element to a
                        // level below
  Bit should_drop;      // whether the current level should receive an element
                        // from a level above
  Bit is_local_deeper;  // whether there's an element in the bucket deeper than
                        // the element from level above

  // a reference to the sub-tree corresponding to the position
  CircuitORAMRefTree* child_entry;

  // the following function sets the op, position, and v_addr in all nodes on
  // the path, and updates child_entry
  // these information can be sent through the path without re-encoded from word
  // into  words
  // therefore, we separate the function from meta_scan
  std::function<FuncOutput(FuncInput)> set_addr = [&](FuncInput inputs) {
    op = inputs[0];
    position = inputs[1];
    v_addr = inputs[2];

    if (remain_level) {
      const Word& child_idx = slice(position, log_way);
      child_entry = children[child_idx];
      child_entry->set_addr({op, position >> log_way, v_addr});
    }
    return FuncOutput{};
  };

  // meta scan to decide evict schedule
  // we perform several optimizations over the pseudo code in the tri-state
  // paper
  // 1. We perform eviction on read path almost for free. Compared to a full
  // eviction, we skip evicting elements from a bucket with the element to read,
  // since we want to scan each bucket only once.
  // 2. In the pseudo code, a level remembers a target level to send element to
  // here, we let each level remember whether it should receive an element.
  // Therefore, we don't need to pass the "dest" from children to parents, and
  // during the actual eviction, we don't need to pass this "dest" from parents
  // to children.
  // 3. We avoid some unnecessary operations for invalidating variables that
  // store level information
  std::function<FuncOutput(FuncInput)> meta_scan = [&](FuncInput inputs) {
    // an element currently residing higher up a path passing through this node
    const Word& deepest_src = inputs[0];  // corresponds to src in metascan 1
    // the maximum tree depth that an element currently residing above the
    // considered node can legally reside
    const Word& goal = inputs[1];

    Bit goal_valid = goal >= level_word;
    if (!remain_level) {
      // last level
      local_has_match = 0;
      local_deepest_one_hot = 0;
      for (uint i = 0; i < bkt_size; ++i) {
        Bit is_valid = get_bit(bkt_valids, i);
        Bit is_match = is_valid & (bkt_v_addrs[i] == v_addr);
        set_bit(local_deepest_one_hot, i, is_match);
        local_has_match |= is_match;
      }
      local_has_match &= (op == READ);
      // can drop an element if there's a vacancy or there's a match for read
      Bit can_drop = (bkt_valids != (1UL << bkt_size) - 1) | local_has_match;
      // there's an element above that can move to this level
      should_drop = can_drop & goal_valid;
      const Word& src_parent = mux(should_drop, level_word, deepest_src);
      const Word& parent_has_dest = should_drop;
      should_send = 0;
      return FuncOutput({src_parent, parent_has_dest});
    }
    FuncOutput local_deepest_info =
        get_local_deepest(position, op == READ, v_addr);
    const Word& local_deepest_level = local_deepest_info[0];

    local_has_match = local_deepest_info[2];
    // if there's a match for read, do not evict from this level
    is_local_deeper = (local_deepest_level > goal) | local_has_match;

    local_deepest_one_hot = local_deepest_info[1];

    // whether to replace the deepest element with the current element
    Bit replace_flag = (!local_has_match) & is_local_deeper;

    const Word& deepest_src_child = mux(replace_flag, deepest_src, level_word);
    const Word& goal_child = mux(replace_flag, goal, local_deepest_level);

    FuncOutput child_output =
        child_entry->meta_scan({deepest_src_child, goal_child});
    // a level which should evict an element
    const Word& src = child_output[0];
    // whether there is a level below wanting to receive an element
    // i.e. src <= level_word
    const Bit& has_dest = child_output[1];

    Bit level_eq_src_flag = level_word == src;

    // whether an element should be sent from this level
    should_send = has_dest & level_eq_src_flag;
    Bit has_vacancy = (bkt_valids != (1UL << bkt_size) - 1);
    // whether an element should be sent from a level above to a level below
    Bit should_skip = has_dest & !level_eq_src_flag;

    // compared to original code, we treat the slot to read as an empty slot
    Bit can_drop =
        ((!should_skip) & (has_vacancy | local_has_match)) | should_send;

    should_drop = can_drop & goal_valid;

    // Here we don't need to set src or deepest_src as invalid like in the
    // tri-state paper. If goal_valid is false, then should_drop is false
    // and deepest_src is not used. If level == src, then src will not match any
    // level above.
    const Word& src_parent = mux(should_drop, src, deepest_src);

    Bit parent_has_dest = should_skip | should_drop;
    const Word& parent_has_dest_word = parent_has_dest;
    return FuncOutput({src_parent, parent_has_dest_word});
  };

  /**
   * @brief Actually evict the data or read the data from the bucket
   *
   */
  std::function<FuncOutput(FuncInput)> access = [&](FuncInput inputs) {
    const Bit& hold_valid = inputs[0];
    const Word& hold_position = inputs[1];
    const Word& hold_v_addr = inputs[2];
    const Word& hold_data = inputs[3];

    // whether we should read the deepest slot
    Bit should_read = should_send | local_has_match;
    // whether we should swap with a dummy slot
    Bit replace_dummy_flag = should_drop & !should_read;

    Bit to_swap_valid = hold_valid & should_drop;
    Word to_swap_position = hold_position;
    Word to_swap_v_addr = hold_v_addr;
    Word to_swap_data = hold_data;
    // whether we should read from the bucket
    Bit should_read_from_bkt = should_read;
    if (level == 0) {
      // to avoid swapping the new element to evict and the element in the stash
      // since they are both in level 0
      should_read_from_bkt &= is_local_deeper;
    }
    for (uint i = 0; i < bkt_size; ++i) {
      Bit is_deepest = get_bit(local_deepest_one_hot, i);
      Bit read_and_remove_flag = is_deepest & should_read_from_bkt;
      Bit write_to_dummy_flag = (!get_bit(bkt_valids, i)) & replace_dummy_flag;
      // we might swap with multiple dummy slots when write to dummy flag is
      // true, but it is fine since swapping dummy with dummy is a effectively
      // no-op
      Bit swap_flag = read_and_remove_flag | write_to_dummy_flag;
      Bit bkt_valid = get_bit(bkt_valids, i);
      cond_swap(swap_flag, bkt_valid, to_swap_valid);
      set_bit(bkt_valids, i, bkt_valid);
      cond_swap(swap_flag, bkt_positions[i], to_swap_position);
      cond_swap(swap_flag, bkt_v_addrs[i], to_swap_v_addr);
      cond_swap(swap_flag, bkt_datas[i], to_swap_data);
    }
    if (!remain_level) {
      to_swap_data *= local_has_match;
      return FuncOutput({to_swap_data});
    }

    // case 1: should_drop = false, should_read = false
    // => valid = hold_valid, addr/data = hold_addr/hold_data

    // case 2: should_drop = true, should_read = false
    // => hold_valid = false

    // case 3: should_drop = true, should_read = true =>
    // => hold_valid = should_send, hold_addr data = to_swap_addr/data

    // case 4: should_drop = false, should_read = true, should_send = true =>
    // => hold_valid = true, hold_addr/data = to_swap_addr/data

    // case 5: should_drop = false, should_read = true, should_send = false =>
    // => hold_valid = hold_valid, hold_addr/data = hold_addr/data
    Bit child_hold_valid = should_send | (hold_valid & !should_drop);
    const Word& child_hold_position =
        mux(should_send, hold_position >> log_way, to_swap_position >> log_way);
    const Word& child_hold_v_addr =
        mux(should_send, hold_v_addr, to_swap_v_addr);
    const Word& child_hold_data = mux(should_send, hold_data, to_swap_data);

    FuncInput child_inputs = {child_hold_valid, child_hold_position,
                              child_hold_v_addr, child_hold_data};
    FuncOutput child_output = child_entry->access(child_inputs);

    const Word& output_data =
        mux(local_has_match, child_output[0], to_swap_data);
    return FuncOutput({output_data});
  };

  void print_tree() {
    std::cout << "level: " << level << std::endl;
    for (uint i = 0; i < bkt_size; ++i) {
      uint64_t bkt_valid_val = get_bit(bkt_valids, i);
      uint64_t bkt_positions_val = bkt_positions[i];
      uint64_t bkt_addrs_val = bkt_v_addrs[i];
      uint64_t bkt_datas_val = bkt_datas[i];
      if (bkt_valid_val) {
        std::cout << " position: " << bkt_positions_val
                  << " v_addr: " << bkt_addrs_val << " data: " << bkt_datas_val
                  << std::endl;
      }
    }
    if (children) {
      for (uint way = 0; way < 1u << log_way; ++way) {
        children[way]->print_tree();
      }
    }
  }

 private:
  /**
   * @brief Fill bkts with blocks and return the overflow blocks
   *
   * @param blocks the blocks to fill the bucket
   * @param position_shift the position in blocks should be right shifted by
   * this
   * @return std::vector<PlainTextBlock>
   */
  std::vector<PlainTextBlock> fill_bkt(const BlocksInPos& blocks,
                                       uint position_shift) {
    for (uint i = 0; i < bkt_size; ++i) {
      if (i < blocks.size()) {
        set_bit(bkt_valids, i, 1);
        bkt_positions[i] =
            slice(blocks[i].position >> position_shift, position_width);
        bkt_v_addrs[i] = blocks[i].v_addr;
        bkt_datas[i] = blocks[i].data;
      } else {
        set_bit(bkt_valids, i, 0);
        bkt_positions[i] = 0;
        bkt_v_addrs[i] = 0;
        bkt_datas[i] = 0;
      }
    }
    if (blocks.size() < bkt_size) {
      return std::vector<PlainTextBlock>();
    }
    // overflows
    return std::vector<PlainTextBlock>(blocks.begin() + bkt_size, blocks.end());
  }

 public:
  /**
   * @brief Initialize the ORAM tree
   *
   * @param begin
   * @param end
   * @param position_shift
   * @return std::vector<PlainTextBlock>
   */
  std::vector<PlainTextBlock> init_ram_g(BlocksInPosIter begin,
                                         BlocksInPosIter end,
                                         uint position_shift) {
    if (std::distance(begin, end) == 1L) {
      // leaf node
      const BlocksInPos& blocks = *begin;
      return fill_bkt(blocks, position_shift);
    }
    Assert_eq(std::distance(begin, end), 1L << position_width);
    std::vector<PlainTextBlock> overflow_blocks;
    uint64_t child_memory_space = 1UL << (position_width - log_way);
    Assert(child_memory_space);
    // recursively initialize children
    for (uint64_t i = 0; i < (1UL << log_way); ++i) {
      std::vector<PlainTextBlock> child_overflow_blocks =
          children[i]->init_ram_g(begin + i * child_memory_space,
                                  begin + (i + 1) * child_memory_space,
                                  position_shift + log_way);
      overflow_blocks.insert(overflow_blocks.end(),
                             child_overflow_blocks.begin(),
                             child_overflow_blocks.end());
    }
    return fill_bkt(overflow_blocks, position_shift);
  }

  /**
   * @brief Find the number of log_way <= max_log_way that minimizes the number
   * of levels and keeps the log_ways in all the remaining levels are evenly
   * distributed.
   *
   * @param max_log_way
   * @param position_width
   * @return uint
   */
  static uint optimal_log_way(uint max_log_way, uint position_width) {
    uint num_level = (position_width + max_log_way - 1) / max_log_way;
    if (num_level == 0) {
      return 1;
    }
    uint log_way = (position_width + num_level - 1) / num_level;
    Assert(log_way > 0);
    return log_way;
  }

  CircuitORAMRefTree(uint64_t T, uint v_addr_width, uint word_width,
                     uint position_width, uint level_width, uint level,
                     uint stash_size, uint bkt_size, uint log_way)
      : v_addr_width(v_addr_width),
        word_width(word_width),
        level_width(level_width),
        position_width(position_width),
        level(level),
        level_word(level),
        remain_level((position_width + log_way - 1) / log_way),
        bkt_size(level == 0 ? stash_size : bkt_size),
        log_way(log_way) {
    bkt_valids = 0;
    for (uint64_t i = 0; i < this->bkt_size; ++i) {
      bkt_positions.push_back(0);
      bkt_v_addrs.push_back(0);
      bkt_datas.push_back(0);
    }
    if (remain_level) {
      uint64_t T_child = (T + (1UL << log_way) - 1) >> log_way;
      if (level == 0) {
        T_child += 1UL << (position_width - log_way);
      }
      uint child_position_width = position_width - log_way;
      uint child_log_way = optimal_log_way(log_way, child_position_width);
      children = new CircuitORAMRefTree*[1UL << log_way];
      for (uint64_t i = 0; i < (1UL << log_way); ++i) {
        children[i] = new CircuitORAMRefTree(
            T_child, v_addr_width, word_width, child_position_width,
            level_width, level + 1, stash_size, bkt_size, child_log_way);
      }
    }
  }

  ~CircuitORAMRefTree() {
    if (children) {
      for (uint i = 0; i < 1UL << log_way; ++i) {
        delete children[i];
      }
      delete[] children;
    }
  }
};

/**
 * @brief The interface of the Circuit ORAM. Connected with the
 * CircuitORAMRefTree through a direct link. For each access request, perform a
 * read request to the ORAM tree followed by one or more additional evictions.
 *
 */
struct CircuitORAMRef {
  CircuitORAMRefTree* oram;
  uint position_width;  // the width of the position word
  uint v_addr_width;    // the width of the virtual address word
  uint word_width;      // the width of the data word
  uint level_width;     // the width of the total level of the tree
  uint evict_freq;      // the frequency of additional eviction per read
  uint log_way;  // the maximum log of the number of children at each level
  uint64_t evict_counter;  // the counter for the number of evictions

  using UpdateFuncType = std::function<Word(const Word&)>;

  /**
   * @brief Perform an access to the ORAM that calls a custom update function.
   *
   * @param position the current position of the word
   * @param v_addr the virtual address of the word
   * @param new_position the new position of the word
   * @param update the update function to apply to the data
   */
  void access(const Word& position, const Word& v_addr,
              const Word& new_position, const UpdateFuncType& update) {
    const Word& read_op = CircuitORAMRefTree::Op::READ;
    const Word& evict_op = CircuitORAMRefTree::Op::EVICT;

    const Word& read_init_deepest = 0;
    const Word& read_init_goal = 0;

    oram->set_addr({read_op, position, v_addr});
    FuncOutput read_meta_scan_output =
        oram->meta_scan({read_init_deepest, read_init_goal});

    const Word& read_init_hold_valid = 0;
    const Word& read_init_hold_position = 0;
    const Word& read_init_hold_v_addr = 0;
    const Word& read_init_hold = 0;
    FuncInput read_access_input = {read_init_hold_valid,
                                   read_init_hold_position,
                                   read_init_hold_v_addr, read_init_hold};
    FuncOutput read_output_data = oram->access(read_access_input);
    for (uint i = 0; i < evict_freq; ++i) {
      const Word& evict_position = slice(evict_counter, position_width);
      const Word& evict_v_addr = 0;
      oram->set_addr({evict_op, evict_position, evict_v_addr});
      if (i == 0) {
        const Word& evict_init_deepest = 0;
        const Word& position_xor = evict_position ^ new_position;
        const Word& common_suffix_len =
            trailing_zeros(position_xor, position_width, log_way);
        const Word& evict_init_goal = common_suffix_len;

        FuncOutput evict_meta_scan_output =
            oram->meta_scan({evict_init_deepest, evict_init_goal});

        const Word& to_write_valid = 1;
        const Word& to_write_data = update(read_output_data[0]);
        // TODO: allow modifying v_addr for more efficient data movement
        FuncInput to_write_input = {to_write_valid, new_position, v_addr,
                                    to_write_data};
        FuncOutput to_write_output = oram->access(to_write_input);
      } else {
        const Word& evict_init_deepest = 0;
        const Word& evict_init_goal = 0;
        FuncOutput evict_meta_scan_output =
            oram->meta_scan({evict_init_deepest, evict_init_goal});
        const Word& evict_init_hold_valid = 0;
        const Word& evict_init_hold_position = 0;
        const Word& evict_init_hold_v_addr = 0;
        const Word& evict_init_hold = 0;
        FuncInput evict_access_input = {
            evict_init_hold_valid, evict_init_hold_position,
            evict_init_hold_v_addr, evict_init_hold};
        FuncOutput evict_output_data = oram->access(evict_access_input);
      }
      ++evict_counter;
    }
  }

  /**
   * @brief Read or write a word from the ORAM
   *
   * @param position the current position of the word
   * @param v_addr the virtual address of the word
   * @param new_position the new position of the word
   * @param is_write whether the operation is a write
   * @param new_data the new data to write
   * @return Word the old data read from the ORAM
   */
  Word read_or_write(const Word& position, const Word& v_addr,
                     const Word& new_position, const Bit& is_write,
                     const Word& new_data) {
    Word old_data_copy;
    UpdateFuncType update = [&](const Word& old_data) {
      old_data_copy = old_data;
      return mux(is_write, old_data, new_data);
    };
    access(position, v_addr, new_position, update);
    return old_data_copy;
  }

  /**
   * @brief Allows the garbler to initialize the ORAM with a set of blocks.
   * The garbler inputs the blocks assigned to each position in the oram tree.
   * The evaluator inputs an empty vector.
   * The function is optional. If not called, the garbler and evaluator agree
   * that the ORAM is empty.
   *
   * @param begin
   * @param end
   */
  void init_ram_g(BlocksInPosIter begin, BlocksInPosIter end) {
    Assert(begin < end);
    std::vector<PlainTextBlock> overflow_blocks =
        oram->init_ram_g(begin, end, 0);
    Assert(overflow_blocks.empty());
  }

  CircuitORAMRef(uint64_t num_access, uint64_t memory_space, uint word_width,
                 uint stash_size = 20, uint bkt_size = 3, uint evict_freq = 2,
                 uint max_log_way = 2)
      : position_width(log2ceil(memory_space)),
        v_addr_width(log2ceil(memory_space)),
        word_width(word_width),
        evict_freq(evict_freq),
        evict_counter(0) {
    uint64_t tree_T = num_access * (1 + evict_freq);
    log_way = CircuitORAMRefTree::optimal_log_way(max_log_way, position_width);
    uint total_level = (position_width + log_way - 1) / log_way;
    level_width = bit_width(total_level);
    oram =
        new CircuitORAMRefTree(tree_T, v_addr_width, word_width, position_width,
                               level_width, 0, stash_size, bkt_size, log_way);
  }

  uint stash_load() const {
    uint64_t valid_bits = oram->bkt_valids;
    // count one
    return __builtin_popcountll(valid_bits);
  }

  ~CircuitORAMRef() { delete oram; }
};
}  // namespace ZebraGRAM