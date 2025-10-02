#include "emp_interface.hpp"

#include <algorithm>

#include "channel_util.hpp"
#include "rec_oram.hpp"
namespace ZebraGRAM {

Label to_label(const BitType& raw_bit) {
  Label label;
  memcpy(label.get_ptr(), raw_bit.label,
         std::min(sizeof(Label), sizeof(raw_bit.label)));
  return label;
}

bool Delta_set_flag = false;

void set_Delta(const BitType& Delta) {
  key_manager.set_Delta(to_label(Delta));
  Delta_set_flag = true;
}

Bit to_bit(Gadget* owner, const BitType& raw_bit) {
  Label label = to_label(raw_bit);
  return Bit(owner, label, 0, false, false);
}

Word to_word(Gadget* owner, const WordType& raw_word) {
  uint word_width = raw_word.bits.size();
  Word word(owner, word_width);
  for (uint i = 0; i < word_width; ++i) {
    word[i] = to_bit(owner, raw_word.bits[i]);
  }
  return word;
}

BitType to_raw_bit(const Label& label) {
  BitType raw_bit;
  memcpy(raw_bit.label, label.get_ptr(), sizeof(raw_bit.label));
  return raw_bit;
}

BitType to_raw_bit(const Bit& bit) { return to_raw_bit(bit.get_label()); }

WordType to_raw_word(const Word& word) {
  WordType raw_word;
  raw_word.bits.reserve(word.width());
  for (uint i = 0; i < word.width(); ++i) {
    raw_word.bits.push_back(to_raw_bit(word[i]));
  }
  return raw_word;
}

// wraps recursive oram into a gadget
struct ORAMGadget : Gadget {
  RecursiveORAM oram;

  // delete copy constructor
  ORAMGadget(const ORAMGadget&) = delete;
  ORAMGadget& operator=(const ORAMGadget&) = delete;

  Word read_or_write(const Word& addr, const Bit& is_write,
                     const Word& new_data) {
    return oram.read_or_write(addr, is_write, new_data);
  }
  ORAMGadget(Mode mode, uint64_t T, uint memory_space, uint word_width)
      : Gadget(mode, T), oram(this, memory_space, word_width, 0, T) {}
};

ORAMType::ORAMType(uint32_t address_width, uint32_t word_width,
                   uint64_t num_accesses, bool is_garbler)
    : is_garbler(is_garbler),
      addr_width(address_width),
      word_width(word_width),
      num_accesses(num_accesses) {
  internal_oram = new ORAMGadget(is_garbler ? GARBLE : EVAL, num_accesses,
                                 1UL << address_width, word_width);
}

// First, translate the input labels into our previously mocked labels.
// For the garbler, directly return the mocked output labels
// For the evaluator, run the access function and return the output encodings

WordType ORAMType::access(const WordType& addr, const BitType& is_write,
                          const WordType& new_data) {
  ORAMGadget* oram = static_cast<ORAMGadget*>(internal_oram);
  Word addr_word = to_word(oram, addr);
  Bit is_write_bit = to_bit(oram, is_write);
  Word new_data_word = to_word(oram, new_data);
  if (is_garbler) {
    Assert_less(acccess_counter, mocked_input_addrs.size());
    Word mocked_addr = to_word(oram, mocked_input_addrs[acccess_counter]);
    Bit mocked_is_write = to_bit(oram, mocked_input_is_writes[acccess_counter]);
    Word mocked_new_data =
        to_word(oram, mocked_input_new_data[acccess_counter]);

    Word addr_diff = mocked_addr ^ addr_word;
    Bit is_write_diff = mocked_is_write ^ is_write_bit;
    Word new_data_diff = mocked_new_data ^ new_data_word;
    // send diffs
    for (uint i = 0; i < addr_width; ++i) {
      io_channel.send_data(addr_diff[i].get_label().get_ptr(), LAMBDA_BYTES);
    }
    io_channel.send_data(is_write_diff.get_label().get_ptr(), LAMBDA_BYTES);
    for (uint i = 0; i < word_width; ++i) {
      io_channel.send_data(new_data_diff[i].get_label().get_ptr(),
                           LAMBDA_BYTES);
    }
    // return labels of the output old data
    const WordType& old_data = mocked_output_old_data[acccess_counter];
    ++acccess_counter;
    return old_data;
  } else {
    // receive diffs
    Word addr_diff = Word::rand_label_word(oram, addr_width);
    for (uint i = 0; i < addr_width; ++i) {
      io_channel.recv_data(addr_diff[i].get_label().get_ptr(), LAMBDA_BYTES);
    }
    Bit is_write_diff = Bit::rand_label_bit(oram);
    io_channel.recv_data(is_write_diff.get_label().get_ptr(), LAMBDA_BYTES);
    Word new_data_diff = Word::rand_label_word(oram, word_width);
    for (uint i = 0; i < word_width; ++i) {
      io_channel.recv_data(new_data_diff[i].get_label().get_ptr(),
                           LAMBDA_BYTES);
    }
    // apply diffs
    addr_word ^= addr_diff;
    is_write_bit ^= is_write_diff;
    new_data_word ^= new_data_diff;

    const Word& old_data_word =
        oram->read_or_write(addr_word, is_write_bit, new_data_word);
    const WordType& old_data = to_raw_word(old_data_word);
    return old_data;
  }
}

std::vector<uint64_t> ORAMType::measure_gadget_gc() const {
  ORAMGadget* measure_oram_gadget =
      new ORAMGadget(MEASURE, num_accesses, 1UL << addr_width, word_width);
  // measure_oram_gadget->oram.init_empty_ram(num_accesses);
  measure_oram_gadget->init_gc_ptr(GCPtr(-1));
  for (uint64_t i = 0; i < num_accesses; ++i) {
    measure_oram_gadget->inc_time();
    Word addr = Word::rand_label_word(measure_oram_gadget, addr_width);
    Bit is_write = Bit::rand_label_bit(measure_oram_gadget);
    Word new_data = Word::rand_label_word(measure_oram_gadget, word_width);
    const Word& old_data =
        measure_oram_gadget->read_or_write(addr, is_write, new_data);
  }
  GCPtr actual_end_ptr = measure_oram_gadget->garble_callees();

  std::vector<GCPtr> gc_ptrs = measure_oram_gadget->get_init_gc_ptrs();
  std::vector<uint64_t> gc_offsets(gc_ptrs.size() + 1);
  std::transform(gc_ptrs.begin(), gc_ptrs.end(), gc_offsets.begin(),
                 [](const GCPtr& gc_ptr) { return gc_ptr.get_offset(); });
  gc_offsets.back() = actual_end_ptr.get_offset();
  delete measure_oram_gadget;
  return gc_offsets;
}

void ORAMType::garble_oram() {
  mocked_input_addrs.resize(num_accesses);
  mocked_input_is_writes.resize(num_accesses);
  mocked_input_new_data.resize(num_accesses);
  mocked_output_old_data.resize(num_accesses);
  ORAMGadget* oram_gadget = static_cast<ORAMGadget*>(internal_oram);
  Assert_eq(oram_gadget->get_mode(), GARBLE);
  for (uint64_t i = 0; i < num_accesses; ++i) {
    oram_gadget->inc_time();
    Word addr = Word::rand_label_word(oram_gadget, addr_width);
    Bit is_write = Bit::rand_label_bit(oram_gadget);
    Word new_data = Word::rand_label_word(oram_gadget, word_width);
    const Word& old_data = oram_gadget->read_or_write(addr, is_write, new_data);
    mocked_input_addrs[i] = to_raw_word(addr);
    mocked_input_is_writes[i] = to_raw_bit(is_write);
    mocked_input_new_data[i] = to_raw_word(new_data);
    mocked_output_old_data[i] = to_raw_word(old_data);
  }
  oram_gadget->garble_callees();
}

template <typename ChannelType>
static void send_uint64(ChannelType io_channel, uint64_t value) {
  uint64_t value_le = htole64(value);
  io_channel.send_data((uint8_t*)&value_le, sizeof(value_le));
}

template <typename ChannelType>
static uint64_t recv_uint64(ChannelType io_channel) {
  uint64_t value_le;
  io_channel.recv_data((uint8_t*)&value_le, sizeof(value_le));
  return le64toh(value_le);
}

// The garbler measures the garbled circuit and first sends the meta
// information to the evaluator. Then, the garbler garbles the ORAM and sends
// the GC to the evaluator at the same time. The evaluator stores the GC
// locally.

void ORAMType::initialize(ChannelType channel) {
  this->io_channel = channel;
  ORAMGadget* oram_gadget = static_cast<ORAMGadget*>(internal_oram);
  if (is_garbler) {
    Assert(num_accesses > 0);
    // measure the number of bytes to be sent
    std::vector<uint64_t> gc_offsets = measure_gadget_gc();
    uint64_t gc_size = gc_offsets.back();
    // send meta information
    std::cout << "prepare to send gc size: " << gc_size << std::endl;
    send_uint64(io_channel, gc_size);
    gc_offsets.pop_back();
    send_uint64(io_channel, gc_offsets.size());
    for (uint64_t gc_offset : gc_offsets) {
      send_uint64(io_channel, gc_offset);
    }
    int fid = Channel::ChannelManager::get_instance().add_channel(io_channel);
    oram_gadget->init_gc_ptr(GCPtr(fid));
    oram_gadget->oram.init_empty_ram();
    garble_oram();

  } else {
    uint64_t gc_size, gc_offsets_count;
    gc_size = recv_uint64(io_channel);
    Assert_less(gc_size, 1UL << 40);
    gc_offsets_count = recv_uint64(io_channel);
    uint64_t file_name_uid = secure_random_uint64();
    std::string file_name = "oram_gc_" + std::to_string(file_name_uid) + ".bin";
    int fid = open_file(file_name.c_str(), (size_t)gc_size);
    Assert(fid >= 0);
    std::vector<GCPtr> gc_ptrs(gc_offsets_count, GCPtr(fid));
    for (uint64_t i = 0; i < gc_offsets_count; ++i) {
      gc_ptrs[i].set_offset(recv_uint64(io_channel));
    }
    // receive the GC through the channel
    const size_t chunk_size = 1UL << 20;
    uint8_t buffer[chunk_size];
    // check current offset of fid is 0
    Assert_eq(lseek(fid, 0, SEEK_CUR), 0);
    for (uint64_t i = 0; i < gc_size; i += chunk_size) {
      uint64_t size = std::min(chunk_size, gc_size - i);
      io_channel.recv_data(buffer, size);
      ssize_t write_size = write(fid, buffer, size);
      Assert_eq(write_size, (ssize_t)size);
    }
    std::cout << "eval finished receiving gc" << std::endl;
    oram_gadget->set_init_gc_ptrs(gc_ptrs);
    oram_gadget->oram.init_empty_ram();
  }
}

ORAMType::~ORAMType() {
  if (internal_oram) {
    delete static_cast<ORAMGadget*>(internal_oram);
  }
  internal_oram = nullptr;
}

}  // namespace ZebraGRAM