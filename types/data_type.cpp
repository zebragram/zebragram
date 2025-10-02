#include "data_type.hpp"

#include "gadget.hpp"

namespace ZebraGRAM {
Mode DataType::get_mode() const { return owner->get_mode(); }
}