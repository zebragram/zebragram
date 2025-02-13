#include "data_type.hpp"

#include "gadget.hpp"

namespace PicoGRAM {
Mode DataType::get_mode() const { return owner->get_mode(); }
}