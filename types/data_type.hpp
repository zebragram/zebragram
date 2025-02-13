#pragma once

#include "label.hpp"
namespace PicoGRAM {

struct Gadget;

/**
 * @brief Base class of any data types.
 * Each data variable holds a pointer to a gadget that owns it.
 * Operations should be performed on data variables with the same owner.
 * These operations may read or write to the GC segment of the owner.
 * Owners can be switched through links (see link.hpp).
 */
struct DataType {
 protected:
  Gadget* owner = NULL;  // the gadget that owns this data variable

 public:
  DataType() {}

  explicit DataType(Gadget* owner) : owner(owner) { Assert(owner); };

  /**
   * @brief Return the mode of the owner gadget
   *
   * @return Mode
   */
  Mode get_mode() const;

  Gadget* get_owner() const { return owner; }

  virtual void set_owner(Gadget* owner) { this->owner = owner; }

  virtual ~DataType() {}
};
}  // namespace PicoGRAM