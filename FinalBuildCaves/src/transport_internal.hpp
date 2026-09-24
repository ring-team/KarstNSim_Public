#pragma once
#include "fbs/caves.hpp"
namespace fbs::caves::detail {
// Only for a value that has just passed generate()'s complete validation.
std::string serialize_validated_region(const Region& region);
}
