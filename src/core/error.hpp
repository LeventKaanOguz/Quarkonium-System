#pragma once
// Numerical failures that a caller may want to survive (a collapsing
// Hamiltonian, an overlap with no usable direction) are returned as
// std::expected; programming and configuration errors throw.
#include <expected>
#include <string>
#include <utility>

namespace qgem {

struct Error {
  std::string message;
};

template <class T>
using Result = std::expected<T, Error>;

inline std::unexpected<Error> fail(std::string message) {
  return std::unexpected<Error>(Error{std::move(message)});
}

}  // namespace qgem
