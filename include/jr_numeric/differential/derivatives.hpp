#pragma once

#include <cmath>
#include <cstdlib>
#include <limits>
#include <tuple>
#include <type_traits>

#include "jr_numeric/utils/concepts.hpp"

namespace differential {
using concepts::FloatingPoint;
using concepts::ScalarField;

namespace implementation {

template <concepts::FloatingPoint T>
consteval auto differentiationEpsilon() -> T {
  return std::sqrt(std::numeric_limits<T>::epsilon());
}

}  // namespace implementation
// N - # of variable to differentiate by in Args
template <std::size_t N, FloatingPoint... Args, FloatingPoint T = std::common_type_t<Args...>>
auto partialDerivative(ScalarField<sizeof...(Args)> auto const& function, Args&&... args) -> T {
  auto tup = std::make_tuple(static_cast<T>(args)...);
  auto& el = std::get<N>(tup);

  constexpr auto kEpsilon = implementation::differentiationEpsilon<T>();

  el -= kEpsilon;

  const auto val = std::apply(function, tup);

  el += 2 * kEpsilon;

  const auto val_h = std::apply(function, tup);

  return (val_h - val) / (2 * kEpsilon);
}

}  // namespace differential
