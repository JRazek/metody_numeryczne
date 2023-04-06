#pragma once

#include <concepts>
#include <cstdint>
#include <functional>
#include <type_traits>

#include "jr_numeric/integrals/utils.hpp"
#include "jr_numeric/utils/meta.hpp"
#include "jr_numeric/utils/utils.hpp"

namespace integrals {

using concepts::R1RealFunction;

template <concepts::Integral IntegralType>
[[nodiscard]] auto simpson(IntegralType const& integral, std::size_t n) -> meta::IntegralFunctionResult<IntegralType> {
  auto dx = (integral.high_ - integral.low_) / n;

  auto const& function = integral.function_;

  auto const low = integral.low_;
  auto const high = integral.high_;

  auto res = function(low) + function(high);

  for (auto i = 0u; i < n / 2; i++) {
    res += 4 * function(low + dx * (i * 2 + 1));
  }
  for (auto i = 0u; i < n / 2 - 1; i++) {
    res += 2 * function(low + dx * (i * 2 + 2));
  }

  return res * dx / 3;
}

}  // namespace integrals
