#pragma once

#include <cassert>
#include <cstdint>

#include "jr_numeric/utils/utils.hpp"

namespace jr_numeric::roots {

template <std::floating_point T>
auto bisection(concepts::R1RealFunction auto const& function, T low, T high, std::uint64_t n) -> T {
  assert(low < high);
  auto f_low = function(low);
  auto f_high = function(high);
  assert(f_low * f_high < 0);

  auto mid = T{};

  for (auto i = 0u; i < n + 1; i++) {
    mid = (low + high) / 2;
    auto f_mid = function(mid);

    if (f_mid * f_high < 0)
      low = mid;
    else
      high = mid;
  }

  return mid;
}

}  // namespace jr_numeric::roots
