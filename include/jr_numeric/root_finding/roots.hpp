#pragma once

#include <cassert>
#include <concepts>
#include <cstdint>

#include "jr_numeric/differential/derivatives.hpp"
#include "jr_numeric/utils/concepts.hpp"
#include "jr_numeric/utils/utils.hpp"

namespace jr_numeric::roots {

using differential::derivative;

template <std::floating_point T>
auto bisection(concepts::R1RealFunction auto const& function, T low, T high, std::uint64_t n) -> T {
  auto f_low = function(low);
  auto f_high = function(high);

  assert(low < high);
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

template <concepts::FloatingPoint T>
auto newtonRaphson(concepts::R1RealFunction auto const& function, T x_0, std::uint64_t n) -> T {
  for (auto i = 0u; i < n; i++) {
    auto dfdx = derivative(function, x_0);
    auto f_x = function(x_0);
    x_0 = (dfdx * x_0 - f_x) / dfdx;
  }

  return x_0;
}

}  // namespace jr_numeric::roots
