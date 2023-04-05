#include <concepts>
#include <cstdint>
#include <functional>

#include "utils/utils.hpp"

namespace simpson {

using utils::Integral;

template <std::floating_point T>
auto simpson(const Integral<T>& integral, std::size_t n) -> T {
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

}  // namespace simpson
