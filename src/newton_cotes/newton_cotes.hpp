#include <fmt/core.h>

#include <cmath>
#include <concepts>
#include <cstdint>
#include <functional>

namespace newton_cotes {

template <std::floating_point T>
using RealFunction = std::function<T(T)>;

template <std::floating_point T>
struct Integral {
  T low_{};
  T high_{};
  RealFunction<T> function_{};
};

template <std::floating_point T>
[[nodiscard]] auto riemannSum(RealFunction<T> const& step_function, T low, T high, T step) -> T {
  auto res = T{};
  const std::int64_t n = (high - low) / step;

  for (auto i = 0u; i < n; ++i) {
    res += step_function(low + i * step);
  }

  return res;
}

template <std::floating_point T>
[[nodiscard]] auto riemannIntegral(Integral<T> const& integral, T dx) -> T {
  return dx * riemannSum(integral.function_, integral.low_, integral.high_, dx);
}

template <std::floating_point T>
[[nodiscard]] auto newtonCotes(Integral<T> integral, T dx) -> T {
  auto transformed = Integral<T>{
      .low_ = integral.low_ + dx,
      .high_ = integral.high_ - dx,
      .function_ = std::move(integral.function_),
  };

  return riemannIntegral(transformed, dx) +
         dx * 0.5 * (transformed.function_(integral.low_) + transformed.function_(integral.high_));
}

}  // namespace newton_cotes
