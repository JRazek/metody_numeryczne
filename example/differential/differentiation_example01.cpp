#include <fmt/core.h>

#include <numbers>

#include "jr_numeric/differential/derivatives.hpp"

auto main() -> int {
  using f = float;
  using d = double;
  using ld = long double;

  auto funciton = [](double x, double y) { return std::sin(x) + y; };

  auto res1 = differential::partialDerivative<0, f>(funciton, 2 * std::numbers::pi, 2.0);

  auto res2 = differential::partialDerivative<0, d>(funciton, 2 * std::numbers::pi, 2.0);

  auto res3 = differential::partialDerivative<0, ld>(funciton, 2 * std::numbers::pi, 2.0);

  fmt::print("differentiation of f(x) = x^2 + y^2 wrt x\n");
  fmt::print("float: {}, double: {}, long double: {}\n", res1, res2, res3);
}
