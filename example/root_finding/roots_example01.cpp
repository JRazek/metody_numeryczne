#include <fmt/core.h>

#include <cmath>
#include <numbers>

#include "jr_numeric/root_finding/roots.hpp"

auto main() -> int {
  using std::numbers::pi;
  auto res = roots::bisection<double>([](double x) { return std::sin(x); }, 2, 4, 100);

  fmt::print("Result: {}\n", res);
}
