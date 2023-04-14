#include <fmt/core.h>

#include <cmath>
#include <numbers>

#include "jr_numeric/root_finding/roots.hpp"

auto main() -> int {
  using jr_numeric::roots::bisection;
  using jr_numeric::roots::newtonRaphson;
  using std::numbers::pi;

  auto function = [](double x) { return (std::pow(x, 3) - 9) / (std::log(x) - 1); };

  auto res1 = bisection<double>(function, 1.8, 2.1, 10);  // approx solution

  auto res2 = newtonRaphson<long double>(function, res1, 10);  // better solution

  fmt::print("Bisection: {}\tNewton-Raphson: {}\t", res1, res2);
}
