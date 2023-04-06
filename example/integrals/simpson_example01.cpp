#include <fmt/core.h>

#include <cmath>
#include <numbers>

#include "jr_numeric/integrals/simpson.hpp"
#include "jr_numeric/integrals/utils.hpp"
#include "jr_numeric/utils/concepts.hpp"

auto main() -> int {
  integrals::Integral integral{0.0, 0.0, [](double x) { return std::cos(x); }};

  auto const res = integrals::simpson(integral, 1000);

  fmt::print("Result: {}\n", res);
}
