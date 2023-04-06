#include <numbers>
#include "jr_numeric/integrals/gauss_quadrature.hpp"

auto main() -> int {
  using jr_numeric::integrals::gaussQuadrature;
  using jr_numeric::integrals::Integral;

  auto square_integral = Integral<double>{0, 2 * std::numbers::pi, [](double x) { return std::sin(x); }};

  auto square_result = gaussQuadrature<double>(square_integral);

  fmt::print("Square result: {}\n", square_result);
}
