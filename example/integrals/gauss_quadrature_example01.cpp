#include "jr_numeric/integrals/gauss_quadrature.hpp"

auto main() -> int {
  auto square_integral = integrals::Integral<double>{-2., 1., [](double x) { return x; }};

  auto square_result = integrals::gaussQuadrature<double>(square_integral);

  fmt::print("Square result: {}\n", square_result);
}
