#include "jr_numeric/integrals/gauss_quadrature.hpp"

auto main() -> int {
  using jr_numeric::integrals::gaussQuadrature;
  using jr_numeric::integrals::Integral;

  auto square_integral = Integral<double>{-2., 1., [](double x) { return x; }};

  auto square_result = gaussQuadrature<double>(square_integral);

  fmt::print("Square result: {}\n", square_result);
}
