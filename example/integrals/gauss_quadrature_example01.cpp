#include "integrals/gauss_quadrature.hpp"

auto main() -> int {
  auto square_integral = integrals::Integral<double>{-2, 1, [](double x) { return x * x; }};
  auto square_result = integrals::gaussQuadrature(square_integral);
  fmt::print("square = {}\n", square_result);

  auto sine_integral = integrals::Integral<double>{0, std::numbers::pi / 2, [](double x) { return std::sin(x); }};
  auto sine_result = integrals::gaussQuadrature(sine_integral);
  fmt::print("sine = {}\n", sine_result);
}
