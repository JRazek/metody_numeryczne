#include "jr_numeric/interpolations/lagrange_polynomial.hpp"

using jr_numeric::interpolations::LagrangePolynomial;
using jr_numeric::interpolations::Method;

inline auto runWithMeasure(LagrangePolynomial const& pol, std::size_t drawing_resolution, Method method) {
  auto start = std::chrono::high_resolution_clock::now();
  auto res = generate(pol, drawing_resolution, method);

  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

  fmt::print("Neville: {} us\n", duration.count());

  return res;
}

auto main() -> int {
  auto samples_count = std::size_t{};
  std::cin >> samples_count;

  auto samples = LagrangePolynomial::SamplesVector(samples_count);

  for (auto& [arg, res] : samples) {
    std::cin >> arg;
  }
  for (auto& [arg, res] : samples) {
    std::cin >> res;
  }

  auto pol = LagrangePolynomial(samples);

  constexpr auto kDrawingResolution = 100;

  auto nev = runWithMeasure(pol, kDrawingResolution, Method::KNeville);
  auto naive = runWithMeasure(pol, kDrawingResolution, Method::KNaive);
}
