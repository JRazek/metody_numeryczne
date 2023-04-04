#include "lagrange_polynomial/lagrange_polynomial.hpp"

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
