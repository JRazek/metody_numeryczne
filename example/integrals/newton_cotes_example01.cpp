#include "jr_numeric/integrals/newton_cotes.hpp"

// increase integral's precision by 2 in each iteration without repeating computations
template <std::floating_point T>
auto excercise(integrals::Integral<T> integral, std::int64_t n) -> void {
  using integrals::Integral;
  using integrals::riemannIntegral;

  // returns integral / delta - evaluations on boundary points
  auto a_0 = integral.low_;
  auto b_0 = integral.high_;

  auto dx_0 = (b_0 - a_0) / n;

  auto res = T{};

  constexpr auto kIterations = 10;

  auto dx = dx_0;
  for (auto i = 0; i < kIterations; i++) {
    res += riemannIntegral(integral, dx) / dx;

    dx *= 0.5;
  }

  res *= dx;

  fmt::print("Excercise: {}\n", res);
}

auto main() -> int {
  using integrals::Integral;
  using integrals::newtonCotes;
  using integrals::riemannIntegral;

  using ld = long double;

  Integral<ld> integral{
      .low_ = 0.,
      .high_ = 10.,
      .function_ = [](auto x) { return x; },
  };

  auto res = riemannIntegral<ld>(integral, 0.00001);
  fmt::print("Riemann integral: {}\n", res);

  res = newtonCotes<ld>(integral, 0.001);
  fmt::print("Newton-Cotes integral: {}\n", res);

  excercise(integral, 100000);
}
