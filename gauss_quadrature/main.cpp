#include <fmt/printf.h>

#include <array>
#include <cassert>
#include <functional>
#include <numbers>
#include <vector>

namespace gauss_quadrature {

using RealFunction = std::function<double(double)>;

struct Integral {
  double low_{};
  double high_{};
  RealFunction function_{};
};

auto compose(RealFunction a, RealFunction b) -> RealFunction {
  return [a = std::move(a), b = std::move(b)](double x) { return a(b(x)); };
}

auto multiply(RealFunction a, RealFunction b) -> RealFunction {
  return [a = std::move(a), b = std::move(b)](double x) { return a(x) * b(x); };
}

// transform integration bounds from [a, b] to [new_low, new_high]
auto transformIntegrationBounds(Integral integral, double new_low, double new_high) -> RealFunction {
  auto c_0 = (integral.high_ - integral.low_) / (new_high - new_low);
  auto c_1 = integral.low_ - new_low * c_0;

  auto g = [a = integral.low_, b = integral.high_, c_0, c_1](double x) { return c_0 * x + c_1; };

  auto dgdx = [c_0](double) { return c_0; };

  return multiply(compose(std::move(integral.function_), g), dgdx);
}

template <std::size_t N>
using Params = std::array<double, N>;

struct Pair {
  double w_;
  double x_;
};

auto getParams(std::size_t i) -> Pair {
  constexpr static std::size_t kN = 10;
  assert(i < kN);

  constexpr static auto kWeights = Params<kN>{
      0.0666713,
      0.149451,
      0.219086,
      0.269267,
      0.295524,
      0.295524,
      0.269267,
      0.219086,
      0.149451,
      0.0666713,
  };

  constexpr static auto kArguments = Params<kN>{
      -0.973907,
      -0.865063,
      -0.679410,
      -0.433395,
      -0.148874,
      0.148874,
      0.433395,
      0.679410,
      0.865063,
      0.973907,
  };

  return {kWeights[i], kArguments[i]};
}

auto gaussQuadrature(Integral integral) -> double {
  constexpr std::size_t kN = 10;
  auto result = 0.0;

  auto new_integral = transformIntegrationBounds(std::move(integral), -1, 1);

  for (auto i = 0; i < 10; ++i) {
    const auto [w, x] = getParams(i);
    result += w * new_integral(x);
  }

  return result;
}

}  // namespace gauss_quadrature

auto cmp(double a, double b, double epsilon = 1e-6) -> bool { return std::abs(a - b) < epsilon; }

auto main() -> int {
  auto square_integral = gauss_quadrature::Integral{-2, 1, [](double x) { return x * x; }};
  auto square_result = gauss_quadrature::gaussQuadrature(square_integral);
  fmt::print("square = {}\n", square_result);

  auto sine_integral = gauss_quadrature::Integral{0, 2 * std::numbers::pi, [](double x) { return std::sin(x); }};
  auto sine_result = gauss_quadrature::gaussQuadrature(sine_integral);
  fmt::print("sine = {}\n", sine_result);
}
