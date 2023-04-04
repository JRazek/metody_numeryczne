#include <fmt/printf.h>

#include <array>
#include <cassert>
#include <concepts>
#include <functional>
#include <numbers>
#include <vector>

namespace gauss_quadrature {

template <std::floating_point T>
using RealFunction = std::function<T(T)>;

template <std::floating_point T>
struct Integral {
  T low_{};
  T high_{};
  RealFunction<T> function_{};
};

template <std::floating_point T>
auto compose(RealFunction<T> a, RealFunction<T> b) -> RealFunction<T> {
  return [a = std::move(a), b = std::move(b)](T x) { return a(b(x)); };
}

template <std::floating_point T>
auto multiply(RealFunction<T> a, RealFunction<T> b) -> RealFunction<T> {
  return [a = std::move(a), b = std::move(b)](T x) { return a(x) * b(x); };
}

// transform integration bounds from [a, b] to [new_low, new_high]
template <std::floating_point T>
auto transformIntegrationBounds(Integral<T> integral, T new_low, T new_high) -> RealFunction<T> {
  auto c_0 = (integral.high_ - integral.low_) / (new_high - new_low);
  auto c_1 = integral.low_ - new_low * c_0;

  auto g = [c_0, c_1](T x) { return c_0 * x + c_1; };

  auto dgdx = [c_0](T) { return c_0; };

  return multiply<T>(compose<T>(std::move(integral.function_), g), dgdx);
}

template <std::floating_point T, std::size_t N>
using Params = std::array<T, N>;

template <std::floating_point T>
struct WeightArgument {
  T x_;
  T w_;
};

template <std::floating_point T>
consteval auto generateParams() {
  constexpr std::size_t kN = 10;

  std::array<WeightArgument<T>, kN> results{{
      {-0.973906528517171720078, 0.0666713443086881375936},
      {-0.8650633666889845107321, 0.149451349150580593146},
      {-0.6794095682990244062343, 0.219086362515982043996},
      {-0.4333953941292471907993, 0.2692667193099963550912},
      {-0.1488743389816312108848, 0.2955242247147528701739},
      {0.1488743389816312108848, 0.295524224714752870174},
      {0.4333953941292471907993, 0.269266719309996355091},
      {0.6794095682990244062343, 0.2190863625159820439955},
      {0.8650633666889845107321, 0.1494513491505805931458},
      {0.973906528517171720078, 0.0666713443086881375936},
  }};

  return results;
}

template <std::floating_point T>
auto gaussQuadrature(Integral<T> integral) -> T {
  constexpr auto kParams = generateParams<T>();

  auto result = 0.0;

  auto new_integral = transformIntegrationBounds<T>(std::move(integral), -1, 1);

  for (const auto [x, w] : kParams) {
    result += w * new_integral(x);
  }

  return result;
}

}  // namespace gauss_quadrature

template <std::floating_point T>
auto abs(T x) -> T {
  return x < 0 ? -x : x;
}

template <std::floating_point T>
auto cmp(T a, T b, T epsilon = 1e-6) -> bool {
  return abs(a - b) < epsilon;
}
