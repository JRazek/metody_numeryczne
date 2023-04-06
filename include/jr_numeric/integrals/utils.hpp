#pragma once

#include <concepts>
#include <functional>

#include "jr_numeric/utils/meta.hpp"
#include "jr_numeric/utils/utils.hpp"

namespace jr_numeric::integrals {

using concepts::R1RealFunction;
using meta::IntegralFunctionResult;
using meta::RealFunctionResult;

template <std::floating_point T>
struct Integral {
  using R1RealFunction = std::function<T(T)>;
  T low_{};
  T high_{};
  R1RealFunction function_{};
};

template <std::floating_point A, std::floating_point B, typename C>
Integral(A a, B b, C c) -> Integral<std::common_type_t<A, B>>;

// note - std::integral is about integer not calculus integral
template <std::integral A, std::integral B, typename C>
Integral(A a, B b, C c) -> Integral<double>;

}  // namespace jr_numeric::integrals
