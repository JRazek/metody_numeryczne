#pragma once

#include <concepts>
#include <functional>

#include "jr_numeric/utils/utils.hpp"

namespace integrals {

using utils::R1RealFunction;
using utils::R1RealFunctionC;

template <std::floating_point T>
struct Integral {
  T low_{};
  T high_{};
  R1RealFunction<T> function_{};
};

template <typename T>
concept IntegralC = requires(T integral) {
                      { integral.low_ } -> std::floating_point;
                      { integral.high_ } -> std::floating_point;
                      { integral.function_ } -> R1RealFunctionC;
                    };

}  // namespace integrals
