#pragma once

#include <concepts>
#include <functional>

#include "jr_numeric/utils/utils.hpp"

namespace integrals {

using utils::RealFunction;

template <std::floating_point T>
struct Integral {
  T low_{};
  T high_{};
  RealFunction<T> function_{};
};

}  // namespace integrals
