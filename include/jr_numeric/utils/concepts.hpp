#pragma once
#include <concepts>

namespace concepts {

template <typename T>
concept FloatingPoint = std::floating_point<std::remove_cvref_t<T>>;

}
