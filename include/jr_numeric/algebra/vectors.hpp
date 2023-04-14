#pragma once

#include <cmath>
#include <cstdlib>
#include <limits>
#include <tuple>
#include <type_traits>

#include "jr_numeric/algebra/matrix.hpp"
#include "jr_numeric/utils/concepts.hpp"

namespace jr_numeric::algebra {

template <std::size_t N, typename T>
struct Vector : public Matrix<N, N, T> {
  using Matrix<N, N, T>::Matrix;
};

template <std::size_t N, concepts::Number T>
struct Polynomial : Vector<N, T> {
  using Vector<N, T>::Vector;
};

}  // namespace jr_numeric::algebra
