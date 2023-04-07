#pragma once

#include <chrono>

#include "jr_numeric/statistics/utils.hpp"
#include "jr_numeric/utils/concepts.hpp"

namespace jr_numeric::statistics {

template <FloatingPoint T>
auto tStudentTest(
    Measurement<T> const& m1, std::size_t sample_size1, Measurement<T> const& m2, std::size_t sample_size2) -> T {
  const auto s_p = std::sqrt(
      ((sample_size1 - 1) * m1.std_deviation_sq_ + (sample_size2 - 1) * m2.std_deviation_sq_) /
      (sample_size1 + sample_size2 - 2));

  return (m1.mean_ - m2.mean_) / (s_p * std::sqrt(T{1} / sample_size1 + T{1} / sample_size2));
}

}  // namespace jr_numeric::statistics
