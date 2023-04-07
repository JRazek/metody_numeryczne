#pragma once

#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/printf.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <concepts>
#include <cstdint>
#include <fstream>
#include <functional>
#include <numeric>
#include <vector>

#include "jr_numeric/differential/derivatives.hpp"
#include "jr_numeric/statistics/utils.hpp"
#include "jr_numeric/utils/concepts.hpp"

namespace jr_numeric::statistics {

namespace implementation {

template <
    std::size_t I = 0,
    FloatingPoint T,
    std::same_as<Quantity<T>>... Quantities,
    ScalarField<sizeof...(Quantities)> Function>
auto addUncertainties(T& uncertainties_combined, Function const& function, Quantities const&... quantities) -> void {
  if constexpr (I < sizeof...(Quantities)) {
    auto uncertainties_sq = std::make_tuple(quantities.uncertainty_sq_...);
    uncertainties_combined +=
        std::pow(differential::partialDerivative<I>(function, quantities.value_...), 2) * std::get<I>(uncertainties_sq);
    addUncertainties<I + 1>(uncertainties_combined, function, quantities...);
  }
}

}  // namespace implementation

template <FloatingPoint T, std::same_as<Quantity<T>>... Quantities, ScalarField<sizeof...(Quantities)> Function>
auto combineQuantities(Function const& function, Quantities const&... quantities) -> Quantity<T> {
  auto values = std::make_tuple(quantities.value_...);

  auto values_combined = std::apply(function, values);

  auto uncertainties_combined_sq = T{};

  implementation::addUncertainties(uncertainties_combined_sq, function, quantities...);

  return Quantity{
      values_combined,
      uncertainties_combined_sq,
  };
}

}  // namespace jr_numeric::statistics

