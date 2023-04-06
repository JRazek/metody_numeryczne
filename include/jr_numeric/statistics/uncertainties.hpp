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
#include "jr_numeric/utils/concepts.hpp"

namespace jr_numeric::uncertainty {

template <typename T>
using Container = std::vector<T>;

using concepts::FloatingPoint;
using concepts::ScalarField;

template <FloatingPoint T>
struct Quantity {
  T value_;
  T uncertainty_sq_;
};

template <FloatingPoint T>
struct Measurement {
  T mean_;
  T variance_;
  T std_uncertainty_of_mean_sq_;
  T generalized_uncertainty_sq_;

  explicit operator Quantity<T>() const {
    return Quantity{
        mean_,
        generalized_uncertainty_sq_,
    };
  }
};

template <FloatingPoint T>
auto calculateMean(Container<T> const& measurements) -> T {
  auto res = std::reduce(measurements.begin(), measurements.end());
  return res / measurements.size();
}

template <FloatingPoint T>
auto calculateVariance(Container<T> const& measurements, T mean) -> T {
  assert(measurements.size() > 1);
  auto res = std::transform_reduce(measurements.begin(), measurements.end(), T{}, std::plus<>{}, [mean](auto const x) {
    return std::pow(x - mean, 2);
  });
  return res / (measurements.size() - 1);
}

template <FloatingPoint T>
auto calculateStdUncertaintyOfMeanSq(T variance, std::uint64_t n) -> T {
  return variance / n;
}

template <FloatingPoint T>
auto calcualteGeneralizedUncertaintySq(T std_uncertainty_of_mean_sq, T uncertainty_of_device) -> T {
  return std_uncertainty_of_mean_sq + std::pow(uncertainty_of_device, 2) / 3;
}

template <FloatingPoint T>
auto setupMeasurement(Container<T> const& measurements, T uncertainty_of_device) -> Measurement<T> {
  auto mean = calculateMean(measurements);
  auto variance = calculateVariance(measurements, mean);
  auto std_uncertainty_of_mean_sq = calculateStdUncertaintyOfMeanSq(variance, measurements.size());
  auto generalized_uncertainty_sq =
      calcualteGeneralizedUncertaintySq(std_uncertainty_of_mean_sq, uncertainty_of_device);

  return Measurement{
      mean,
      variance,
      std_uncertainty_of_mean_sq,
      generalized_uncertainty_sq,
  };
}

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

template <FloatingPoint T>
inline auto meanQuantity(Container<Quantity<T>> const& quantities) -> Quantity<T> {
  auto mean = T{};
  auto combined_uncertainty_sq = T{};
  for (auto const& [value, uncertainty_sq] : quantities) {
    mean += value;
    combined_uncertainty_sq += uncertainty_sq;
  }

  mean /= quantities.size();
  combined_uncertainty_sq /= (quantities.size() * quantities.size());

  return Quantity{
      mean,
      combined_uncertainty_sq,
  };
}

}  // namespace jr_numeric::uncertainty

template <jr_numeric::concepts::FloatingPoint T>
struct fmt::formatter<jr_numeric::uncertainty::Measurement<T>> {
  template <typename ParseContext>
  constexpr auto parse(ParseContext& ctx) {
    return ctx.begin();
  }

  template <typename FormatContext>
  auto format(jr_numeric::uncertainty::Measurement<T> const& quantity, FormatContext& ctx) {
    return format_to(
        ctx.out(),
        "mean: {:.3}, std_deviation: {:.3}, std_uncertainty_of_mean: {:.3}, generalized_uncertainty: {:.3}",
        quantity.mean_,
        std::sqrt(quantity.variance_),
        std::sqrt(quantity.std_uncertainty_of_mean_sq_),
        std::sqrt(quantity.generalized_uncertainty_sq_));
  }
};

template <jr_numeric::concepts::FloatingPoint T>
struct fmt::formatter<jr_numeric::uncertainty::Quantity<T>> {
  template <typename ParseContext>
  constexpr auto parse(ParseContext& ctx) {
    return ctx.begin();
  }

  template <typename FormatContext>
  auto format(jr_numeric::uncertainty::Quantity<T> const& quantity, FormatContext& ctx) {
    return format_to(ctx.out(), "value: {:.3} \\pm {:.3}", quantity.value_, std::sqrt(quantity.uncertainty_sq_));
  }
};
