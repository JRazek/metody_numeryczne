#pragma once

#include <fmt/core.h>

#include <chrono>
#include <cmath>
#include <numeric>
#include <vector>

#include "jr_numeric/utils/concepts.hpp"

namespace jr_numeric::statistics {

using concepts::FloatingPoint;
using concepts::ReadOnlyRange;
using concepts::ScalarField;

template <FloatingPoint T>
struct Quantity {
  T value_;
  T uncertainty_sq_;
};

template <FloatingPoint T>
struct Measurement {
  T mean_;
  T std_deviation_sq_;
  T std_uncertainty_of_mean_sq_;
  T generalized_uncertainty_sq_;

  explicit operator Quantity<T>() const {
    return Quantity{
        mean_,
        generalized_uncertainty_sq_,
    };
  }
};

template <FloatingPoint T, ReadOnlyRange<T> Range>
auto calculateMean(Range const& samples) -> T {
  auto res = std::reduce(samples.begin(), samples.end());
  return res / samples.size();
}

/**
 * @brief Calculates the sample standard deviation squared
 *
 * @tparam T - type representing a precision to calculate with
 * @param samples - samples to calculate the sample standard deviation squared from
 * @param mean - mean of the samples (may be calculated from jr_numeric::statistics::calculateMean function)
 *
 * note that this does not calculate the sample standard deviation, but the standard deviation of sample.
 * see: https://en.wikipedia.org/wiki/Standard_deviation#Sample_standard_deviation
 *
 */
template <FloatingPoint T, ReadOnlyRange<T> Range>
auto calculateSampleStdDeviationSq(Range const& samples, T mean) -> T {
  assert(samples.size() > 1);
  auto res = std::transform_reduce(
      samples.begin(), samples.end(), T{}, std::plus<>{}, [mean](auto const x) { return std::pow(x - mean, 2); });
  return res / (samples.size() - 1);
}

template <FloatingPoint T>
auto calculateStdUncertaintyOfMeanSq(T variance, std::uint64_t n) -> T {
  return variance / n;
}

template <FloatingPoint T>
auto calcualteGeneralizedUncertaintySq(T std_uncertainty_of_mean_sq, T uncertainty_of_device) -> T {
  return std_uncertainty_of_mean_sq + std::pow(uncertainty_of_device, 2) / 3;
}

template <FloatingPoint T, ReadOnlyRange<T> Range>
auto setupMeasurement(Range const& samples, T uncertainty_of_device) -> Measurement<T> {
  auto mean = calculateMean<T>(samples);
  auto variance = calculateSampleStdDeviationSq(samples, mean);
  auto std_uncertainty_of_mean_sq = calculateStdUncertaintyOfMeanSq(variance, samples.size());
  auto generalized_uncertainty_sq =
      calcualteGeneralizedUncertaintySq(std_uncertainty_of_mean_sq, uncertainty_of_device);

  return Measurement{
      mean,
      variance,
      std_uncertainty_of_mean_sq,
      generalized_uncertainty_sq,
  };
}

template <FloatingPoint T, ReadOnlyRange<Quantity<T>> Range>
auto meanQuantity(Range const& quantities) -> Quantity<T> {
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

}  // namespace jr_numeric::statistics

template <jr_numeric::concepts::FloatingPoint T>
struct fmt::formatter<jr_numeric::statistics::Measurement<T>> {
  template <typename ParseContext>
  constexpr auto parse(ParseContext& ctx) {
    return ctx.begin();
  }

  template <typename FormatContext>
  auto format(jr_numeric::statistics::Measurement<T> const& quantity, FormatContext& ctx) {
    return format_to(
        ctx.out(),
        "mean: {:.3}, std_deviation: {:.3}, std_uncertainty_of_mean: {:.3}, generalized_uncertainty: {:.3}",
        quantity.mean_,
        std::sqrt(quantity.std_deviation_sq_),
        std::sqrt(quantity.std_uncertainty_of_mean_sq_),
        std::sqrt(quantity.generalized_uncertainty_sq_));
  }
};

template <jr_numeric::concepts::FloatingPoint T>
struct fmt::formatter<jr_numeric::statistics::Quantity<T>> {
  template <typename ParseContext>
  constexpr auto parse(ParseContext& ctx) {
    return ctx.begin();
  }

  template <typename FormatContext>
  auto format(jr_numeric::statistics::Quantity<T> const& quantity, FormatContext& ctx) {
    return format_to(
        ctx.out(),
        "value: {:.3}, generalized_uncertainty: {:.3}",
        quantity.value_,
        std::sqrt(quantity.uncertainty_sq_));
  }
};
