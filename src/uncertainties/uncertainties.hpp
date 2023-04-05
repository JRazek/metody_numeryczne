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

namespace utils {

using ld = long double;
using Container = std::vector<ld>;

inline auto readDataset(std::string const& path) -> Container {
  std::fstream file(fmt::format("{}/{}", "data", path));
  if (file.fail()) {
    throw std::runtime_error("Failed to open file");
  }
  Container res;
  ld x;
  while (file >> x) {
    res.push_back(x);
  }
  return res;
}

}  // namespace utils

namespace uncertainty {

using utils::Container;
using utils::ld;

// N - # of variable to differentiate by in Args
template <std::size_t N, std::convertible_to<ld>... Args>
auto partialDerivative(auto const& function, Args&&... args) -> ld {
  auto tup = std::make_tuple(static_cast<ld>(args)...);
  auto& el = std::get<N>(tup);

  constexpr auto kEpsilon = ld{1e-10};

  el -= kEpsilon;

  const auto val = std::apply(function, tup);

  el += 2 * kEpsilon;

  const auto val_h = std::apply(function, tup);

  return (val_h - val) / (2 * kEpsilon);
}

struct Quantity {
  ld value_;
  ld uncertainty_sq_;
};

struct Measurement {
  ld mean_;
  ld variance_;
  ld std_uncertainty_of_mean_sq_;
  ld generalized_uncertainty_sq_;

  explicit operator Quantity() const {
    return Quantity{
        mean_,
        generalized_uncertainty_sq_,
    };
  }
};

inline auto calculateMean(Container const& measurements) -> ld {
  auto res = std::reduce(measurements.begin(), measurements.end());
  return res / measurements.size();
}

inline auto calculateVariance(Container const& measurements, ld mean) -> ld {
  assert(measurements.size() > 1);
  auto res = std::transform_reduce(measurements.begin(), measurements.end(), ld{}, std::plus<>{}, [mean](auto const x) {
    return std::pow(x - mean, 2);
  });
  return res / (measurements.size() - 1);
}

inline auto calculateStdUncertaintyOfMeanSq(ld variance, std::uint64_t n) -> ld { return variance / n; }

inline auto calcualteGeneralizedUncertaintySq(ld std_uncertainty_of_mean_sq, ld uncertainty_of_device) -> ld {
  return std_uncertainty_of_mean_sq + std::pow(uncertainty_of_device, 2) / 3;
}

inline auto setupMeasurement(Container const& measurements, ld uncertainty_of_device) -> Measurement {
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

template <std::size_t I = 0, std::same_as<Quantity>... Quantities>
auto addUncertainties(ld& uncertainties_combined, auto const& function, Quantities const&... quantities) -> void {
  if constexpr (I < sizeof...(Quantities)) {
    auto uncertainties_sq = std::make_tuple(quantities.uncertainty_sq_...);
    uncertainties_combined +=
        std::pow(partialDerivative<I>(function, quantities.value_...), 2) * std::get<I>(uncertainties_sq);
    addUncertainties<I + 1>(uncertainties_combined, function, quantities...);
  }
}

}  // namespace implementation

template <std::same_as<Quantity>... Quantities>
auto combineQuantities(auto const& function, Quantities const&... quantities) -> Quantity {
  auto values = std::make_tuple(quantities.value_...);

  auto values_combined = std::apply(function, values);

  auto uncertainties_combined_sq = ld{};

  implementation::addUncertainties(uncertainties_combined_sq, function, quantities...);

  return Quantity{
      values_combined,
      uncertainties_combined_sq,
  };
}

inline auto meanQuantity(std::vector<Quantity> const& quantities) -> Quantity {
  auto mean = ld{};
  auto combined_uncertainty_sq = ld{};
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

}  // namespace uncertainty

template <>
struct fmt::formatter<uncertainty::Measurement> {
  template <typename ParseContext>
  constexpr auto parse(ParseContext& ctx) {
    return ctx.begin();
  }

  template <typename FormatContext>
  auto format(uncertainty::Measurement const& quantity, FormatContext& ctx) {
    return format_to(
        ctx.out(),
        "mean: {:.3}, std_deviation: {:.3}, std_uncertainty_of_mean: {:.3}, generalized_uncertainty: {:.3}",
        quantity.mean_,
        std::sqrt(quantity.variance_),
        std::sqrt(quantity.std_uncertainty_of_mean_sq_),
        std::sqrt(quantity.generalized_uncertainty_sq_));
  }
};

template <>
struct fmt::formatter<uncertainty::Quantity> {
  template <typename ParseContext>
  constexpr auto parse(ParseContext& ctx) {
    return ctx.begin();
  }

  template <typename FormatContext>
  auto format(uncertainty::Quantity const& quantity, FormatContext& ctx) {
    return format_to(ctx.out(), "value: {:.3} \\pm {:.3}", quantity.value_, std::sqrt(quantity.uncertainty_sq_));
  }
};
