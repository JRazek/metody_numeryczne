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
auto combineMeasurements(auto const& function, Quantities const&... quantities) -> Quantity {
  auto values = std::make_tuple(quantities.value_...);

  auto values_combined = std::apply(function, values);

  auto uncertainties_combined_sq = ld{};

  implementation::addUncertainties(uncertainties_combined_sq, function, quantities...);

  return Quantity{
      values_combined,
      uncertainties_combined_sq,
  };
}

}  // namespace uncertainty

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

//auto main() -> int {  // NOLINT
//  using uncertainty::combineMeasurements;
//  using uncertainty::Measurement;
//  using uncertainty::partialDerivative;
//  using uncertainty::Quantity;
//  using uncertainty::setupMeasurement;
//  using utils::ld;
//
//  auto a = static_cast<Quantity>(setupMeasurement(utils::readDataset("a.txt"), 0.01));
//  auto d = static_cast<Quantity>(setupMeasurement(utils::readDataset("d.txt"), 0.00002));
//
//  auto h1 = static_cast<Quantity>(setupMeasurement(utils::readDataset("h1.txt"), 0.01));
//  auto l1 = combineMeasurements([](ld a, ld h, ld d) { return a - h - d; }, a, h1, d);
//  auto tx1h1 = static_cast<Quantity>(setupMeasurement(utils::readDataset("x1h1.txt"), 0.01));
//  auto tx2h1 = static_cast<Quantity>(setupMeasurement(utils::readDataset("x2h1.txt"), 0.01));
//  auto tx5h1 = static_cast<Quantity>(setupMeasurement(utils::readDataset("x5h1.txt"), 0.01));
//  auto tx10h1 = static_cast<Quantity>(setupMeasurement(utils::readDataset("x10h1.txt"), 0.01));
//
//  auto h2 = static_cast<Quantity>(setupMeasurement(utils::readDataset("h2.txt"), 0.01));
//  auto l2 = combineMeasurements([](ld a, ld h, ld d) { return a - h - d; }, a, h2, d);
//  auto tx10h2 = static_cast<Quantity>(setupMeasurement(utils::readDataset("x10h2.txt"), 0.01));
//
//  auto h3 = static_cast<Quantity>(setupMeasurement(utils::readDataset("h3.txt"), 0.01));
//  auto l3 = combineMeasurements([](ld a, ld h, ld d) { return a - h - d; }, a, h3, d);
//  auto tx10h3 = static_cast<Quantity>(setupMeasurement(utils::readDataset("x10h3.txt"), 0.01));
//
//  auto h4 = static_cast<Quantity>(setupMeasurement(utils::readDataset("h4.txt"), 0.01));
//  auto l4 = combineMeasurements([](ld a, ld h, ld d) { return a - h - d; }, a, h4, d);
//  auto tx10h4 = static_cast<Quantity>(setupMeasurement(utils::readDataset("x10h4.txt"), 0.01));
//
//  auto h5 = static_cast<Quantity>(setupMeasurement(utils::readDataset("h5.txt"), 0.01));
//  auto l5 = combineMeasurements([](ld a, ld h, ld d) { return a - h - d; }, a, h5, d);
//  auto tx10h5 = static_cast<Quantity>(setupMeasurement(utils::readDataset("x10h5.txt"), 0.01));
//
//  //  fmt::print("h1: {}\n", h1);
//  //  fmt::print("tx10h1: {}\n\n", tx10h1);
//  //
//  //  fmt::print("h2: {}\n", h2);
//  //  fmt::print("tx10h2: {}\n\n", tx10h2);
//  //
//  //  fmt::print("h3: {}\n", h3);
//  //  fmt::print("tx10h3: {}\n\n", tx10h3);
//  //
//  //  fmt::print("h4: {}\n", h4);
//  //  fmt::print("tx10h4: {}\n\n", tx10h4);
//  //
//  //  fmt::print("h5: {}\n", h5);
//  //  fmt::print("tx10h5: {}\n\n", tx10h5);
//
//  fmt::print("h1: {}\n", h1);
//  fmt::print("l1: {}\n\n", l1);
//
//  fmt::print("h2: {}\n", h2);
//  fmt::print("l2: {}\n\n", l2);
//
//  fmt::print("h3: {}\n", h3);
//  fmt::print("l3: {}\n\n", l3);
//
//  fmt::print("h4: {}\n", h4);
//  fmt::print("l4: {}\n\n", l4);
//
//  fmt::print("h5: {}\n", h5);
//  fmt::print("l5: {}\n\n", l5);
//}
