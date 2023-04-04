#include <fmt/core.h>

#include <algorithm>
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

auto readDataset(std::string const& path) -> Container {
  std::fstream file(path);
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

  const auto val = std::apply(function, tup);

  el += kEpsilon;

  auto val_h = std::apply(function, tup);

  return (val_h - val) / kEpsilon;
}

struct Measurement {
  Container measurements_;
  ld mean_;
  ld variance_;
  ld std_uncertainty_of_mean_sq_;
  ld generalized_uncertainty_;
};

auto calculateMean(Container const& measurements) -> ld {
  auto res = std::reduce(measurements.begin(), measurements.end());
  return res / measurements.size();
}

auto calculateVariance(Container const& measurements, ld mean) -> ld {
  auto res = std::transform_reduce(measurements.begin(), measurements.end(), ld{}, std::plus<>{}, [mean](auto const x) {
    return (x - mean) * (x - mean);
  });
  return res / measurements.size();
}

auto calculateStdUncertaintyOfMeanSq(ld variance, std::uint64_t n) -> ld { return variance / (n + 1); }

auto calcualteGeneralizedUncertaintySq(ld std_uncertainty_of_mean_sq, ld uncertainty_of_device) -> ld {
  return std_uncertainty_of_mean_sq + uncertainty_of_device * uncertainty_of_device / 3;
}

auto calculateMeasurement(Container measurements, ld uncertainty_of_device) -> Measurement {
  auto mean = calculateMean(measurements);
  auto variance = calculateVariance(measurements, mean);
  auto std_uncertainty_of_mean_sq = calculateStdUncertaintyOfMeanSq(variance, measurements.size());
  auto generalized_uncertainty = calcualteGeneralizedUncertaintySq(std_uncertainty_of_mean_sq, uncertainty_of_device);

  return Measurement{
      std::move(measurements),
      mean,
      variance,
      std_uncertainty_of_mean_sq,
      generalized_uncertainty,
  };
}

}  // namespace uncertainty

auto main() -> int {
  using uncertainty::partialDerivative;
  using utils::ld;

  std::fstream file("data.txt");
  auto res = partialDerivative<0>([](ld x, ld y) -> ld { return std::sin(x * y) + y; }, 0, 1);

  fmt::print("res: {}\n", res);
}
