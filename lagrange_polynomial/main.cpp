#include <fmt/printf.h>

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <iostream>
#include <numeric>
#include <ranges>
#include <span>
#include <thread>
#include <vector>

namespace rg = std::ranges;

struct LagrangePolynomial {
  // param, value
  using SamplesVector = std::vector<std::pair<double, double>>;

  SamplesVector samples_;

  explicit LagrangePolynomial(SamplesVector samples) : samples_(std::move(samples)) {}

  [[nodiscard]] auto neville(const double t) const noexcept -> double {
    using CoeffVec = std::vector<double>;
    CoeffVec values(samples_.size());
    rg::transform(samples_, values.begin(), &SamplesVector::value_type::second);

    auto merges = values.size();

    for (auto i = 0; i < merges; i++) {
      CoeffVec tmp(samples_.size() - 1);
      for (auto j = 0; j < values.size() - 1; j++) {
        tmp[j] = ((t - samples_[j].first) * values[j + 1] - (t - samples_[j + 1].first) * values[j + 1]) /
                 (samples_[j + i + 1].first - samples_[j].first);
      }

      values = std::move(tmp);
    }

    assert(!values.empty());

    return values[0];
  }

  [[nodiscard]] auto interpolate(const double t) const noexcept -> double {
    auto res = double{};

    for (auto i = 0u; i < samples_.size(); i++) {  // sum
      auto [t_i, res_i] = samples_[i];
      for (auto j = 0u; j < samples_.size(); j++) {  // pi
        if (i != j) {
          auto t_j = samples_[j].first;
          res_i *= (t - t_j) / (t_i - t_j);
        }
      }

      res += res_i;
    }

    return res;
  }
};

enum class Method {
  KNeville,
  KNaive,
};

[[nodiscard]] auto generate(LagrangePolynomial const& pol, std::size_t resulotion, Method const method) noexcept
    -> auto{
  auto const& samples = pol.samples_;

  auto diff = (samples.back().first - samples.front().first) / static_cast<double>(resulotion);

  LagrangePolynomial::SamplesVector interpolated_values(resulotion);

  auto arg = samples.front().first;

  for (auto& [t, res] : interpolated_values) {
    t = arg;
    res = method == Method::KNeville ? pol.neville(t) : pol.interpolate(t);
    arg += diff;
  }

  return interpolated_values;
}

auto runWithMeasure(LagrangePolynomial const& pol, std::size_t drawing_resolution, Method method) {
  auto start = std::chrono::high_resolution_clock::now();
  auto res = generate(pol, drawing_resolution, method);

  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

  fmt::print("Neville: {} us\n", duration.count());

  return res;
}

auto main() -> int {
  auto samples_count = std::size_t{};
  std::cin >> samples_count;

  auto samples = LagrangePolynomial::SamplesVector(samples_count);

  for (auto& [arg, res] : samples) {
    std::cin >> arg;
  }
  for (auto& [arg, res] : samples) {
    std::cin >> res;
  }

  auto pol = LagrangePolynomial(samples);

  constexpr auto kDrawingResolution = 100;

  auto nev = runWithMeasure(pol, kDrawingResolution, Method::KNeville);
  auto naive = runWithMeasure(pol, kDrawingResolution, Method::KNaive);
}
