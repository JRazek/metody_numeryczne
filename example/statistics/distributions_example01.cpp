#include <fmt/core.h>

#include <cmath>
#include <cstdint>
#include <map>
#include <numbers>
#include <random>
#include <sciplot/sciplot.hpp>

#include "jr_numeric/statistics/distributions.hpp"

using Rolls = std::map<std::uint8_t, std::uint8_t>;
using RollSum = std::map<std::uint8_t, std::uint64_t>;

namespace sc = sciplot;

auto dice() -> std::uint8_t {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(1, 6);
  return dis(gen);
}

auto roll(std::uint32_t n) -> Rolls {
  Rolls rolls;
  for (auto i = 0u; i < n; ++i) {
    ++rolls[dice()];
  }
  return rolls;
}

auto sumRolls(std::uint32_t n, std::uint64_t dice_num) -> RollSum {
  RollSum roll_sums;
  for (auto i = 0u; i < n; ++i) {
    auto sum = 0u;
    for (auto j = 0u; j < dice_num; ++j) {
      sum += dice();
    }
    ++roll_sums[sum];
  }

  return roll_sums;
}

auto drawDistribution(sc::Plot2D& plot, RollSum const& rolls) -> void {
  auto x = std::vector<std::uint64_t>();
  auto y = std::vector<std::uint64_t>();

  for (auto const& [key, value] : rolls) {
    x.push_back(key);
    y.push_back(value);
  }

  plot.drawBoxes(x, y).fillSolid().fillColor("green").fillIntensity(0.5).borderShow().labelNone();
}

auto showRolls(RollSum const& rolls) -> void {
  auto plot = sc::Plot2D();
  plot.legend().atTopLeft();

  // Set the y label and its range
  plot.ylabel("occurences");
  plot.yrange(0., 1000);
  plot.xlabel("sum of dice");
  plot.xtics().fontSize(12);

  drawDistribution(plot, rolls);

  sc::Figure fig = {{plot}};
  sc::Canvas canvas = {{fig}};

  canvas.size(1000, 500);
  canvas.show();
}

auto main() -> int {
  using jr_numeric::statistics::tStudentTest;
  using std::numbers::pi;

  auto rolls = sumRolls(10000, 100);

  showRolls(rolls);
}
