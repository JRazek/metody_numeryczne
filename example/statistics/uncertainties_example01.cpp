#include <cstdint>
#include <ranges>
#include <sciplot/Canvas.hpp>
#include <sciplot/Figure.hpp>
#include <sciplot/Plot2D.hpp>
#include <sciplot/sciplot.hpp>
#include <valarray>

#include "statistics/uncertainties.hpp"
#include "utils/utils.hpp"

using sciplot::Canvas;
using sciplot::Figure;
using sciplot::Plot2D;
using sciplot::Strings;
using sciplot::Vec;
using uncertainty::combineQuantities;
using uncertainty::meanQuantity;
using uncertainty::Measurement;
using uncertainty::Quantity;
using uncertainty::setupMeasurement;
using utils::ld;
namespace rg = std::ranges;

auto readDatasset(std::string const& path) -> std::vector<ld> { return utils::readDataset<ld>(path); }

auto generateRange(ld low, ld high, ld increment) -> std::valarray<uncertainty::ld> {
  assert(high > low);
  auto n = static_cast<std::size_t>((high - low) / increment);

  auto res = std::valarray<uncertainty::ld>(n);
  for (auto i = 0u; i < n; i++) {
    res[i] = low + i * increment;
  }

  return res;
}

auto drawDistribution(Plot2D& plot, utils::Container const& data, double resolution) {
  auto [min, max] = rg::minmax(data);
  std::map<std::int64_t, std::uint64_t> distribution;

  for (const auto x : data) {
    auto key = static_cast<std::int64_t>(std::floor(x / resolution));
    distribution[key]++;
  }

  auto n = (max - min) / resolution;
  for (auto i = 0; i < n; i++) {
    auto key = static_cast<std::int64_t>(std::floor(min / resolution + i));

    if (!distribution.contains(key)) {
      distribution[key] = 0;
    }
  }

  auto tmp_x = std::vector<double>();
  auto tmp_y = std::vector<double>();

  for (auto [k, v] : distribution) {
    tmp_x.push_back(static_cast<double>(k));
    tmp_y.push_back(static_cast<double>(v));
  }

  auto x = Vec(tmp_x.data(), tmp_x.size());
  auto y = Vec(tmp_y.data(), tmp_y.size());

  x *= resolution;

  return plot.drawBoxes(x, y).fillSolid().fillColor("green").fillIntensity(0.5).borderShow().labelNone();
}

auto draw(std::string const& dataset_name) -> void {
  auto data = readDatasset(dataset_name);

  auto plot = Plot2D();
  plot.legend().atTopLeft();

  // Set the y label and its range
  plot.ylabel("liczba wystąpień");
  plot.yrange(0., 20);
  plot.xlabel("czas [s]");
  plot.xtics().fontSize(12);

  drawDistribution(plot, data, 0.01);

  Figure fig = {{plot}};
  Canvas canvas = {{fig}};

  canvas.size(1000, 500);
  //  canvas.show();

  auto trimmed = dataset_name.substr(0, dataset_name.find('.'));

  auto path = fmt::format("/home/user/pracownia_wstepna/wahadlo/images/dystrybucja_{}.png", trimmed);

  canvas.save(path);
}

constexpr ld kAun = 0.01;
constexpr ld kDun = 0.01;
constexpr ld kHun = 0.002;
constexpr ld kTun = 0.51;

auto handleMeasuredPeriods() -> void {
  using MeasurementTuple = std::tuple<std::size_t, std::size_t, Measurement>;

  std::vector<MeasurementTuple> periods_measurements = {
      {1, 1, {}},
      {2, 1, {}},
      {5, 1, {}},
      {10, 1, {}},
      {10, 2, {}},
      {10, 3, {}},
      {10, 4, {}},
      {10, 5, {}},
  };
  auto generate_name = [](std::size_t period_multiple, std::size_t height_id) {
    return fmt::format("x{}h{}.txt", period_multiple, height_id);
  };

  for (auto& [period_multiple, height_id, period_measurement] : periods_measurements) {
    auto name = generate_name(period_multiple, height_id);
    auto dataset = readDatasset(name);
    period_measurement = setupMeasurement(dataset, kTun);
  }

  for (auto const& [period_multiple, height_id, period_measurement] : periods_measurements) {
    auto name = generate_name(period_multiple, height_id);

    auto period_quantity = static_cast<Quantity>(period_measurement);

    auto transform_to_x1 = [n = period_multiple](ld t_x_n) -> ld { return t_x_n / n; };

    auto x1_transformed = combineQuantities(transform_to_x1, period_quantity);

    fmt::print("{}: {{ t_{}: {{ {} }}, t_1: {{ {} }} }}\n", name, period_multiple, period_measurement, x1_transformed);
  }
}

auto calculateExpectedPeriods() -> void {
  constexpr auto kG = 9.80665;

  auto a = static_cast<Quantity>(setupMeasurement(readDatasset("a.txt"), kAun));
  auto d = static_cast<Quantity>(setupMeasurement(readDatasset("d.txt"), kDun));

  using MeasurementTuple = std::tuple<std::size_t, Measurement>;

  std::vector<MeasurementTuple> h_measurements = {
      {1, {}},
      {2, {}},
      {3, {}},
      {4, {}},
      {5, {}},
  };

  auto generate_name = [](std::size_t height_id) { return fmt::format("h{}.txt", height_id); };

  for (auto& [height_id, measurement] : h_measurements) {
    auto name = generate_name(height_id);
    auto dataset = readDatasset(name);
    measurement = setupMeasurement(dataset, kHun);
  }

  for (auto const& [height_id, measurement] : h_measurements) {
    auto h = static_cast<Quantity>(measurement);
    auto l = combineQuantities([](ld a, ld h, ld d) { return a - h - d; }, a, h, d);

    auto t = combineQuantities([](ld l) { return 2 * std::numbers::pi * std::sqrt(l / kG); }, l);

    fmt::print("l: {}: {{ {} }}\n", l.value_, t);
  }
}

auto calculateGravity() -> void {
  struct MeasurementMetadata {
    std::size_t period_multiple_;
    std::size_t height_id_;
    Measurement period_measurement_;
    Measurement h_measurement_;
  };

  std::vector<MeasurementMetadata> periods_measurements = {
      {1, 1},
      {2, 1},
      {5, 1},
      {10, 1},
      {10, 2},
      {10, 3},
      {10, 4},
      {10, 5},
  };
  auto generate_name_period = [](std::size_t period_multiple, std::size_t height_id) {
    return fmt::format("x{}h{}.txt", period_multiple, height_id);
  };

  auto generate_name_h = [](std::size_t height_id) { return fmt::format("h{}.txt", height_id); };

  auto a = static_cast<Quantity>(setupMeasurement(readDatasset("a.txt"), kAun));
  auto d = static_cast<Quantity>(setupMeasurement(readDatasset("d.txt"), kDun));

  for (auto& [period_multiple, height_id, period_measurement, h_measurement] : periods_measurements) {
    auto period_dataset = readDatasset(generate_name_period(period_multiple, height_id));
    auto h_dataset = readDatasset(generate_name_h(height_id));

    period_measurement = setupMeasurement(period_dataset, kTun);
    h_measurement = setupMeasurement(h_dataset, kHun);
  }

  for (auto const& [period_multiple, height_id, period_measurement, h_measurement] : periods_measurements) {
    auto period_quantity = static_cast<Quantity>(period_measurement);

    auto transform_to_x1 = [n = period_multiple](ld t_x_n) -> ld { return t_x_n / n; };

    auto tx1 = combineQuantities(transform_to_x1, period_quantity);

    auto h = static_cast<Quantity>(h_measurement);

    auto l = combineQuantities([](ld a, ld h, ld d) { return a - h - d; }, a, h, d);

    auto g =
        combineQuantities([](ld l, ld t) { return 4 * std::numbers::pi * std::numbers::pi * l / (t * t); }, l, tx1);

    fmt::print("g_{}: {{ {} }}\n", height_id, g);
  }
}

auto calculateLengths() -> void {
  //  auto a = static_cast<Quantity>(setupMeasurement(readDataset("a.txt"), kAun));
  //  auto d = static_cast<Quantity>(setupMeasurement(readDataset("d.txt"), kDun));
}

auto main() -> int {  // NOLINT
  calculateGravity();
}
