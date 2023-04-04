#include "uncertainties/uncertainties.hpp"

auto main() -> int {  // NOLINT
  using uncertainty::combineMeasurements;
  using uncertainty::Measurement;
  using uncertainty::partialDerivative;
  using uncertainty::Quantity;
  using uncertainty::setupMeasurement;
  using utils::ld;

  auto a = static_cast<Quantity>(setupMeasurement(utils::readDataset("a.txt"), 0.01));
  auto d = static_cast<Quantity>(setupMeasurement(utils::readDataset("d.txt"), 0.00002));

  auto h1 = static_cast<Quantity>(setupMeasurement(utils::readDataset("h1.txt"), 0.01));
  auto l1 = combineMeasurements([](ld a, ld h, ld d) { return a - h - d; }, a, h1, d);
  auto tx1h1 = static_cast<Quantity>(setupMeasurement(utils::readDataset("x1h1.txt"), 0.01));
  auto tx2h1 = static_cast<Quantity>(setupMeasurement(utils::readDataset("x2h1.txt"), 0.01));
  auto tx5h1 = static_cast<Quantity>(setupMeasurement(utils::readDataset("x5h1.txt"), 0.01));
  auto tx10h1 = static_cast<Quantity>(setupMeasurement(utils::readDataset("x10h1.txt"), 0.01));

  auto h2 = static_cast<Quantity>(setupMeasurement(utils::readDataset("h2.txt"), 0.01));
  auto l2 = combineMeasurements([](ld a, ld h, ld d) { return a - h - d; }, a, h2, d);
  auto tx10h2 = static_cast<Quantity>(setupMeasurement(utils::readDataset("x10h2.txt"), 0.01));

  auto h3 = static_cast<Quantity>(setupMeasurement(utils::readDataset("h3.txt"), 0.01));
  auto l3 = combineMeasurements([](ld a, ld h, ld d) { return a - h - d; }, a, h3, d);
  auto tx10h3 = static_cast<Quantity>(setupMeasurement(utils::readDataset("x10h3.txt"), 0.01));

  auto h4 = static_cast<Quantity>(setupMeasurement(utils::readDataset("h4.txt"), 0.01));
  auto l4 = combineMeasurements([](ld a, ld h, ld d) { return a - h - d; }, a, h4, d);
  auto tx10h4 = static_cast<Quantity>(setupMeasurement(utils::readDataset("x10h4.txt"), 0.01));

  auto h5 = static_cast<Quantity>(setupMeasurement(utils::readDataset("h5.txt"), 0.01));
  auto l5 = combineMeasurements([](ld a, ld h, ld d) { return a - h - d; }, a, h5, d);
  auto tx10h5 = static_cast<Quantity>(setupMeasurement(utils::readDataset("x10h5.txt"), 0.01));

  //  fmt::print("h1: {}\n", h1);
  //  fmt::print("tx10h1: {}\n\n", tx10h1);
  //
  //  fmt::print("h2: {}\n", h2);
  //  fmt::print("tx10h2: {}\n\n", tx10h2);
  //
  //  fmt::print("h3: {}\n", h3);
  //  fmt::print("tx10h3: {}\n\n", tx10h3);
  //
  //  fmt::print("h4: {}\n", h4);
  //  fmt::print("tx10h4: {}\n\n", tx10h4);
  //
  //  fmt::print("h5: {}\n", h5);
  //  fmt::print("tx10h5: {}\n\n", tx10h5);

  fmt::print("h1: {}\n", h1);
  fmt::print("l1: {}\n\n", l1);

  fmt::print("h2: {}\n", h2);
  fmt::print("l2: {}\n\n", l2);

  fmt::print("h3: {}\n", h3);
  fmt::print("l3: {}\n\n", l3);

  fmt::print("h4: {}\n", h4);
  fmt::print("l4: {}\n\n", l4);

  fmt::print("h5: {}\n", h5);
  fmt::print("l5: {}\n\n", l5);
}
