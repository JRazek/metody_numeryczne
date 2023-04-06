#pragma once

#include <fmt/core.h>

#include <concepts>
#include <fstream>
#include <functional>
#include <string>

#include "jr_numeric/utils/concepts.hpp"

namespace jr_numeric::utils {

/**
 * @brief Reads a dataset from a file line by line
 *
 * @tparam T Type of the dataset
 * @param path Path to the file
 * @return std::vector<T> Dataset
 */
template <typename T>
auto readDataset(std::string const& path) -> std::vector<T> {
  std::fstream file(path);
  if (file.fail()) {
    throw std::runtime_error("Failed to open file");
  }
  std::vector<T> res;
  T x;
  while (file >> x) {
    res.push_back(x);
  }
  return res;
}

}  // namespace jr_numeric::utils
