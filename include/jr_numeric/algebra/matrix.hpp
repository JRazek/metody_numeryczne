#pragma once

#include <algorithm>
#include <array>
#include <bitset>
#include <cmath>
#include <concepts>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <iterator>
#include <numeric>

#include "jr_numeric/utils/concepts.hpp"

namespace jr_numeric::algebra {

using concepts::FloatingPoint;

namespace rg = std::ranges;

struct MatrixExtent {
  std::size_t rows_;
  std::size_t cols_;
};

template <std::size_t N, std::size_t M, typename T>
class Matrix {
 public:
  using Row = std::array<T, M>;

 private:
  std::array<Row, N> data_;

 public:
  constexpr Matrix() noexcept : data_({}){};

  constexpr explicit Matrix(std::array<std::array<T, M>, N> const& arr) noexcept : data_(arr) {}

  constexpr explicit Matrix(std::array<T, N * M> const& arr) noexcept {
    for (auto y = 0; y < N; y++) {
      for (auto x = 0; x < M; x++) {
        data_[y][x] = arr[y * M + x];
      }
    }
  }

  auto data() const noexcept -> std::array<std::array<T, M>, N> const& { return data_; }

  auto cdata() const noexcept -> std::array<std::array<T, M>, N> const& { return data_; }

  auto data() noexcept -> std::array<std::array<T, M>, N>& { return data_; }

  template <std::size_t X>
  constexpr auto operator*(Matrix<M, X, T> const& rhs) const noexcept -> Matrix<N, X, T> {
    Matrix<N, X, T> res;
    multiply(*this, rhs, res);
    return res;
  }

  constexpr auto operator*(const T scalar) const noexcept -> Matrix {
    auto res = *this;
    multiply(res, scalar);
    return res;
  }

  constexpr auto operator*=(Matrix<M, M, T> const& rhs) noexcept -> Matrix {
    multiply(*this, rhs, *this);
    return *this;
  }

  constexpr auto operator*=(const T scalar) noexcept -> Matrix& {
    multiply(*this, scalar);
    return *this;
  }

  constexpr auto operator+(Matrix const& rhs) const noexcept -> Matrix {
    auto res = *this;
    add(*this, rhs, res);
    return res;
  }

  constexpr auto operator+=(Matrix const& rhs) noexcept -> Matrix& {
    add(*this, rhs, *this);
    return *this;
  }

  constexpr auto operator-(Matrix const& rhs) const noexcept -> Matrix {
    auto res = *this;
    subtract(*this, rhs, res);
    return res;
  }

  constexpr auto operator-=(Matrix const& rhs) noexcept -> Matrix& {
    subtract(*this, rhs, *this);
    return *this;
  }

  constexpr auto operator[](const std::size_t y) const noexcept -> std::array<T, M> const& { return data_[y]; }

  constexpr auto operator[](const std::size_t y) noexcept -> std::array<T, M>& { return data_[y]; }

  constexpr auto operator[](const std::pair<std::size_t, std::size_t> cell) const noexcept -> T const& {
    return data_[cell.first][cell.second];
  }

  constexpr auto operator[](const std::pair<std::size_t, std::size_t> cell) noexcept -> T& {
    return data_[cell.first][cell.second];
  }

  constexpr auto operator==(Matrix const& m) const noexcept -> bool {
    for (auto y = 0; y < N; y++) {
      for (auto x = 0; x < M; x++) {
        if (data_[y][x] != m[y][x]) return false;
      }
    }
    return true;
  }

  /**
   * @brief takes O(N^2*M) time and O(N^2*M) memory
   *
   */
  constexpr auto determinant() const noexcept
    requires(N == M && N > 3)
  {
    /**
     * @brief copies matrix with erased row and column
     *
     */
    auto make_matrix = [](const auto y_erased, const auto x_erased, Matrix const& m) -> Matrix<N - 1, M - 1, T> {
      Matrix<N - 1, M - 1, T> res;
      for (auto y_dest = 0, y_from = 0; y_dest < N - 1; y_dest++, y_from++) {
        if (y_from == y_erased) y_from++;
        for (auto x_dest = 0, x_from = 0; x_dest < N - 1 && y_from < N; x_dest++, x_from++) {
          if (x_from == x_erased) x_from++;
          res[y_dest][x_dest] = m[y_from][x_from];
        }
      }
      return res;
    };

    auto det = T();
    const auto x = static_cast<std::size_t>(0);  // fixed column for row reducing. Why not
    for (auto y = 0; y < N; y++) {
      det += data_[y][x] * std::pow(-1, y + x) * make_matrix(y, x, *this).determinant();
    }
    return det;
  }

  /**
   * @brief takes O(3*3) time and O(1) memory
   *
   */
  constexpr auto determinant() const noexcept
    requires(N == M && N == 3)
  {
    auto det = T();
    for (auto x = 0; x < N; x++) {
      auto tmp1 = T(1);
      auto tmp2 = T(1);
      for (auto y = 0; y < N; y++) {
        tmp1 *= data_[y][(x + y) % M];
        tmp2 *= data_[y][(N - 1 + x - y) % M];
      }
      det += tmp1 - tmp2;
    }
    return det;
  }

  /**
   * @brief takes O(2*2) time and O(1) memory
   *
   */
  constexpr auto determinant() const noexcept
    requires(N == M && N == 2)
  {
    return data_[0][0] * data_[1][1] - data_[0][1] * data_[1][0];
  }

  constexpr auto determinant() const noexcept
    requires(N == M && N == 1)
  {
    return data_[0][0];
  }

  constexpr auto gaussElimination() noexcept -> Matrix<N, M, T>& { return *this; }

  constexpr auto elementWiseXor(Matrix const& m) noexcept -> Matrix& {
    for (auto x = 0; x < M; ++x) {
      for (auto y = 0; y < N; ++y) {
        data_[y][x] ^= m[y][x];
      }
    }
    return *this;
  }

  constexpr auto transpose() const noexcept -> Matrix<M, N, T> {
    Matrix<M, N, T> res;
    for (auto y = 0; y < N; ++y) {
      for (auto x = 0; x < M; ++x) {
        res[x][y] = data_[y][x];
      }
    }
    return res;
  }

  constexpr auto inverse() noexcept -> Matrix<N, M, T>& {}

 private:
  template <std::size_t X>
  constexpr static inline auto multiply(Matrix const& lhs, Matrix<M, X, T> const& rhs, Matrix<N, X, T>& res) noexcept
      -> void {  // TODO (jrazek) implement Strassens
    for (auto y = 0; y < N; y++) {
      for (auto x = 0; x < X; x++) {
        auto sum = T();
        for (auto i = 0; i < M; i++) sum += lhs[y][i] * rhs[i][x];
        res[y][x] = sum;
      }
    }
  }
  constexpr static inline auto multiply(Matrix& res, const T scalar) noexcept -> void {
    for (auto y = 0; y < N; y++)
      for (auto x = 0; x < M; x++) res[y][x] *= scalar;
  }
  constexpr static inline auto add(Matrix const& lhs, Matrix const& rhs, Matrix& res) noexcept -> void {
    for (auto y = 0; y < N; y++)
      for (auto x = 0; x < M; x++) res[y][x] = lhs[y][x] + rhs[y][x];
  }
  constexpr static inline auto subtract(Matrix const& lhs, Matrix const& rhs, Matrix& res) noexcept -> void {
    for (auto y = 0; y < N; y++)
      for (auto x = 0; x < M; x++) res[y][x] = lhs[y][x] - rhs[y][x];
  }

  template <bool Constant>
  struct Iterator {
    using iterator_category = std::forward_iterator_tag;
    using differece_type = std::ptrdiff_t;
    using value_type = std::array<T, M>;
    using pointer = std::conditional_t<Constant, const value_type*, value_type*>;
    using reference = std::conditional_t<Constant, value_type const&, value_type&>;

    explicit Iterator(pointer ptr) : ptr_(ptr) {}

    auto operator*() -> reference { return *ptr_; };

    auto operator->() -> pointer { return ptr_; };
    auto operator++() -> Iterator& {
      ++ptr_;
      return *this;
    };  // preincrement ++it
    auto operator++(int) -> Iterator {
      auto tmp = *this;
      this->operator++;
      return tmp;
    };  // postincrement it++

    friend auto operator==(Iterator const& rhs, Iterator const& lhs) { return rhs.ptr_ == lhs.ptr_; };

    friend auto operator!=(Iterator const& rhs, Iterator const& lhs) { return rhs.ptr_ != lhs.ptr_; };

   private:
    pointer ptr_;
  };
  using iterator = Iterator<false>;
  using const_iterator = Iterator<true>;

 public:
  constexpr auto begin() noexcept -> iterator { return iterator(&data_.front()); }

  constexpr auto end() noexcept -> iterator { return iterator(&data_.back()); }

  constexpr auto cbegin() const noexcept -> const_iterator { return const_iterator(&data_.front()); }

  constexpr auto cend() const noexcept -> const_iterator { return const_iterator(&data_.back()); }
};

template <std::size_t N, std::size_t M, typename T>
auto operator<<(std::ostream& os, Matrix<N, M, T> const& m) -> std::ostream& {
  os << "{\n";
  for (auto y = 0; y < N; y++) {
    os << "\t{";
    for (auto x = 0; x < M; x++) {
      os << m[y][x];
      if (x + 1 != M) os << ", ";
    }
    os << "},\n";
  }
  os << '}';
  return os;
}

template <std::size_t N, std::size_t M, typename T>
struct ReducedMatrix {
  static constexpr auto kReducedSize = std::min(N, M);
  std::array<T, kReducedSize> data_;  // each element contains a non zero value from each row

  struct Proxy {
    T* v_;
    auto operator[](std::size_t /*j*/) noexcept -> T& { return *v_; }
  };

  auto operator[](std::size_t i) -> Proxy& { return Proxy{&data_[i]}; }
};

template <std::size_t N, std::size_t M, std::size_t K, typename T>
auto operator*(Matrix<N, M, T> const& lhs, ReducedMatrix<M, K, T> const& rhs) noexcept -> Matrix<N, K, T> {
  Matrix<N, K, T> res;
  for (auto i = 0u; i < N; i++) {
    for (auto j = 0u; j < K; j++) {
      res[i][j] = lhs[i][j] * rhs[j][i];
    }
  }
  return res;
}

template <std::size_t N, std::size_t M, std::size_t K, typename T>
auto operator*(ReducedMatrix<N, M, T> const& lhs, ReducedMatrix<M, K, T> const& rhs) noexcept -> Matrix<N, K, T> {
  ReducedMatrix<N, K, T> res;
  for (auto i = 0u; i < N; i++) {
    res[i][i] = lhs[i][i] * rhs[i][i];
  }
  return res;
}

template <std::size_t N, typename T>
auto operator*=(ReducedMatrix<N, N, T>& lhs, ReducedMatrix<N, N, T> const& rhs) noexcept -> Matrix<N, N, T> {
  for (auto i = 0u; i < N; i++) {
    lhs[i][i] *= rhs[i][i];
  }
  return lhs;
}

template <std::size_t N, std::size_t M, typename T>
auto operator+(ReducedMatrix<N, M, T> const& lhs, ReducedMatrix<N, M, T> const& rhs) noexcept
    -> ReducedMatrix<N, M, T> {
  ReducedMatrix<N, M, T> res;
  for (auto i = 0u; i < N; i++) {
    res[i][i] = lhs[i][i] + rhs[i][i];
  }
  return res;
}

namespace implementation {

template <FloatingPoint T>
auto equal(T lhs, T rhs) noexcept -> bool {
  return std::abs(lhs - rhs) < std::numeric_limits<T>::epsilon();
}

// sorts all rows by the value at the given column (descending order)
// all the rows from start_row to N will be sorted
template <std::size_t N, std::size_t M, FloatingPoint T>
constexpr static auto sortRows(Matrix<N, M, T>& mat, std::size_t start_row, std::size_t col) noexcept -> void {
  auto rows = std::span(mat.data()).last(N - start_row);

  rg::sort(rows, [col](auto const& row1, auto const& row2) { return std::abs(row1[col]) > std::abs(row2[col]); });
}

// extracts r1 multiplied by factor from r2 (starting from col)
// if element is almost 0, it will be set to 0
template <std::size_t N, std::size_t M, FloatingPoint T>
constexpr static auto extractRow(
    typename Matrix<N, M, T>::Row& r1, typename Matrix<N, M, T>::Row& r2, T factor, std::size_t col = 0) noexcept
    -> void {
  for (auto j = col; j < M; j++) {
    r1[j] -= factor * r2[j];
    if (equal(r2[j], 0.)) r2[j] = 0;
  }
}

template <std::size_t N, std::size_t M, FloatingPoint T>
constexpr static auto extractRowMatrix(Matrix<N, M, T>& mat, std::size_t row, std::size_t col = 0) noexcept -> void {
  for (auto i = row + 1; i < N; i++) {
    auto denom = mat[row][col];
    if (denom == 0) continue;
    auto factor = mat[i][col] / mat[row][col];
    extractRow<N, M, T>(mat[i], mat[row], factor, col);
  }
}

}  // namespace implementation

/**
 * @param mat matrix with sorted rows in non ascending order.
 */
template <std::size_t N, std::size_t M, FloatingPoint T>
constexpr static auto rowEchelon(Matrix<N, M, T>& mat) noexcept -> void {
  auto j = 0;
  for (auto i = 0u; i < N; i++) {
    implementation::sortRows(mat, i, j);
    while (implementation::equal(mat[i][j], 0.)) {
      j++;
      if (j == M) return;
    }
    implementation::extractRowMatrix(mat, i, j);
  }
  // assume
}

/**
 * @param rowEcholon form matrix
 */
template <std::size_t N, std::size_t M, FloatingPoint T>
constexpr static auto rowReduce(Matrix<N, M, T>& mat) noexcept -> void {
  // TODO(jrazek)
  for (auto i = 0u; i < N; i++) {
    auto col = static_cast<std::int64_t>(M - 1 - i);
    auto row = N - i - 1;

    while (implementation::equal(mat[i][col], 0.) && col >= 0) {
      col--;
    }
  }
}

/**
 * @param mat matrix with sorted rows in non ascending order.
 */
template <std::size_t N, std::size_t M, FloatingPoint T>
constexpr static auto reducedRowEchelon(Matrix<N, M, T>& mat) noexcept -> void {}

template <std::size_t N, std::size_t M, FloatingPoint T>
constexpr static auto multiplyRow(std::array<T, M>& lhs, const T rhs) noexcept -> void {
  for (auto i = 0; i < M; i++) lhs[i] *= rhs;
}

// unoptimized and unstable yet
template <std::size_t N, std::size_t M, FloatingPoint T>
auto gaussElimination(Matrix<N, M, T> matrix) {}

}  // namespace jr_numeric::algebra
