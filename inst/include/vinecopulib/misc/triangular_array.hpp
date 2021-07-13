// Copyright © 2016-2021 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <iostream>
#include <stdexcept>
#include <vector>

namespace vinecopulib {

//! @brief Triangular arrays.
//!
//! A triangular array behaves like a matrix with the structure
//! ```
//! x x x x x
//! x x x x
//! x x x
//! x x
//! x
//! ```
//! and all other elements omitted. This structure appears naturally in the
//! representation of a vine copula model and related algorithms. Each row
//! corresponds to one tree in the vine, starting from the top. In each tree
//! (=row), each column represents an edge.
//!
//! For truncated vine models the last few trees are omitted. For example, a
//! 3-truncated version of the above array contains the elements
//! ```
//! x x x x x
//! x x x x
//! x x x
//! ```
//! Only the elements indicated by `x`s are stored and can be accessed.
//!
//! The data structure is templated and any type or class can be used to fill
//! the entries (`x`s) of the triangular array.
template<typename T>
class TriangularArray
{
public:
  TriangularArray() = default;
  explicit TriangularArray(size_t d);
  TriangularArray(size_t d, size_t trunc_lvl);
  explicit TriangularArray(const std::vector<std::vector<T>>& rows);

  T& operator()(size_t row, size_t column);
  T operator()(size_t row, size_t column) const;
  bool operator==(const TriangularArray<T>& rhs) const;
  void truncate(size_t trunc_lvl);

  size_t get_trunc_lvl() const;
  size_t get_dim() const;

  std::string str() const;

private:
  size_t d_{ 0 };
  size_t trunc_lvl_{ 0 };
  std::vector<std::vector<T>> arr_;
};

//! @brief Construct a triangular array of dimension `d`.
//!
//! The array has `d-1` columns and `d-1` rows.
//! @param d The dimension of the underlying vine.
template<typename T>
TriangularArray<T>::TriangularArray(size_t d)
  : TriangularArray(d, d - 1)
{}

//! @brief Construct a truncated triangular array.
//!
//! The array has `d-1` columns and `min(trunc_lvl, d-1)` rows.
//! @param d The dimension of the vine.
//! @param trunc_lvl The truncation level.
template<typename T>
TriangularArray<T>::TriangularArray(size_t d, size_t trunc_lvl)
  : d_(d)
  , trunc_lvl_(std::min(d - 1, trunc_lvl))
{
  if (d < 1)
    throw std::runtime_error("d should be greater than 0");

  arr_ = std::vector<std::vector<T>>(trunc_lvl_);
  for (size_t i = 0; i < trunc_lvl_; i++)
    arr_[i] = std::vector<T>(d_ - i);
}

//! @brief Construct a truncated triangular array from nested vector.
//!
//! An arrax of dimension `d` has `d-1` columns and `min(trunc_lvl, d-1)`
//! rows.
//! @param rows A vector of rows; the length of the first row defines
//! the dimension of the triangular array. The number of rows defines
//! the truncation level.
template<typename T>
TriangularArray<T>::TriangularArray(const std::vector<std::vector<T>>& rows)
  : trunc_lvl_(rows.size())
{
  if (trunc_lvl_ == 0) {
    return;
  } else {
    d_ = rows[0].size() + 1;
  }
  if (trunc_lvl_ > d_) {
    throw std::runtime_error("Not a triangular array: more rows than columns.");
  }
  for (size_t i = 0; i < trunc_lvl_; i++) {
    if (rows[i].size() != d_ - 1 - i) {
      throw std::runtime_error(
        "Not a triangular array: row i must have (d - 1 - i) entries.");
    }
  }

  arr_ = rows;
}

//! @brief Access one element of the trapezoid (writable).
//! @param row The row level.
//! @param column The column in this row.
template<typename T>
T&
TriangularArray<T>::operator()(size_t row, size_t column)
{
  assert(row < trunc_lvl_);
  assert(column < d_ - 1 - row);
  return arr_[row][column];
}

//! @brief Access one element of the trapezoid (non-writable).
//! @param row The row level.
//! @param column The column in this row.
template<typename T>
T
TriangularArray<T>::operator()(size_t row, size_t column) const
{
  assert(row < trunc_lvl_);
  assert(column < d_ - 1 - row);
  return arr_[row][column];
}

//! @brief Truncates the trapezoid.
//! If the trapezoid is already truncated at a level
//! less than `trunc_lvl`, the function does nothing.
//! @param trunc_lvl The truncation level.
template<typename T>
void
TriangularArray<T>::truncate(size_t trunc_lvl)
{
  if (trunc_lvl < trunc_lvl_) {
    trunc_lvl_ = trunc_lvl;
    arr_.resize(trunc_lvl);
  }
}

//! @brief Equality operator to compare two TriangularArray objects.
//! @param rhs Right-hand-side of the equality operator.
template<typename T>
bool
TriangularArray<T>::operator==(const TriangularArray<T>& rhs) const
{
  if ((d_ != rhs.get_dim()) | (trunc_lvl_ != rhs.get_trunc_lvl()))
    return false;

  for (size_t i = 0; i < trunc_lvl_; i++) {
    for (size_t j = 0; j < d_ - 1 - i; j++) {
      if ((*this)(i, j) != rhs(i, j))
        return false;
    }
  }
  return true;
}

//! Get the truncation level of the underlying vine..
template<typename T>
size_t
TriangularArray<T>::get_trunc_lvl() const
{
  return trunc_lvl_;
}

//! Get the dimension of the underlying vine (the matrix has `d-1` columns and.
//! `min(trunv_lvl, d-1)` rows).
template<typename T>
size_t
TriangularArray<T>::get_dim() const
{
  return d_;
}

//! represent triangular array as a string.
template<typename T>
std::string
TriangularArray<T>::str() const
{
  std::stringstream str;
  for (size_t i = 0; i < trunc_lvl_; i++) {
    for (size_t j = 0; j < d_ - 1 - i; j++) {
      str << arr_[i][j] << " ";
    }
    str << std::endl;
  }
  return str.str();
}

//! @brief Ostream method for TriangularArray, to be used with `std::cout`.
//! @param os An output stream.
//! @param tri_array A triangular array.
template<typename T>
std::ostream&
operator<<(std::ostream& os, const TriangularArray<T>& tri_array)
{
  os << tri_array.str();
  return os;
}

} // end of namespace vinecopulib!
