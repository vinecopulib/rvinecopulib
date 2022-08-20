// Copyright © 2016-2021 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/misc/tools_serialization.hpp>
#include <vinecopulib/misc/tools_stats.hpp>
#include <vinecopulib/misc/tools_stl.hpp>

namespace vinecopulib {

//! @brief Instantiates an RVineStructure object from a matrix representing an
//! R-vine array.
//!
//! The matrix must contain zeros in the lower right triangle and
//! the upper left triangle must be a valid R-vine array. Truncated vines can
//! be encoded by putting zeros above the digonal in all rows below the
//! truncation level. Example of a 1-truncated matrix:
//! ```
//! 4 4 4 4
//! 0 0 3 0
//! 0 2 0 0
//! 1 0 0 0
//! ```
//! @param mat A matrix representing a valid R-vine array.
//! @param check Whether `mat` shall be checked for validity.
inline RVineStructure::RVineStructure(
  const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& mat,
  bool check)
{
  d_ = mat.cols();
  if (check) {
    check_if_quadratic(mat);
    check_lower_tri(mat);
  }

  order_ = get_order(mat);
  if (check)
    check_antidiagonal();

  trunc_lvl_ = find_trunc_lvl(mat);
  struct_array_ = to_rvine_array(mat);

  if (check)
    check_upper_tri();

  struct_array_ = to_natural_order();
  if (check)
    check_columns();

  min_array_ = compute_min_array();
  if (check)
    check_proximity_condition();

  needed_hfunc1_ = compute_needed_hfunc1();
  needed_hfunc2_ = compute_needed_hfunc2();
}

//! @brief Instantiates an RVineStructure object to a D-vine for a given
//! dimension.
//! @param d The dimension.
//! @param trunc_lvl The truncation level. By default, it is dim - 1.
inline RVineStructure::RVineStructure(const size_t& d, const size_t& trunc_lvl)
  : RVineStructure(tools_stl::seq_int(1, d), std::min(d - 1, trunc_lvl), false)
{}

//! @brief Instantiates an RVineStructure object to a D-vine with a given
//! ordering of the variables.
//! @param order The order of variables in the D-vine (diagonal entries in the
//!    R-vine array); must be a permutation of 1, ..., d.
//! @param trunc_lvl The truncation level. By default, it is d - 1.
//! @param check Whether `order shall be checked for validity.
inline RVineStructure::RVineStructure(const std::vector<size_t>& order,
                                      const size_t& trunc_lvl,
                                      bool check)
  : RVineStructure(
      order,
      make_dvine_struct_array(order.size(),
                              std::min(order.size() - 1, trunc_lvl)),
      true,
      false)
{
  if (check)
    check_antidiagonal();
}

//! @brief Instantiates an RVineStructure object from the variable order
//! (diagonal elements of the R-vine array) and a triangular structure array
//! (all elements above the diagonal).
//!
//! @param order The order of variables (diagonal entries in the
//!    R-vine array); must be a permutation of 1, ..., d.
//! @param struct_array The structure array  (all elements
//!    above the diagonal in the R-vine array). For truncated vines, all rows
//!    below the truncation level are omitted.
//! @param natural_order Whether `struct_array` is already in natural order.
//! @param check Whether `order` and `struct_array` shall be checked for
//! validity.
inline RVineStructure::RVineStructure(
  const std::vector<size_t>& order,
  const TriangularArray<size_t>& struct_array,
  bool natural_order,
  bool check)
  : order_(order)
  , d_(order.size())
  , trunc_lvl_(struct_array.get_trunc_lvl())
  , struct_array_(struct_array)
{
  if (check) {
    if ((trunc_lvl_ > 0) & (struct_array.get_dim() != d_)) {
      throw std::runtime_error("order and struct_array have "
                               "incompatible dimensions");
    }
    check_antidiagonal();
  }

  if (trunc_lvl_ > 0) {
    if (check)
      check_upper_tri();
    if (!natural_order)
      struct_array_ = to_natural_order();
    if (check)
      check_columns();

    min_array_ = compute_min_array();
    if (check)
      check_proximity_condition();

    needed_hfunc1_ = compute_needed_hfunc1();
    needed_hfunc2_ = compute_needed_hfunc2();
  } else {
    struct_array_ = TriangularArray<size_t>(d_, trunc_lvl_);
    min_array_ = TriangularArray<size_t>(d_, trunc_lvl_);
    needed_hfunc1_ = TriangularArray<short unsigned>(d_, trunc_lvl_);
    needed_hfunc2_ = TriangularArray<short unsigned>(d_, trunc_lvl_);
  }
}

//! @brief Instantiates from a nlohmann::json object.
//! @param input The nlohmann::json object to convert from
//! (see `to_json()` for the structure of the input).
//! @param check Whether to check if the input represents
//!      a valid R-vine structure.
inline RVineStructure::RVineStructure(const nlohmann::json& input,
                                      const bool check)
  : RVineStructure(
      tools_serialization::json_to_vector<size_t>(input["order"]),
      tools_serialization::json_to_triangular_array<size_t>(input["array"]),
      check)
{}

//! @brief Instantiates an RVineStructure from a JSON file.
//!
//! The file needs to contain two values: `"array"` for the structure
//! triangular array and `"order"` for the order vector.
//!
//! @param filename The name of the JSON file to read.
//! @param check Whether to check if the input represents
//!      a valid R-vine matrix.
inline RVineStructure::RVineStructure(const std::string& filename,
                                      const bool check)
  : RVineStructure(tools_serialization::file_to_json(filename), check)
{}

//! @brief Converts the structure into a nlohmann::json object.
//!
//! The `nlohmann::json` object contains two nodes: `"array"` for the structure
//! triangular array and `"order"` for the order vector.
//!
//! @return the nlohmann::json object containing the structure.
inline nlohmann::json
RVineStructure::to_json() const
{
  nlohmann::json output;
  auto array_json =
    tools_serialization::triangular_array_to_json(struct_array_);
  output["array"] = array_json;
  auto order_json = tools_serialization::vector_to_json(order_);
  output["order"] = order_json;

  return output;
}

//! @brief Write the structure into a JSON file.
//!
//! The written file contains two values: `"array"` for the structure
//! triangular array and `"order"` for the order vector.
//!
//! @param filename The name of the file to write.
inline void
RVineStructure::to_file(const std::string& filename) const
{
  tools_serialization::json_to_file(filename, this->to_json());
}

//! Gets the dimension of the vine.
inline size_t
RVineStructure::get_dim() const
{
  return d_;
}

//! Gets the truncation level of the vine.
inline size_t
RVineStructure::get_trunc_lvl() const
{
  return trunc_lvl_;
}

//! @brief Extract the order of variables in the vine (diagonal entries in the
//! R-vine array).
inline std::vector<size_t>
RVineStructure::get_order() const
{
  return order_;
}

//! @brief Extract structure array (all elements above the diagonal in the
//! R-vine array).
//! @param natural_order Whether indices correspond to natural order.
inline TriangularArray<size_t>
RVineStructure::get_struct_array(bool natural_order) const
{
  if (natural_order) {
    return struct_array_;
  }

  // convert all labels to original order
  auto new_array = struct_array_;
  for (size_t tree = 0; tree < trunc_lvl_; tree++) {
    for (size_t edge = 0; edge < d_ - 1 - tree; edge++) {
      new_array(tree, edge) = this->struct_array(tree, edge, false);
    }
  }
  return new_array;
}

//! @brief Gets the minimum array.
//!
//! The minimum array is derived from an R-vine array by
//! iteratively computing the (elementwise) minimum of two subsequent rows
//! (starting from the top). It is used in estimation and evaluation
//! algorithms to find the two edges in the previous tree that are joined by
//! the current edge.
inline TriangularArray<size_t>
RVineStructure::get_min_array() const
{
  return min_array_;
}

//! @brief Gets an array indicating which of the first h-functions are needed.
//!
//! It is usually not necessary to compute both h-functions for each
//! pair-copula.
inline TriangularArray<short unsigned>
RVineStructure::get_needed_hfunc1() const
{
  return needed_hfunc1_;
}

//! @brief Gets an array indicating which of the second h-functions are needed.
//!
//! It is usually not necessary to compute both h-functions for each
//! pair-copula.
inline TriangularArray<short unsigned>
RVineStructure::get_needed_hfunc2() const
{
  return needed_hfunc2_;
}

//! @brief Accesses elements of the structure array.
//! @param tree Tree index.
//! @param edge Edge index.
//! @param natural_order Whether indices correspond to natural order.
inline size_t
RVineStructure::struct_array(size_t tree, size_t edge, bool natural_order) const
{
  if (natural_order) {
    return struct_array_(tree, edge);
  }
  // convert label back to original order
  return order_[struct_array_(tree, edge) - 1];
}

//! @brief Access elements of the minimum array.
//! @param tree Tree index.
//! @param edge Edge index.
inline size_t
RVineStructure::min_array(size_t tree, size_t edge) const
{
  return min_array_(tree, edge);
}

//! @brief Access elements of the needed_hfunc1 array.
//! @param tree Tree index.
//! @param edge Edge index.
inline bool
RVineStructure::needed_hfunc1(size_t tree, size_t edge) const
{
  return needed_hfunc1_(tree, edge);
}

//! @brief Access elements of the needed_hfunc2 array.
inline bool
RVineStructure::needed_hfunc2(size_t tree, size_t edge) const
{
  return needed_hfunc2_(tree, edge);
}

//! @brief Truncates the R-vine structure.
//! @param trunc_lvl The truncation level.
//!
//! If the structure is already truncated at a level
//! less than `trunc_lvl`, the function does nothing.
inline void
RVineStructure::truncate(size_t trunc_lvl)
{
  if (trunc_lvl < trunc_lvl_) {
    struct_array_.truncate(trunc_lvl);
    min_array_.truncate(trunc_lvl);
    needed_hfunc1_.truncate(trunc_lvl);
    needed_hfunc2_.truncate(trunc_lvl);
    trunc_lvl_ = struct_array_.get_trunc_lvl();
  }
}

//! @brief Converts the structure to a string representation (most useful for
//! printing).
inline std::string
RVineStructure::str() const
{
  std::stringstream str;
  for (size_t i = 0; i < d_ - 1; i++) {
    for (size_t j = 0; j < d_ - i - 1; j++) {
      if (i < trunc_lvl_) {
        str << this->struct_array(i, j, false) << " ";
      } else {
        str << "  ";
      }
    }
    str << order_[d_ - 1 - i] << " " << std::endl;
  }
  str << order_[0] << " " << std::endl;

  return str.str();
}

//! @brief Randomly sample a regular vine structure.
//! @param d The dimension.
//! @param natural_order Should the sampled structure be in natural order?
//! @param seeds Seeds of the random number generator; if empty (default),
//!   the random number generator is seeded randomly.
//! @note Implementation of Algorithm 13 in Harry Joe's 2014 book (p. 288),
//! but there's a typo: the end of line 6 in the book should be
//! 'column j' instead of 'column k'.
inline RVineStructure
RVineStructure::simulate(size_t d, bool natural_order, std::vector<int> seeds)
{
  auto U = tools_stats::simulate_uniform(d, d, false, seeds);

  // A is the R-vine matrix we want to create (upper right-triag format).
  // B is a random binary representation that we need to convert.
  Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> A(d, d), B(d, d);
  A.setZero();
  B = (U.leftCols(d).array() > 0.5).cast<size_t>();

  for (size_t i = 0; i < d; i++) {
    A(i, i) = i + 1;
    B(i, i) = 1;
    if (i > 0) {
      A(i - 1, i) = i;
      B(0, i) = 1;
      B(i - 1, i) = 1;
    }
  }
  if (d > 2) {
    A(0, 2) = 1;
  }

  for (size_t j = 3; j < d; j++) {
    size_t ac = j - 2;
    auto to_assign = tools_stl::seq_int(1, j - 1);
    for (ptrdiff_t k = j - 2; k >= 0; k--) {
      if (B(k, j) == 1) {
        A(k, j) = ac + 1;
        to_assign = tools_stl::set_diff(to_assign, { A(k, j) });
        if (k > 0) {
          // to_assign is always ordered ascendingly -> we pick largest
          ac = to_assign[to_assign.size() - 1] - 1;
        }
      } else {
        A(k, j) = A(k - 1, ac);
        to_assign = tools_stl::set_diff(to_assign, { A(k, j) });
      }
    }
  }

  // need to convert to upper left triangular form (our notation)
  auto rvm =
    RVineStructure(Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>(
      A.rowwise().reverse()));

  // sampling the variable order randomly the first column of U has not been
  // used to construct B, hence it is stochastically independent of B. Calling
  // get_order gives us a permutation of (1, ..., d) that is independent of B.
  std::vector<size_t> order(d);
  if (!natural_order) {
    std::vector<double> u(U.data(), U.data() + d); // first column of U
    auto osim = tools_stl::get_order(u);
    for (size_t k = 0; k < d; k++) {
      order[k] = osim[k] + 1;
    }
  } else {
    for (size_t i = 0; i < d; ++i) {
      order[i] = i + 1;
    }
  }

  return RVineStructure(order, rvm.get_struct_array(true), true, false);
}

//! @brief Gets the R-vine matrix representation.
inline Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>
RVineStructure::get_matrix() const
{
  Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> mat(d_, d_);
  mat.fill(0);
  for (size_t i = 0; i < trunc_lvl_; ++i) {
    for (size_t j = 0; j < d_ - i - 1; ++j) {
      mat(i, j) = this->struct_array(i, j, false);
    }
  }
  for (size_t i = 0; i < d_; ++i) {
    mat(d_ - i - 1, i) = order_[i];
  }
  return mat;
}

//! @brief Find the truncation level in an R-vine array.
//!
//! The truncation level is
//! determined by the first row (starting from the bottom) that contains only
//! zeros above the diagonal.
//!
//! @param mat An array representing the R-vine array.
inline size_t
RVineStructure::find_trunc_lvl(
  const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& mat) const
{
  size_t trunc_lvl;
  size_t d = mat.cols();

  std::stringstream problem;
  problem << "not a valid R-vine array: "
          << "a row with a 0 above the diagonal contains one or more "
          << "non-zero values.";

  for (trunc_lvl = d - 1; trunc_lvl > 0; trunc_lvl--) {
    std::vector<size_t> row_vec(d - trunc_lvl);
    Eigen::Matrix<size_t, Eigen::Dynamic, 1>::Map(&row_vec[0], d - trunc_lvl) =
      mat.row(trunc_lvl - 1).head(d - trunc_lvl);

    if (*(std::min_element(row_vec.begin(), row_vec.end())) != 0)
      break;
  }

  return trunc_lvl;
}

//! @brief Find the order of an R-vine array.
//!
//! The order is contained in the counter-diagonal of the R-vine array.
//! @param mat A matrix representing the R-vine array.
inline std::vector<size_t>
RVineStructure::get_order(
  const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& mat) const
{
  std::vector<size_t> order(d_);
  for (size_t i = 0; i < d_; i++)
    order[i] = mat(d_ - i - 1, i);

  return order;
}

//! @brief Gets the structure array (entries above the diagonal in R-vine.
//! array).
//! @param mat A array representing the R-vine array.
inline TriangularArray<size_t>
RVineStructure::to_rvine_array(
  const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& mat) const
{
  // copy upper triangle
  TriangularArray<size_t> struct_array(d_, trunc_lvl_);
  for (size_t j = 0; j < d_ - 1; j++) {
    for (size_t i = 0; i < std::min(d_ - 1 - j, trunc_lvl_); i++) {
      struct_array(i, j) = mat(i, j);
    }
  }

  return struct_array;
}

//! @brief Converts `struct_array_` to natural order.
inline TriangularArray<size_t>
RVineStructure::to_natural_order() const
{
  // create vector of new variable labels
  auto order = tools_stl::get_order(get_order());

  // relabel to natural order
  TriangularArray<size_t> struct_array(d_, trunc_lvl_);
  for (size_t j = 0; j < d_ - 1; j++) {
    for (size_t i = 0; i < std::min(d_ - 1 - j, trunc_lvl_); i++) {
      struct_array(i, j) = order[struct_array_(i, j) - 1] + 1;
    }
  }

  return struct_array;
}

//! @brief Creates a structure array corresponding to a D-vine in natural order.
inline TriangularArray<size_t>
RVineStructure::make_dvine_struct_array(size_t d, size_t trunc_lvl)
{
  TriangularArray<size_t> struct_array(d, trunc_lvl);
  for (size_t j = 0; j < d - 1; j++) {
    for (size_t i = 0; i < std::min(d - 1 - j, trunc_lvl); i++) {
      struct_array(i, j) = i + j + 2;
    }
  }

  return struct_array;
}

//! @brief Creates a structure array corresponding to a D-vine in natural order.
inline TriangularArray<size_t>
RVineStructure::make_cvine_struct_array(size_t d, size_t trunc_lvl)
{
  TriangularArray<size_t> struct_array(d, trunc_lvl);
  for (size_t i = 0; i < std::min(d - 1, trunc_lvl); i++) {
    for (size_t j = 0; j < d - 1 - i; j++) {
      struct_array(i, j) = d - i;
    }
  }

  return struct_array;
}

inline TriangularArray<size_t>
RVineStructure::compute_min_array() const
{
  TriangularArray<size_t> min_array = struct_array_;
  for (size_t j = 0; j < d_ - 1; j++) {
    for (size_t i = 1; i < std::min(d_ - 1 - j, trunc_lvl_); i++) {
      min_array(i, j) = std::min(struct_array_(i, j), min_array(i - 1, j));
    }
  }

  return min_array;
}

inline TriangularArray<short unsigned>
RVineStructure::compute_needed_hfunc1() const
{
  TriangularArray<short unsigned> needed_hfunc1(d_, trunc_lvl_);
  if (d_ == 1) {
    return needed_hfunc1;
  }

  for (size_t i = 0; i < std::min(d_ - 2, trunc_lvl_ - 1); i++) {
    for (size_t j = 0; j < d_ - 2 - i; j++) {
      if (struct_array_(i + 1, j) != min_array_(i + 1, j))
        needed_hfunc1(i, min_array_(i + 1, j) - 1) = 1;
    }
  }

  return needed_hfunc1;
}

inline TriangularArray<short unsigned>
RVineStructure::compute_needed_hfunc2() const
{
  TriangularArray<short unsigned> needed_hfunc2(d_, trunc_lvl_);
  if (d_ == 1) {
    return needed_hfunc2;
  }

  for (size_t i = 0; i < std::min(d_ - 2, trunc_lvl_ - 1); i++) {
    for (size_t j = 0; j < d_ - 2 - i; j++) {
      needed_hfunc2(i, j) = 1;
      if (struct_array_(i + 1, j) == min_array_(i + 1, j))
        needed_hfunc2(i, min_array_(i + 1, j) - 1) = 1;
    }
  }

  return needed_hfunc2;
}

inline void
RVineStructure::check_if_quadratic(
  const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& mat) const
{
  std::string problem = "must be quadratic.";
  if (mat.rows() != mat.cols()) {
    throw std::runtime_error("not a valid R-vine array: " + problem);
  }
}

inline void
RVineStructure::check_lower_tri(
  const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& mat) const
{
  std::string problem = "the lower right triangle must only contain zeros";
  size_t sum_lwr = 0;
  for (size_t j = 1; j < d_; ++j) {
    sum_lwr += mat.block(d_ - j, j, j, 1).array().sum();
    if (sum_lwr != 0) {
      throw std::runtime_error("not a valid R-vine array: " + problem);
    }
  }
}

inline void
RVineStructure::check_upper_tri() const
{
  std::string problem;
  problem += "the upper left triangle can only contain numbers ";
  problem += "between 1 and d (number of variables).";

  for (size_t i = 0; i < trunc_lvl_; ++i) {
    for (size_t j = 0; j < d_ - 1 - i; ++j) {
      if ((struct_array_(i, j) < 1) | (struct_array_(i, j) > d_)) {
        throw std::runtime_error("not a valid R-vine array: " + problem);
      }
    }
  }
}

inline void
RVineStructure::check_columns() const
{
  std::string problem = "";
  for (size_t j = 0; j < d_ - 1; j++) {
    // read column into vector so we can use stl methods
    std::vector<size_t> col(std::min(trunc_lvl_, d_ - 1 - j));
    for (size_t i = 0; i < col.size(); i++) {
      col[i] = struct_array_(i, j);
    }

    std::sort(col.begin(), col.end());
    if (col[0] <= 1 + j) {
      problem += "the antidiagonal entry of a column must not be ";
      problem += "contained in any column further to the right.";
    }

    size_t unique_in_col = std::unique(col.begin(), col.end()) - col.begin();
    if (unique_in_col != col.size()) {
      problem = "a column must not contain duplicate entries.";
    }
    if (problem != "") {
      throw std::runtime_error("not a valid R-vine array: " + problem);
    }
  }
}

inline void
RVineStructure::check_antidiagonal() const
{
  std::string problem;
  problem += "the order/antidiagonal must contain the numbers ";
  problem += "1, ..., d (the number of variables)";
  if (!tools_stl::is_same_set(order_, tools_stl::seq_int(1, d_))) {
    throw std::runtime_error("not a valid R-vine array: " + problem);
  }
}

inline void
RVineStructure::check_proximity_condition() const
{
  for (size_t t = 1; t < trunc_lvl_; ++t) {
    for (size_t e = 0; e < d_ - t - 1; ++e) {
      std::vector<size_t> target_set(t + 1), test_set(t + 1);
      // conditioning set
      for (size_t i = 0; i < t; i++) {
        target_set[i] = struct_array_(i, e);
        test_set[i] = struct_array_(i, min_array_(t, e) - 1);
      }

      // non-diagonal conditioned variable
      target_set[t] = struct_array_(t, e);
      // diagonal conditioned variable in other column
      test_set[t] = min_array_(t, e);

      if (!tools_stl::is_same_set(target_set, test_set)) {
        std::stringstream problem;
        problem << "not a valid R-vine array: "
                << "proximity condition violated; "
                << "cannot extract conditional distribution (" << target_set[t]
                << " | ";
        for (size_t i = 0; i < t - 1; ++i) {
          problem << order_[target_set[i] - 1] << ", ";
        }
        problem << order_[target_set[t - 1] - 1] << ") from pair-copulas.";
        throw std::runtime_error(problem.str().c_str());
      }
    }
  }
}

//! @brief Ostream method for RVineStructure, to be used with `std::cout`.
//! @param os Output stream.
//! @param rvs R-vine structure array.
inline std::ostream&
operator<<(std::ostream& os, const RVineStructure& rvs)
{
  os << rvs.str();
  return os;
}

//! @param order The order of variables in the D-vine (diagonal entries in the
//!    R-vine array); must be a permutation of 1, ..., d.
inline DVineStructure::DVineStructure(const std::vector<size_t>& order)
  : RVineStructure(order,
                   make_dvine_struct_array(order.size(), order.size() - 1),
                   true,
                   false)
{}

//! @param order The order of variables in the D-vine (diagonal entries in the
//!    R-vine array); must be a permutation of 1, ..., d.
//! @param trunc_lvl The truncation level.
inline DVineStructure::DVineStructure(const std::vector<size_t>& order,
                                      size_t trunc_lvl)
  : RVineStructure(order,
                   make_dvine_struct_array(order.size(), trunc_lvl),
                   true,
                   false)
{}

//! @param order The order of variables in the C-vine (diagonal entries in the
//!    R-vine array); must be a permutation of 1, ..., d.
inline CVineStructure::CVineStructure(const std::vector<size_t>& order)
  : RVineStructure(order,
                   make_cvine_struct_array(order.size(), order.size() - 1),
                   true,
                   false)
{}

//! @param order The order of variables in the C-vine (diagonal entries in the
//!    R-vine array); must be a permutation of 1, ..., d.
//! @param trunc_lvl The truncation level.
inline CVineStructure::CVineStructure(const std::vector<size_t>& order,
                                      size_t trunc_lvl)
  : RVineStructure(order,
                   make_cvine_struct_array(order.size(), trunc_lvl),
                   true,
                   false)
{}
}
