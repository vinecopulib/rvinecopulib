// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <algorithm>
#include <fstream>
#include <iostream>
#include <utility>
#include <vinecopulib/misc/tools_stl.hpp>

namespace vinecopulib {

namespace tools_eigen {

//! remove rows of a matrix which contain nan values
//! @param x The matrix.
inline void
remove_nans(Eigen::MatrixXd& x)
{
  // if a row has nan, move it to the end
  size_t last = x.rows() - 1;
  for (size_t i = 0; i < last + 1; i++) {
    if (x.row(i).array().isNaN().any())
      x.row(i--).swap(x.row(last--));
  }
  // remove nan rows
  x.conservativeResize(last + 1, x.cols());
}

//! remove rows of a matrix which contain nan values or have zero weight
//! @param x The matrix.
//! @param a Vector of weights that is either empty or whose size is equal to
//!   the number of columns of x.
inline void
remove_nans(Eigen::MatrixXd& x, Eigen::VectorXd& weights)
{
  if ((weights.size() > 0) && (weights.size() != x.rows()))
    throw std::runtime_error("sizes of x and weights don't match.");

  // if a row has nan or weight is zero, move it to the end
  size_t last = x.rows() - 1;
  for (size_t i = 0; i < last + 1; i++) {
    bool row_has_nan = x.row(i).array().isNaN().any();
    if (weights.size() > 0) {
      row_has_nan = row_has_nan || (std::isnan)(weights(i));
      row_has_nan = row_has_nan || (weights(i) == 0.0);
    }
    if (row_has_nan) {
      if (weights.size() > 0)
        std::swap(weights(i), weights(last));
      x.row(i--).swap(x.row(last--));
    }
  }

  // remove nan rows
  x.conservativeResize(last + 1, x.cols());
  if (weights.size() > 0)
    weights.conservativeResize(last + 1);
}

//! trims all elements in the matrix to the interval `[lower, upper]`.
//! @param x Data matrix.
//! @param lower Lower bound of the interval.
//! @param upper Upper bound of the interval.
inline void
trim(Eigen::MatrixXd& x, const double& lower, const double& upper)
{
  // code of std::for_each (save some compile time by not including <algorithm>)
  auto it = x.data();
  auto last = x.data() + x.size();
  for (; it != last; ++it) {
    if (!std::isnan(*it))
      *it = std::min(std::max(*it, lower), upper);
  }
}

//! trims all elements in the matrix to the interval `[lower, upper]`.
//! @param x Data matrix.
//! @param lower Lower bound of the interval.
//! @param upper Upper bound of the interval.
inline void
trim(Eigen::VectorXd& x, const double& lower, const double& upper)
{
  // code of std::for_each (save some compile time by not including <algorithm>)
  auto it = x.data();
  auto last = x.data() + x.size();
  for (; it != last; ++it) {
    if (!std::isnan(*it))
      *it = std::min(std::max(*it, lower), upper);
  }
}

//! check if all elements are contained in the unit cube.
//! @param u Copula data.
//! @return `true` if all data lie in the unit cube; throws an error otherwise.
inline bool
check_if_in_unit_cube(const Eigen::MatrixXd& u)
{
  bool any_outside = (u.array() < 0.0).any() || (u.array() > 1.0).any();
  if (any_outside) {
    throw std::runtime_error("all data must be contained in [0, 1]^d.");
  }
  return !any_outside;
}

//! swap the columns of a two-column matrix
//! @param u The matrix.
//! @return a new matrix v with `v.col(0) = u.col(1)`, `v.col(1) = u.col(0)`.
inline Eigen::MatrixXd
swap_cols(Eigen::MatrixXd u)
{
  u.col(0).swap(u.col(1));
  return u;
}

inline Eigen::VectorXd
unique(const Eigen::VectorXd& x)
{
  std::vector<double> v(x.data(), x.data() + x.size());
  std::sort(v.begin(), v.end()); // 1 1 2 2 3 3 3 4 4 5 5 6 7
  auto last = std::unique(v.begin(), v.end());
  v.erase(last, v.end());
  return Eigen::Map<Eigen::VectorXd>(&v[0], v.size());
}

//! computes the inverse \f$ f^{-1} \f$ of a function \f$ f \f$ by the
//! bisection method.
//!
//! @param x Evaluation points.
//! @param f The function to invert.
//! @param lb Lower bound.
//! @param ub Upper bound.
//! @param n_iter The number of iterations for the bisection (defaults to 35,
//! guaranteeing an accuracy of 0.5^35 ~= 6e-11).
//!
//! @return \f$ f^{-1}(x) \f$.
inline Eigen::VectorXd
invert_f(const Eigen::VectorXd& x,
         std::function<Eigen::VectorXd(const Eigen::VectorXd&)> f,
         const double lb,
         const double ub,
         int n_iter)
{
  Eigen::VectorXd xl = Eigen::VectorXd::Constant(x.size(), lb);
  Eigen::VectorXd xh = Eigen::VectorXd::Constant(x.size(), ub);
  Eigen::VectorXd x_tmp = x;
  Eigen::VectorXd fm(x.size());
  for (int iter = 0; iter < n_iter; ++iter) {
    x_tmp = (xh + xl) / 2.0;
    fm = f(x_tmp) - x;
    xl = (fm.array() < 0).select(x_tmp, xl);
    xh = (fm.array() < 0).select(xh, x_tmp);
  }
  if (fm.array().isNaN().any()) {
    size_t n = x.size();
    for (size_t j = 0; j < n; j++) {
      if ((std::isnan)(fm(j))) {
        x_tmp(j) = std::numeric_limits<double>::quiet_NaN();
      }
    }
  }

  return x_tmp;
}

//! computes the inverse \f$ f^{-1} \f$ of a monotone increasing function
//! \f$ f \f$ by a vectorized safeguarded Newton method (the `rtsafe`
//! algorithm of Numerical Recipes, run per element).
//!
//! Each element keeps its own bracket \f$ [x_l, x_h] \f$ with
//! \f$ f(x_l) < x < f(x_h) \f$. A Newton step is taken only when it stays
//! inside the bracket *and* shrinks the step by at least a factor of two;
//! otherwise a bisection step is used. This guards against the Newton limit
//! cycles that occur when the derivative on the flat wings badly overshoots
//! a steep interior (e.g. Tawn h-functions). Elements are frozen once the
//! step is within `tol`, so a result does not depend on the other rows in
//! the batch. Assumes (as bisection does) that \f$ f \f$ is increasing with
//! \f$ f(\text{lb}) \le x \le f(\text{ub}) \f$.
//!
//! @param x Evaluation points (the target values of \f$ f \f$).
//! @param eval Fills `f` and `f'` at the currently-active (unconverged) rows;
//!   evaluating only these rows keeps the expensive transcendental work
//!   proportional to the average number of iterations, not the worst case.
//! @param lb Lower bound.
//! @param ub Upper bound.
//! @param tol Convergence tolerance on the step size.
//! @param n_iter Maximum number of iterations (the pure-bisection fallback
//!   reaches `tol` well within the default).
//!
//! @return \f$ f^{-1}(x) \f$.
inline Eigen::VectorXd
invert_f_newton(const Eigen::VectorXd& x,
                const NewtonEval& eval,
                const double lb,
                const double ub,
                const double tol,
                int n_iter)
{
  const Eigen::Index n = x.size();
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const Eigen::ArrayXd xa = x.array();
  Eigen::ArrayXd xl = Eigen::ArrayXd::Constant(n, lb);
  Eigen::ArrayXd xh = Eigen::ArrayXd::Constant(n, ub);
  Eigen::ArrayXd v = 0.5 * (xl + xh);
  Eigen::ArrayXd dx = (xh - xl).abs();
  Eigen::ArrayXd dxold = dx;
  // fv/dfv carry the last evaluation per row; frozen rows keep stale values
  // that the (cheap) full-width arithmetic below never acts on
  Eigen::ArrayXd fv = Eigen::ArrayXd::Zero(n);
  Eigen::ArrayXd dfv = Eigen::ArrayXd::Ones(n);
  std::vector<char> done(static_cast<std::size_t>(n), 0);

  for (Eigen::Index j = 0; j < n; ++j) {
    if ((std::isnan)(xa(j))) { // NaN inputs finalized immediately
      v(j) = nan;
      done[j] = 1;
    }
  }

  // evaluate f/f' at the still-active rows and scatter them back; a NaN
  // result on a live row cannot be bracketed, so freeze it as NaN
  std::vector<Eigen::Index> active;
  active.reserve(static_cast<std::size_t>(n));
  auto evaluate_active = [&]() {
    active.clear();
    for (Eigen::Index j = 0; j < n; ++j) {
      if (!done[j]) {
        active.push_back(j);
      }
    }
    if (active.empty()) {
      return;
    }
    Eigen::VectorXd v_active(active.size());
    for (std::size_t k = 0; k < active.size(); ++k) {
      v_active(k) = v(active[k]);
    }
    Eigen::VectorXd f_out, fp_out;
    eval(active, v_active, f_out, fp_out);
    for (std::size_t k = 0; k < active.size(); ++k) {
      const Eigen::Index j = active[k];
      if (!std::isfinite(f_out(k))) {
        // the function value cannot be bracketed; finalize as NaN
        v(j) = nan;
        done[j] = 1;
      } else {
        fv(j) = f_out(k) - xa(j);
        // a non-finite derivative (overflow on the flat/steep wings) forces
        // a bisection step (dfv == 0 trips both rtsafe conditions), so the
        // solver stays as robust as plain bisection there
        dfv(j) = std::isfinite(fp_out(k)) ? fp_out(k) : 0.0;
      }
    }
  };

  evaluate_active();

  for (int iter = 0; iter < n_iter; ++iter) {
    bool any = false;
    for (Eigen::Index j = 0; j < n && !any; ++j) {
      any = !done[j];
    }
    if (!any) {
      break;
    }
    // bisect where the Newton step would leave the bracket or would not at
    // least halve the previous step (the rtsafe safeguard)
    auto out_of_bracket = ((v - xh) * dfv - fv) * ((v - xl) * dfv - fv) > 0.0;
    auto too_slow = (2.0 * fv).abs() > (dxold * dfv).abs();
    auto bisect = out_of_bracket || too_slow;

    dxold = dx;
    Eigen::ArrayXd dx_bisect = 0.5 * (xh - xl);
    Eigen::ArrayXd dx_newton = fv / dfv;
    dx = bisect.select(dx_bisect, dx_newton);
    Eigen::ArrayXd v_new = bisect.select(xl + dx_bisect, v - dx_newton);
    for (Eigen::Index j = 0; j < n; ++j) {
      if (!done[j]) {
        v(j) = v_new(j);
        if (std::fabs(dx(j)) < tol) {
          done[j] = 1; // step within tolerance
        }
      }
    }

    evaluate_active();
    // maintain the bracket for rows still live after this evaluation
    for (Eigen::Index j : active) {
      if (!done[j]) {
        if (fv(j) < 0.0) {
          xl(j) = v(j);
        } else {
          xh(j) = v(j);
        }
      }
    }
  }

  return v.matrix();
}

//! expand a vector into a matrix with two columns where each row
//! contains one combination of the vector elements
//!
//! @param grid_points The vector to expand.
inline Eigen::MatrixXd
expand_grid(const Eigen::VectorXd& grid_points)
{
  ptrdiff_t m = grid_points.size();
  Eigen::MatrixXd grid_2d(m * m, 2);
  ptrdiff_t k = 0;
  for (ptrdiff_t i = 0; i < m; ++i) {
    for (ptrdiff_t j = 0; j < m; ++j) {
      grid_2d(k, 0) = grid_points(i);
      grid_2d(k, 1) = grid_points(j);
      ++k;
    }
  }
  return grid_2d;
}

//! reads data from a file to an Eigen matrix of integers.
//!
//! The function is currently **not safe** and may cause crashes when the
//! arguments are specified incorrectly.
//!
//! @param filename The name of the file to read from.
//! @param max_buffer_size The maximal buffer size.
inline Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>
read_matxs(const char* filename, int max_buffer_size)
{
  Eigen::MatrixXd temp = read_matxd(filename, max_buffer_size);
  Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> output =
    temp.cast<size_t>();
  return output;
}

//! reads data from a file to an Eigen matrix of doubles.
//!
//! The function is currently **not safe** and may cause crashes when the
//! arguments are specified incorrectly.
//!
//! @param filename The name of the file to read from.
//! @param max_buffer_size The maximal buffer size.
inline Eigen::MatrixXd
read_matxd(const char* filename, int max_buffer_size)
{
  using namespace std;

  int cols = 0, rows = 0;
  double* buff = new double[max_buffer_size];

  // Read numbers from file into buffer.
  ifstream infile;
  infile.open(filename);
  while (!infile.eof()) {
    string line;
    getline(infile, line);

    int temp_cols = 0;
    stringstream stream(line);
    while (!stream.eof()) {
      stream >> buff[cols * rows + temp_cols++];
    }
    if (temp_cols == 0) {
      continue;
    }
    if (cols == 0) {
      cols = temp_cols;
    }
    rows++;
  }

  infile.close();

  rows--;

  // Populate matrix with numbers.
  Eigen::MatrixXd result(rows, cols);
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      result(i, j) = buff[cols * i + j];
    }
  }

  delete[] buff;
  return result;
}

//! @}
}
}
