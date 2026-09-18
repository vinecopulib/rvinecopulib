// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <Eigen/Dense>
#include <utility>
#include <vector>
#include <vinecopulib/misc/tools_eigen.hpp>

namespace vinecopulib {

namespace tools_interpolation {
//! A class for bilinear interpolation of bivariate copulas
//!
//! The class is used for implementing kernel estimators. It makes storing the
//! observations obsolete and allows for fast numerical integration. The
//! interpolant is piecewise bilinear, so its mass over a grid-aligned
//! rectangle is available in closed form; see `rect_mass()`.
class InterpolationGrid
{
public:
  InterpolationGrid() = default;

  InterpolationGrid(const Eigen::VectorXd& grid_points,
                    const Eigen::MatrixXd& values,
                    int norm_maxiter = 25);

  Eigen::MatrixXd get_values() const;

  void set_values(const Eigen::MatrixXd& values, int norm_maxiter = 25);

  void flip();

  Eigen::VectorXd interpolate(const tools_eigen::ConstMatRef& x);

  Eigen::VectorXd integrate_1d(const tools_eigen::ConstMatRef& u,
                               size_t cond_var);

  //! solves `integrate_1d(...) = p` for the non-conditioning coordinate;
  //! direct inversion of the piecewise-quadratic conditional cdf (used for
  //! the h-function inverses of kernel copulas).
  Eigen::VectorXd inverse_integrate_1d(const tools_eigen::ConstMatRef& u,
                                       size_t cond_var);

  Eigen::VectorXd integrate_2d(const tools_eigen::ConstMatRef& u);

  //! @brief probability of one rectangle, without the cancellation a
  //! difference of four `integrate_2d()` values carries.
  double rect_mass(double a1, double b1, double a2, double b2) const;

  //! @brief probability of one interval in the free coordinate, given the
  //! other.
  double cond_interval_mass(double u_cond,
                            double lo,
                            double hi,
                            size_t cond_var) const;

private:
  // the grid line at a fixed conditioning coordinate; `cond_knot` evaluates a
  // knot of it on demand, so no caller has to materialize the line
  //! @brief A grid line at a fixed conditioning coordinate.
  struct CondLine
  {
    ptrdiff_t cell;
    double x2x, xx1, x2x1;
    size_t cond_var;
  };
  CondLine cond_line(double u_cond, size_t cond_var) const;
  double cond_knot(const CondLine& line, ptrdiff_t j) const;

  // the weights of the two nodes of a cell `[g0, g1]`, integrating the linear
  // basis over the cell's overlap with `[a, b]`; the quadrature every integral
  // here is built from
  static std::pair<double, double> cell_weights(double g0,
                                                double g1,
                                                double a,
                                                double b);
  ptrdiff_t interval_weights(double lo, double hi, Eigen::VectorXd& w) const;
  void row_integrals(double u, Eigen::VectorXd& out) const;
  // normalizes the grid margins; internal only (callers must refresh the
  // cached integrals afterwards, as the ctor and set_values do)
  void normalize_margins(int max_iter);
  void update_weights();
  ptrdiff_t binary_search(double x);
  ptrdiff_t find_cell(double x) const;
  void update_cell_lookup();
  void update_cached_integrals();
  double cond_quantile(double u_cond,
                       double p,
                       size_t cond_var,
                       Eigen::VectorXd& knots) const;
  double int_on_grid(double upr, const Eigen::VectorXd& vals) const;

  Eigen::VectorXd grid_points_;
  Eigen::MatrixXd values_;
  // bucket acceleration table for cell searches; built once (the grid is
  // immutable after construction)
  std::vector<ptrdiff_t> cell_lookup_;
  // trapezoid weights of `grid_points_`, so that `weights_.dot(v)` integrates
  // the piecewise linear function through (grid_points_, v) over [0, 1];
  // built once alongside `cell_lookup_`
  Eigen::VectorXd weights_;
  // cumulative row integrals R(k, j) = int_0^{grid_j} values_(k, .);
  // refreshed eagerly whenever values_ changes (lazy caching would race
  // when a shared grid is evaluated from multiple threads)
  Eigen::MatrixXd row_cum_int_;
};
}
}

#include <vinecopulib/misc/implementation/tools_interpolation.ipp>
