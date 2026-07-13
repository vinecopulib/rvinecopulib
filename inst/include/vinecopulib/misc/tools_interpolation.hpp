// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <Eigen/Dense>
#include <vector>
#include <vinecopulib/misc/tools_eigen.hpp>

namespace vinecopulib {

namespace tools_interpolation {
//! A class for cubic spline interpolation of bivariate copulas
//!
//! The class is used for implementing kernel estimators. It makes storing the
//! observations obsolete and allows for fast numerical integration.
class InterpolationGrid
{
public:
  InterpolationGrid() {}

  InterpolationGrid(const Eigen::VectorXd& grid_points,
                    const Eigen::MatrixXd& values,
                    int norm_times = 3);

  Eigen::MatrixXd get_values() const;

  void set_values(const Eigen::MatrixXd& values, int norm_times = 3);

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

private:
  // normalizes the grid margins; internal only (callers must refresh the
  // cached integrals afterwards, as the ctor and set_values do)
  void normalize_margins(int times);
  Eigen::Matrix<ptrdiff_t, 1, 2> get_indices(double x0, double x1);
  ptrdiff_t binary_search(double x);
  ptrdiff_t find_cell(double x) const;
  void update_cell_lookup();
  void update_cached_integrals();
  double cond_cdf(double u_cond, double u, size_t cond_var) const;
  double cond_quantile(double u_cond, double p, size_t cond_var) const;
  double bilinear_interpolation(double z11,
                                double z12,
                                double z21,
                                double z22,
                                double x1,
                                double x2,
                                double y1,
                                double y2,
                                double x,
                                double y);
  double int_on_grid(const double& upr,
                     const Eigen::VectorXd& vals,
                     const Eigen::VectorXd& grid);

  Eigen::VectorXd grid_points_;
  Eigen::MatrixXd values_;
  // bucket acceleration table for cell searches; built once (the grid is
  // immutable after construction)
  std::vector<ptrdiff_t> cell_lookup_;
  // cumulative row integrals R(k, j) = int_0^{grid_j} values_(k, .);
  // refreshed eagerly whenever values_ changes (lazy caching would race
  // when a shared grid is evaluated from multiple threads)
  Eigen::MatrixXd row_cum_int_;
};
}
}

#include <vinecopulib/misc/implementation/tools_interpolation.ipp>
