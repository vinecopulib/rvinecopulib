// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <Eigen/Dense>
#include <functional>
#include <vector>

namespace vinecopulib {

//! Tools for working with Eigen types
namespace tools_eigen {
//! An `Eigen::Matrix` containing `bool`s (similar to `Eigen::MatrixXd`).
using MatrixXb = Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>;

//! A reference to a constant `Eigen::MatrixXd`; binds blocks (e.g.
//! `u.leftCols(2)`) and full matrices without copying.
using ConstMatRef = Eigen::Ref<const Eigen::MatrixXd>;

template<typename T>
Eigen::MatrixXd
unaryExpr_or_nan(const tools_eigen::ConstMatRef& x, const T& func)
{
  // raw-pointer loop: coefficient access through the Ref's dynamic strides
  // is measurably slower in the hot per-element paths
  const Eigen::Index n = x.rows();
  Eigen::MatrixXd out(n, x.cols());
  for (Eigen::Index j = 0; j < x.cols(); ++j) {
    const double* p = x.data() + j * x.outerStride();
    double* q = out.data() + j * n;
    for (Eigen::Index i = 0; i < n; ++i) {
      if ((std::isnan)(p[i])) {
        q[i] = std::numeric_limits<double>::quiet_NaN();
      } else {
        q[i] = func(p[i]);
      }
    }
  }
  return out;
}

template<typename T>
Eigen::VectorXd
binaryExpr_or_nan(const tools_eigen::ConstMatRef& u, const T& func)
{
  // raw-pointer loop: coefficient access through the Ref's dynamic strides
  // is measurably slower in the hot per-element paths
  const Eigen::Index n = u.rows();
  const double* p0 = u.data();
  const double* p1 = u.data() + u.outerStride();
  Eigen::VectorXd out(n);
  for (Eigen::Index i = 0; i < n; ++i) {
    if ((std::isnan)(p0[i]) || (std::isnan)(p1[i])) {
      out(i) = std::numeric_limits<double>::quiet_NaN();
    } else {
      out(i) = func(p0[i], p1[i]);
    }
  }
  return out;
}

//! @brief Applies a bivariate function row-wise with per-observation
//! parameters, propagating NaNs.
//!
//! @param u An \f$ n \times 2 \f$ matrix of evaluation points.
//! @param parameters An \f$ m \times p \f$ matrix of parameters with
//!   \f$ m \in \{1, n\} \f$. Row \f$ i \f$ holds the \f$ p \f$ parameters
//!   for observation \f$ i \f$; a single row is broadcast to all observations.
//! @param func A callable
//!   `(double u1, double u2, const Eigen::Ref<const Eigen::VectorXd>& par) ->
//!   double`.
template<typename T>
Eigen::VectorXd
binaryExpr_or_nan(const tools_eigen::ConstMatRef& u,
                  const tools_eigen::ConstMatRef& parameters,
                  const T& func)
{
  const Eigen::Index n = u.rows();
  const bool broadcast = (parameters.rows() == 1);
  // hoist the broadcast parameter set out of the loop so the common
  // single-parameter (state-based) path does not rebuild it on every row
  const Eigen::VectorXd par0 =
    broadcast ? parameters.row(0).transpose() : Eigen::VectorXd();
  const double* p0 = u.data();
  const double* p1 = u.data() + u.outerStride();
  Eigen::VectorXd out(n);
  for (Eigen::Index i = 0; i < n; ++i) {
    const double u1 = p0[i];
    const double u2 = p1[i];
    if ((std::isnan)(u1) || (std::isnan)(u2)) {
      out(i) = std::numeric_limits<double>::quiet_NaN();
    } else if (broadcast) {
      out(i) = func(u1, u2, par0);
    } else {
      const Eigen::VectorXd par_i = parameters.row(i).transpose();
      out(i) = func(u1, u2, par_i);
    }
  }
  return out;
}

//! @brief Returns the `k`-th parameter as a length-`n` vector, broadcasting
//! a single parameter set.
//!
//! @param parameters An \f$ m \times p \f$ matrix of parameters with
//!   \f$ m \in \{1, n\} \f$ (see `binaryExpr_or_nan`).
//! @param k The index of the parameter to extract.
//! @param n The desired output length.
inline Eigen::VectorXd
parameter_as_vector(const tools_eigen::ConstMatRef& parameters,
                    const Eigen::Index k,
                    const Eigen::Index n)
{
  if (parameters.rows() == 1) {
    return Eigen::VectorXd::Constant(n, parameters(0, k));
  }
  return parameters.col(k);
}

void
remove_nans(Eigen::MatrixXd& x);

void
remove_nans(Eigen::MatrixXd& x, Eigen::VectorXd& weights);

void
trim(Eigen::MatrixXd& x,
     const double& lower = 1e-10,
     const double& upper = 1 - 1e-10);

void
trim(Eigen::VectorXd& x,
     const double& lower = 1e-10,
     const double& upper = 1 - 1e-10);

bool
check_if_in_unit_cube(const Eigen::MatrixXd& u);

Eigen::MatrixXd
swap_cols(Eigen::MatrixXd u);

Eigen::VectorXd
unique(const Eigen::VectorXd& x);

Eigen::VectorXd
invert_f(const Eigen::VectorXd& x,
         std::function<Eigen::VectorXd(const Eigen::VectorXd&)> f,
         const double lb = 1e-20,
         const double ub = 1 - 1e-20,
         int n_iter = 35);

//! evaluates the function and its derivative at the active rows only.
//! Given the still-unconverged row indices and their current values, it fills
//! the (same-length) outputs with \f$ f \f$ and \f$ f' \f$ at those rows.
using NewtonEval =
  std::function<void(const std::vector<Eigen::Index>& active_rows,
                     const Eigen::VectorXd& v_active,
                     Eigen::VectorXd& f_out,
                     Eigen::VectorXd& fprime_out)>;

Eigen::VectorXd
invert_f_newton(const Eigen::VectorXd& x,
                const NewtonEval& eval,
                const double lb = 1e-20,
                const double ub = 1 - 1e-20,
                const double tol = 1e-10,
                int n_iter = 50);

Eigen::MatrixXd
expand_grid(const Eigen::VectorXd& grid_points);

Eigen::MatrixXd
read_matxd(const char* filename, int max_buffer_size = static_cast<int>(1e6));

Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>
read_matxs(const char* filename, int max_buffer_size = static_cast<int>(1e6));
}
}

#include <vinecopulib/misc/implementation/tools_eigen.ipp>
