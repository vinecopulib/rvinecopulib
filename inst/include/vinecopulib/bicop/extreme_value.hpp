// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <vinecopulib/bicop/parametric.hpp>

namespace vinecopulib {
//! @brief An abstract class for extreme value copula families.
//!
//! This class is used in the implementation underlying the Bicop class.
//! Users should not use AbstractBicop or derived classes directly, but
//! always work with the Bicop interface.
//!
//! @literature
//! Joe, Harry. Dependence modeling with copulas. CRC Press, 2014.
class ExtremeValueBicop : public ParBicop
{
private:
  // pdf, cdf, hfunctions and inverses
  Eigen::VectorXd cdf(const Eigen::MatrixXd& u);

  Eigen::VectorXd pdf_raw(const Eigen::MatrixXd &u);

  Eigen::VectorXd hfunc1_raw(const Eigen::MatrixXd& u);

  Eigen::VectorXd hfunc2_raw(const Eigen::MatrixXd& u);

  Eigen::VectorXd hinv1_raw(const Eigen::MatrixXd& u);

  Eigen::VectorXd hinv2_raw(const Eigen::MatrixXd& u);

  // pickands dependence functions and its derivatives
  virtual double pickands(const double& t) = 0;

  virtual double pickands_derivative(const double& t) = 0;

  virtual double pickands_derivative2(const double& t) = 0;

  // link between Kendall's tau and the par_bicop parameter
  double parameters_to_tau(const Eigen::MatrixXd& par);
};
}

#include <vinecopulib/bicop/implementation/extreme_value.ipp>
