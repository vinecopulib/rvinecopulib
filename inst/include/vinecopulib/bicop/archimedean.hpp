// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <vinecopulib/bicop/parametric.hpp>

namespace vinecopulib {
//! @brief An abstract class for Archimedean copula families.
//!
//! This class is used in the implementation underlying the Bicop class.
//! Users should not use AbstractBicop or derived classes directly, but
//! always work with the Bicop interface.
//!
//! @literature
//! Joe, Harry. Dependence modeling with copulas. CRC Press, 2014.
class ArchimedeanBicop : public ParBicop
{
private:
  // cdf, hfunctions and inverses
  // Eigen::VectorXd pdf(const Eigen::MatrixXd &u);

  Eigen::VectorXd cdf(const Eigen::MatrixXd& u);

  Eigen::VectorXd hfunc1_raw(const Eigen::MatrixXd& u);

  Eigen::VectorXd hfunc2_raw(const Eigen::MatrixXd& u);

  Eigen::VectorXd hinv1_raw(const Eigen::MatrixXd& u);

  Eigen::VectorXd hinv2_raw(const Eigen::MatrixXd& u);

  // generator, its inverse and derivative
  virtual double generator(const double& u) = 0;

  virtual double generator_inv(const double& u) = 0;

  virtual double generator_derivative(const double& u) = 0;

  // virtual double generator_derivative2(const double &u) = 0;

  Eigen::VectorXd get_start_parameters(const double tau);
};
}

#include <vinecopulib/bicop/implementation/archimedean.ipp>
