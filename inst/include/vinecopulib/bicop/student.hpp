// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <vinecopulib/bicop/elliptical.hpp>

namespace vinecopulib {
//! @brief The Student t copula.
//!
//! This class is used in the implementation underlying the Bicop class.
//! Users should not use AbstractBicop or derived classes directly, but
//! always work with the Bicop interface.
//!
//! @literature
//! Joe, Harry. Dependence modeling with copulas. CRC Press, 2014.
class StudentBicop : public EllipticalBicop
{
public:
  // constructor
  StudentBicop();

private:
  // evaluation leaves (`parameters` is m x 2, m in {1, n}); these loop per row
  // because the t-distribution helpers take scalar parameters
  Eigen::VectorXd pdf_raw(const Eigen::MatrixXd& u,
                          const Eigen::MatrixXd& parameters);

  Eigen::VectorXd cdf(const Eigen::MatrixXd& u,
                      const Eigen::MatrixXd& parameters);

  Eigen::VectorXd hfunc1_raw(const Eigen::MatrixXd& u,
                             const Eigen::MatrixXd& parameters);

  Eigen::VectorXd hinv1_raw(const Eigen::MatrixXd& u,
                            const Eigen::MatrixXd& parameters);

  // single source of truth for the math (shared by the state-based and
  // parameter-aware leaves)
  static Eigen::VectorXd pdf_impl(const Eigen::MatrixXd& u,
                                  double rho,
                                  double nu);
  static Eigen::VectorXd cdf_impl(const Eigen::MatrixXd& u,
                                  double rho,
                                  double nu);
  static Eigen::VectorXd hfunc1_impl(const Eigen::MatrixXd& u,
                                     double rho,
                                     double nu);
  static Eigen::VectorXd hinv1_impl(const Eigen::MatrixXd& u,
                                    double rho,
                                    double nu);

  Eigen::MatrixXd tau_to_parameters(const double& tau);

  Eigen::MatrixXd parameters_to_taildep(const Eigen::MatrixXd& parameters);

  Eigen::VectorXd get_start_parameters(const double tau);
};
}

#include <vinecopulib/bicop/implementation/student.ipp>
