// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
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
  // cdf, hfunctions and inverses (`parameters` is m x p, m in {1, n}; a single
  // row is broadcast to all observations)
  Eigen::VectorXd cdf(const Eigen::MatrixXd& u,
                      const Eigen::MatrixXd& parameters) override;

  Eigen::VectorXd hfunc1_raw(const Eigen::MatrixXd& u,
                             const Eigen::MatrixXd& parameters) override;

  Eigen::VectorXd hfunc2_raw(const Eigen::MatrixXd& u,
                             const Eigen::MatrixXd& parameters) override;

  Eigen::VectorXd hinv1_raw(const Eigen::MatrixXd& u,
                            const Eigen::MatrixXd& parameters) override;

  Eigen::VectorXd hinv2_raw(const Eigen::MatrixXd& u,
                            const Eigen::MatrixXd& parameters) override;

  // Archimedean copulas are exchangeable: the second h-function derivatives
  // are the first ones at swapped arguments/selectors
  Eigen::VectorXd hfunc2_deriv_raw(const Eigen::MatrixXd& u,
                                   const Eigen::MatrixXd& parameters,
                                   const std::string& deriv) override;

  Eigen::VectorXd hfunc2_deriv2_raw(const Eigen::MatrixXd& u,
                                    const Eigen::MatrixXd& parameters,
                                    const std::string& deriv) override;

  // generator, its inverse and derivative; `parameters` is a single parameter
  // set (a p x 1 column)
  virtual double generator(
    const double& u,
    const Eigen::Ref<const Eigen::VectorXd>& parameters) = 0;

  virtual double generator_inv(
    const double& u,
    const Eigen::Ref<const Eigen::VectorXd>& parameters) = 0;

  virtual double generator_derivative(
    const double& u,
    const Eigen::Ref<const Eigen::VectorXd>& parameters) = 0;

  // virtual double generator_derivative2(const double &u) = 0;

  Eigen::VectorXd get_start_parameters(const double tau) override;
};
}

#include <vinecopulib/bicop/implementation/archimedean.ipp>
