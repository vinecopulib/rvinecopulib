// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <vinecopulib/bicop/parametric.hpp>
#include <vinecopulib/bicop/gumbel.hpp>

namespace vinecopulib {
//! @brief The Gaussian copula.
//!
//! This class is used in the implementation underlying the Bicop class.
//! Users should not use AbstractBicop or derived classes directly, but
//! always work with the Bicop interface.
//!
//! @literature
//! Joe, Harry. Dependence modeling with copulas. CRC Press, 2014.
class XtdGumbelBicop : public ParBicop
{
public:
  // constructor
  XtdGumbelBicop();

  void set_parameters(const Eigen::MatrixXd& parameters) override;

private:
  // PDF
  Eigen::VectorXd pdf_raw(const Eigen::MatrixXd& u);

  // CDF
  Eigen::VectorXd cdf(const Eigen::MatrixXd& u);

  // hfunction
  Eigen::VectorXd hfunc1_raw(const Eigen::MatrixXd& u);
  Eigen::VectorXd hfunc2_raw(const Eigen::MatrixXd& u);
  Eigen::VectorXd hinv1_raw(const Eigen::MatrixXd& u);
  Eigen::VectorXd hinv2_raw(const Eigen::MatrixXd& u);

  Eigen::MatrixXd tau_to_parameters(const double& tau);
  double parameters_to_tau(const Eigen::MatrixXd& par);

  Eigen::VectorXd get_start_parameters(const double tau);

  std::shared_ptr<GumbelBicop> gumbel_bicop_;
};
}

#include <vinecopulib/bicop/implementation/xtdgumbel.ipp>
