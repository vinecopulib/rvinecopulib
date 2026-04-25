// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <vinecopulib/bicop/elliptical.hpp>

namespace vinecopulib {
//! @brief A fixed 30/70 mixture of two Gaussian copulas.
//!
//! This class is used in the implementation underlying the Bicop class.
//! Users should not use AbstractBicop or derived classes directly, but
//! always work with the Bicop interface.
class GaussianMixBicop : public EllipticalBicop
{
public:
  // constructor
  GaussianMixBicop();

private:
  // PDF
  Eigen::VectorXd pdf_raw(const Eigen::MatrixXd& u);

  // CDF
  Eigen::VectorXd cdf(const Eigen::MatrixXd& u);

  // hfunction
  Eigen::VectorXd hfunc1_raw(const Eigen::MatrixXd& u);

  // inverse hfunction
  Eigen::VectorXd hinv1_raw(const Eigen::MatrixXd& u);

  // link between Kendall's tau and copula parameters
  double parameters_to_tau(const Eigen::MatrixXd& parameters);

  Eigen::MatrixXd tau_to_parameters(const double& tau);

  Eigen::VectorXd get_start_parameters(const double tau);
};
}

#include <vinecopulib/bicop/implementation/gaussian_mix.ipp>
