// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <vinecopulib/bicop/elliptical.hpp>

namespace vinecopulib {
//! @brief The Gaussian copula.
//!
//! This class is used in the implementation underlying the Bicop class.
//! Users should not use AbstractBicop or derived classes directly, but
//! always work with the Bicop interface.
//!
//! @literature
//! Joe, Harry. Dependence modeling with copulas. CRC Press, 2014.
class GaussianBicop : public EllipticalBicop
{
public:
  // constructor
  GaussianBicop();

private:
  // evaluation leaves (`parameters` is m x 1, m in {1, n})
  Eigen::VectorXd pdf_raw(const Eigen::MatrixXd& u,
                          const Eigen::MatrixXd& parameters);

  Eigen::VectorXd cdf(const Eigen::MatrixXd& u,
                      const Eigen::MatrixXd& parameters);

  Eigen::VectorXd hfunc1_raw(const Eigen::MatrixXd& u,
                             const Eigen::MatrixXd& parameters);

  Eigen::VectorXd hinv1_raw(const Eigen::MatrixXd& u,
                            const Eigen::MatrixXd& parameters);

  // analytic derivative leaves (see tools_deriv for the selector encoding)
  Eigen::VectorXd pdf_deriv_raw(const Eigen::MatrixXd& u,
                                const Eigen::MatrixXd& parameters,
                                const std::string& deriv);

  Eigen::VectorXd pdf_deriv2_raw(const Eigen::MatrixXd& u,
                                 const Eigen::MatrixXd& parameters,
                                 const std::string& deriv);

  Eigen::VectorXd hfunc1_deriv_raw(const Eigen::MatrixXd& u,
                                   const Eigen::MatrixXd& parameters,
                                   const std::string& deriv);

  Eigen::VectorXd hfunc1_deriv2_raw(const Eigen::MatrixXd& u,
                                    const Eigen::MatrixXd& parameters,
                                    const std::string& deriv);

  Eigen::VectorXd logpdf_deriv_raw(const Eigen::MatrixXd& u,
                                   const Eigen::MatrixXd& parameters,
                                   const std::string& deriv);

  Eigen::VectorXd logpdf_deriv2_raw(const Eigen::MatrixXd& u,
                                    const Eigen::MatrixXd& parameters,
                                    const std::string& deriv);

  Eigen::MatrixXd tau_to_parameters(const double& tau);

  // the Gaussian copula has no tail dependence in any corner
  Eigen::MatrixXd parameters_to_taildep(const Eigen::MatrixXd& parameters);

  Eigen::VectorXd get_start_parameters(const double tau);
};
}

#include <vinecopulib/bicop/implementation/gaussian.ipp>
