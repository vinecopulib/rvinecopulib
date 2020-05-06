// Copyright © 2016-2020 Thomas Nagler and Thibault Vatter
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
  // PDF
  Eigen::VectorXd pdf_raw(const Eigen::MatrixXd& u);

  // CDF
  Eigen::VectorXd cdf(const Eigen::MatrixXd& u);

  // hfunction
  Eigen::VectorXd hfunc1_raw(const Eigen::MatrixXd& u);

  // inverse hfunction
  Eigen::VectorXd hinv1_raw(const Eigen::MatrixXd& u);

  Eigen::MatrixXd tau_to_parameters(const double& tau);

  Eigen::VectorXd get_start_parameters(const double tau);
};
}

#include <vinecopulib/bicop/implementation/student.ipp>
