// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <vinecopulib/bicop/extreme_value.hpp>

namespace vinecopulib {
//! @brief The Tawn copula.
//!
//! This class is used in the implementation underlying the Bicop class.
//! Users should not use AbstractBicop or derived classes directly, but
//! always work with the Bicop interface.
//!
//! @literature
//! Joe, Harry. Dependence modeling with copulas. CRC Press, 2014.
class TawnBicop : public ExtremeValueBicop
{
public:
  // constructor
  TawnBicop();

private:
  // pickands dependence functions and its derivatives
  double pickands(const double& t) override;

  double pickands_derivative(const double& t) override;

  double pickands_derivative2(const double& t) override;

  Eigen::MatrixXd tau_to_parameters(const double& tau) override;

  Eigen::VectorXd get_start_parameters(const double) override;

  void flip() override;
};
}

#include <vinecopulib/bicop/implementation/tawn.ipp>
