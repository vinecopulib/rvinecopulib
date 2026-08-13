// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <cmath>
#include <vinecopulib/misc/tools_eigen.hpp>

namespace vinecopulib {

inline Eigen::VectorXd
EllipticalBicop::hfunc2_raw(const Eigen::MatrixXd& u,
                            const Eigen::MatrixXd& parameters)
{
  return hfunc1_raw(tools_eigen::swap_cols(u), parameters);
}

inline Eigen::VectorXd
EllipticalBicop::hinv2_raw(const Eigen::MatrixXd& u,
                           const Eigen::MatrixXd& parameters)
{
  return hinv1_raw(tools_eigen::swap_cols(u), parameters);
}

inline double
EllipticalBicop::parameters_to_tau(const Eigen::MatrixXd& parameters)
{
  double tau = (2 / constant::pi) * asin(parameters(0));
  return tau;
}
}
