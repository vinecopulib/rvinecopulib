// Copyright © 2018 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <cmath>

#ifndef M_PI
#define M_PI       3.14159265358979323846
#endif

namespace vinecopulib {
inline Eigen::VectorXd EllipticalBicop::hfunc2(
    const Eigen::Matrix<double, Eigen::Dynamic, 2> &u
)
{
    return hfunc1(tools_eigen::swap_cols(u));
}

inline Eigen::VectorXd EllipticalBicop::hinv2(
    const Eigen::Matrix<double, Eigen::Dynamic, 2> &u
)
{
    return hinv1(tools_eigen::swap_cols(u));
}

inline double
EllipticalBicop::parameters_to_tau(const Eigen::MatrixXd &parameters)
{
    double tau = (2 / M_PI) * asin(parameters(0));
    return tau;
}
}
