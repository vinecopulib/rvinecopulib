// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <vinecopulib/misc/tools_stl.hpp>
#include <vector>
#include <Eigen/Dense>

namespace vinecopulib {

//! Tools for working with Eigen types
namespace tools_eigen {
//! An `Eigen::Matrix` containing `bool`s (similar to `Eigen::MatrixXd`).
typedef Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> MatrixXb;

template<typename T>
Eigen::MatrixXd unaryExpr_or_nan(const Eigen::MatrixXd &x, T func)
{
    return x.unaryExpr([&func](double y) {
        return tools_stl::unaryFunc_or_nan(func, y);
    });
}

template<typename T>
Eigen::VectorXd binaryExpr_or_nan(
    const Eigen::Matrix<double, Eigen::Dynamic, 2> &u,
    T func)
{
    return u.col(0).binaryExpr(u.col(1),
                               [&func](double u1, double u2) {
                                   return tools_stl::binaryFunc_or_nan(func,
                                                                       u1,
                                                                       u2);
                               });
}

Eigen::MatrixXd nan_omit(const Eigen::MatrixXd &x);

Eigen::Matrix<double, Eigen::Dynamic, 2> swap_cols(
    Eigen::Matrix<double, Eigen::Dynamic, 2> u);

Eigen::VectorXd invert_f(
    const Eigen::VectorXd &x,
    std::function< Eigen::VectorXd(const Eigen::VectorXd &)

> f,
const double lb = 1e-20,
const double ub = 1 - 1e-20,
int n_iter = 35
);

Eigen::Matrix<double, Eigen::Dynamic, 2> expand_grid(
    const Eigen::VectorXd &grid_points
);

Eigen::MatrixXd read_matxd(const char *filename,
                           int max_buffer_size = static_cast<int>(1e6));

Eigen::Matrix <size_t, Eigen::Dynamic, Eigen::Dynamic> read_matxs(
    const char *filename, int max_buffer_size = static_cast<int>(1e6));
}

}

#include <vinecopulib/misc/implementation/tools_eigen.ipp>
