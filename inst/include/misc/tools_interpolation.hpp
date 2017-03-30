// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include "misc/tools_eigen.hpp"

namespace vinecopulib
{
    //! A class for cubic spline interpolation of bivariate copulas
    //!
    //! The class is used for implementing kernel estimators. It makes storing the
    //! observations obsolete and allows for fast numerical integration.
    class InterpolationGrid
    {
    public:
        InterpolationGrid() {}
        InterpolationGrid(const Eigen::VectorXd& grid_points, const Eigen::MatrixXd& values);

        void flip();

        Eigen::VectorXd interpolate(const Eigen::MatrixXd& x);
        Eigen::VectorXd intergrate_1d(const Eigen::MatrixXd& u, size_t cond_var);

    private:
        // Utility functions for spline Interpolation
        double cubic_poly(const double& x, const Eigen::VectorXd& a);
        double cubic_indef_integral(const double& x, const Eigen::VectorXd& a);
        double cubic_integral(const double& lower, const double& upper, const Eigen::VectorXd& a);
        Eigen::VectorXd find_coefs(const Eigen::VectorXd& vals, const Eigen::VectorXd& grid);
        double interp_on_grid(const double& x, const Eigen::VectorXd& vals, const Eigen::VectorXd& grid);

        // Utility functions for integration
        double int_on_grid(const double& upr, const Eigen::VectorXd& vals, const Eigen::VectorXd& grid);

        Eigen::VectorXd grid_points_;
        Eigen::MatrixXd values_;
    };
}
