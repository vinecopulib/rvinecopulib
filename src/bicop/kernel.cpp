// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#include "bicop/kernel.hpp"
#include "misc/tools_stats.hpp"

namespace vinecopulib
{
    KernelBicop::KernelBicop()
     {
         // construct default grid (equally spaced on Gaussian scale)
         size_t m = 30;
         Eigen::VectorXd grid_points(m);
         for (size_t i = 0; i < m; ++i)
             grid_points(i) = - 3.25 + i * (6.25 / (double) m);
         interp_grid_ = InterpolationGrid(
             tools_stats::pnorm(grid_points), 
             Eigen::MatrixXd::Constant(30, 30, 1.0)  // independence
         );
     }

    Eigen::VectorXd KernelBicop::pdf(
        const Eigen::Matrix<double, Eigen::Dynamic, 2>& u
    )
    {
        return interp_grid_.interpolate(u);
    }
    Eigen::VectorXd KernelBicop::hfunc1(
        const Eigen::Matrix<double, Eigen::Dynamic, 2>& u
    )
    {
        return interp_grid_.intergrate_1d(u, 1);
    }
    Eigen::VectorXd KernelBicop::hfunc2(
        const Eigen::Matrix<double, Eigen::Dynamic, 2>& u
    )
    {
        return interp_grid_.intergrate_1d(u, 2);
    }
    Eigen::VectorXd KernelBicop::hinv1(
        const Eigen::Matrix<double, Eigen::Dynamic, 2>& u
    )
    {
        return hinv1_num(u);
    }
    Eigen::VectorXd KernelBicop::hinv2(
        const Eigen::Matrix<double, Eigen::Dynamic, 2>& u
    )
    {
        return hinv2_num(u);
    }

    // TODO
    double KernelBicop::parameters_to_tau(const Eigen::VectorXd&)
    {
        throw std::runtime_error(
                "parameters_to_tau not yet implemented for kernel estimator"
        );
    }

    double KernelBicop::calculate_npars()
    {
        return npars_;
    }

    void KernelBicop::flip()
    {
        interp_grid_.flip();
    }

    Eigen::MatrixXd KernelBicop::tau_to_parameters(const double& tau)
    {
        return vinecopulib::no_tau_to_parameters(tau);
    }
}
