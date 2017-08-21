// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <cstddef>
#include <vinecopulib/misc/tools_interpolation.hpp>
#include <vinecopulib/bicop/abstract.hpp>

namespace vinecopulib
{
    //! @brief An abstract class for kernel copulas
    //! 
    //! Evaluation functions of kernel estimators are implemented efficiently 
    //! using spline interpolation, see Nagler (2016). 
    //!
    //! This class is used in the implementation underlying the Bicop class. 
    //! Users should not use AbstractBicop or derived classes directly, but 
    //! always work with the Bicop interface.
    //! 
    //! @literature
    //! Nagler, Thomas. *kdecopula: An R Package for the Kernel Estimation of 
    //! Copula Densities*. arXiv:1603.04229 [stat.CO], 2016
    class KernelBicop : public AbstractBicop
    {
    public:
        KernelBicop();

    protected:
        Eigen::VectorXd pdf(
            const Eigen::Matrix<double, Eigen::Dynamic, 2>& u
        );
        Eigen::VectorXd cdf(
                const Eigen::Matrix<double, Eigen::Dynamic, 2>& u
        );
        Eigen::VectorXd hfunc1(
            const Eigen::Matrix<double, Eigen::Dynamic, 2>& u
        );
        Eigen::VectorXd hfunc2(
            const Eigen::Matrix<double, Eigen::Dynamic, 2>& u
        );
        Eigen::VectorXd hinv1(
            const Eigen::Matrix<double, Eigen::Dynamic, 2>& u
        );
        Eigen::VectorXd hinv2(
            const Eigen::Matrix<double, Eigen::Dynamic, 2>& u
        );

        double parameters_to_tau(const Eigen::MatrixXd & parameters);
        Eigen::MatrixXd tau_to_parameters(const double& tau);
        double calculate_npars();

        Eigen::MatrixXd get_parameters() const;
        void set_parameters(const Eigen::MatrixXd& parameters);
        void flip();

        tools_interpolation::InterpolationGrid interp_grid_;
        double npars_;
    };
}
