// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/bicop/archimedean.hpp>
#include <vinecopulib/misc/tools_stl.hpp>

namespace vinecopulib
{
    Eigen::VectorXd ArchimedeanBicop::pdf(
        const Eigen::Matrix<double, Eigen::Dynamic, 2>& u
    )
    {
        auto f = [this](const double& u, const double& v) {
            double temp = generator_inv(generator(u) + generator(v));
            temp = generator_derivative2(temp)/std::pow(generator_derivative(temp), 3.0);
            return (-1)*generator_derivative(u)*generator_derivative(v)*temp;
        };
        return u.col(0).binaryExpr(u.col(1), f);
    }

    Eigen::VectorXd ArchimedeanBicop::cdf(
            const Eigen::Matrix<double, Eigen::Dynamic, 2>& u
    )
    {
        auto f = [this](const double& u, const double& v) {
            return generator_inv(generator(u) + generator(v));
        };
        return u.col(0).binaryExpr(u.col(1), f);
    }

    Eigen::VectorXd ArchimedeanBicop::hfunc1(
        const Eigen::Matrix<double, Eigen::Dynamic, 2>& u
    )
    {
        auto f = [this](const double& u, const double& v) {
            double temp = generator_inv(generator(u) + generator(v));
            return generator_derivative(u)/generator_derivative(temp);
        };
        return u.col(0).binaryExpr(u.col(1), f);
    }

    Eigen::VectorXd ArchimedeanBicop::hfunc2(
        const Eigen::Matrix<double, Eigen::Dynamic, 2>& u
    )
    {
        return hfunc1(tools_eigen::swap_cols(u));
    }

    Eigen::VectorXd ArchimedeanBicop::hinv1(
        const Eigen::Matrix<double, Eigen::Dynamic, 2>& u
    )
    {
        Eigen::VectorXd hinv = hinv1_num(u);
        return hinv;
    }

    Eigen::VectorXd ArchimedeanBicop::hinv2(
        const Eigen::Matrix<double, Eigen::Dynamic, 2>& u
    )
    {
        return hinv1(tools_eigen::swap_cols(u));
    }
    
    Eigen::VectorXd ArchimedeanBicop::get_start_parameters(const double)
    {
        Eigen::MatrixXd lb = this->get_parameters_lower_bounds();
        Eigen::VectorXd parameters = lb + Eigen::VectorXd::Constant(2, 0.1);
        return parameters;
    }
}
