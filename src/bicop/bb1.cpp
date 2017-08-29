// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/bicop/bb1.hpp>
#include <vinecopulib/misc/tools_integration.hpp>

namespace vinecopulib
{
    Bb1Bicop::Bb1Bicop()
    {
        family_ = BicopFamily::bb1;
        parameters_ = Eigen::VectorXd(2);
        parameters_lower_bounds_ = Eigen::VectorXd(2);
        parameters_upper_bounds_ = Eigen::VectorXd(2);
        parameters_ << 0, 1;
        parameters_lower_bounds_ << 0, 1;
        parameters_upper_bounds_ << 200, 200;
    }


    double Bb1Bicop::generator(const double& u)
    {
        return std::pow(std::pow(u, -parameters_(0)) - 1, parameters_(1));
    }

    double Bb1Bicop::generator_inv(const double& u)
    {
        return std::pow(std::pow(u, 1/parameters_(1)) + 1, -1/parameters_(0));
    }

    double Bb1Bicop::generator_derivative(const double& u)
    {
        double theta = double(parameters_(0));
        double delta = double(parameters_(1));
        double res = -delta * theta * std::pow(u, -(1 + theta));
        return res*std::pow(std::pow(u, -theta) - 1, delta - 1);
    }

    double Bb1Bicop::generator_derivative2(const double& u)
    {
        double theta = double(parameters_(0));
        double delta = double(parameters_(1));
        double res = delta * theta * std::pow(std::pow(u, -theta) - 1, delta);
        res /= (std::pow(std::pow(u, theta) - 1, 2)*std::pow(u, 2));
        return res * (1 + delta * theta - (1 + theta) * std::pow(u, theta));
    }

    double Bb1Bicop::parameters_to_tau(const Eigen::MatrixXd& parameters)
    {
        return 1 - 2 / (parameters(1) * (parameters(0) + 2));
    }

    Eigen::MatrixXd Bb1Bicop::tau_to_parameters(const double& tau)
    {
        return vinecopulib::no_tau_to_parameters(tau);
    }
}
