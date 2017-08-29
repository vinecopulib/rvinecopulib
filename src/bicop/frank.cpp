// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/bicop/frank.hpp>
#include <vinecopulib/misc/tools_c.h>

namespace vinecopulib
{
    FrankBicop::FrankBicop()
    {
        family_ = BicopFamily::frank;
        parameters_ = Eigen::VectorXd(1);
        parameters_lower_bounds_ = Eigen::VectorXd(1);
        parameters_upper_bounds_ = Eigen::VectorXd(1);
        parameters_ << 0;
        parameters_lower_bounds_ << -200;
        parameters_upper_bounds_ << 200;
    }

    double FrankBicop::generator(const double& u)
    {
        double theta = double(this->parameters_(0));
        return (-1)*std::log((std::exp(-theta*u)-1)/(std::exp(-theta)-1));
    }
    double FrankBicop::generator_inv(const double& u)
    {
        double theta = double(this->parameters_(0));
        return (-1/theta)*std::log(1+std::exp(-theta-u)-std::exp(-u));
    }

    double FrankBicop::generator_derivative(const double& u)
    {
        double theta = double(this->parameters_(0));
        return theta/(1-std::exp(theta*u));
    }

    double FrankBicop::generator_derivative2(const double& u)
    {
        double theta = double(this->parameters_(0));
        return std::pow(theta,2)/std::pow(std::exp(theta*u/2) - std::exp(-theta*u/2), 2);
    }

    Eigen::MatrixXd FrankBicop::tau_to_parameters(const double& tau)
    {
        Eigen::VectorXd tau2 = Eigen::VectorXd::Constant(1, std::fabs(tau));
        auto f = [&](const Eigen::VectorXd &v) {
            return Eigen::VectorXd::Constant(1, std::fabs(parameters_to_tau(v)));
        };
        return tools_eigen::invert_f(tau2, f, -100+1e-6, 100);
    }

    double FrankBicop::parameters_to_tau(const Eigen::MatrixXd& parameters)
    {
        double par = parameters(0);
        double tau = 1 - 4/par;
        double d = debyen(std::fabs(par), 1) / std::fabs(par);
        if (par < 0) {
            d = d - par/2;
        }
        tau = tau + (4/par) * d;
        return tau;
    }

    Eigen::VectorXd FrankBicop::get_start_parameters(const double tau)
    {
        return tau_to_parameters(tau);
    }
}
