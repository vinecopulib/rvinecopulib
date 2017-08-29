// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/bicop/bb8.hpp>
#include <vinecopulib/misc/tools_integration.hpp>

namespace vinecopulib
{
    Bb8Bicop::Bb8Bicop()
    {
        family_ = BicopFamily::bb8;
        parameters_ = Eigen::VectorXd(2);
        parameters_lower_bounds_ = Eigen::VectorXd(2);
        parameters_upper_bounds_ = Eigen::VectorXd(2);
        parameters_ << 1, 1;
        parameters_lower_bounds_ << 1, 0;
        parameters_upper_bounds_ << 200, 1;
    }

    double Bb8Bicop::generator(const double& u)
    {
        double theta = double(parameters_(0));
        double delta = double(parameters_(1));
        double res = (1-std::pow(1-delta*u,theta));
        return -std::log(res/(1-std::pow(1-delta,theta)));
    }

    double Bb8Bicop::generator_inv(const double& u)
    {
        double theta = double(parameters_(0));
        double delta = double(parameters_(1));
        double res = std::exp(-u)*(std::pow(1-delta,theta)-1);
        return (1-std::pow(1+res,1/theta))/delta;
    }

    double Bb8Bicop::generator_derivative(const double& u)
    {
        double theta = double(parameters_(0));
        double delta = double(parameters_(1));
        double res = delta*theta*std::pow(1-delta*u,theta-1);
        return -res/(1-std::pow(1-delta*u,theta));
    }

    double Bb8Bicop::generator_derivative2(const double& u)
    {
        double theta = double(parameters_(0));
        double delta = double(parameters_(1));
        double tmp = std::pow(1-delta*u,theta);
        double res = std::pow(delta,2)*theta*std::pow(1-delta*u,theta-2);
        return res*(theta-1+tmp)/std::pow(tmp-1,2);
    }

    double Bb8Bicop::parameters_to_tau(const Eigen::MatrixXd& parameters)
    {
        double theta = parameters(0);
        double delta = parameters(1);
        auto f = [theta, delta](const double t) {
            double tmp = std::pow(1-t*delta,theta);
            double res = std::log((tmp-1)/(std::pow(1-delta,theta)-1));
            return res*(1-t*delta-std::pow(1-t*delta,1-theta));
        };
        return 1 - 4/(delta*theta)*tools_integration::integrate_zero_to_one(f);
    }

    Eigen::MatrixXd Bb8Bicop::tau_to_parameters(const double& tau)
    {
        return vinecopulib::no_tau_to_parameters(tau);
    }
}
