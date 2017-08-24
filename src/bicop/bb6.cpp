// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/bicop/bb6.hpp>
#include <vinecopulib/misc/tools_integration.hpp>

#include <boost/math/special_functions/expm1.hpp>
#include <boost/math/special_functions/log1p.hpp>

namespace vinecopulib
{
    Bb6Bicop::Bb6Bicop()
    {
        family_ = BicopFamily::bb6;
        parameters_ = Eigen::VectorXd(2);
        parameters_lower_bounds_ = Eigen::VectorXd(2);
        parameters_upper_bounds_ = Eigen::VectorXd(2);
        parameters_ << 1, 1;
        parameters_lower_bounds_ << 1, 1;
        parameters_upper_bounds_ << 200, 200;
    }

    double Bb6Bicop::generator(const double& u)
    {
        double res = boost::math::log1p(-std::pow(1-u,parameters_(0)));
        return std::pow((-1)*res, parameters_(1));
    }

    double Bb6Bicop::generator_inv(const double& u)
    {
        double res = boost::math::expm1(-std::pow(u, 1/parameters_(1)));
        return 1-std::pow(-res, 1/parameters_(0));
    }

    double Bb6Bicop::generator_derivative(const double& u)
    {
        double theta = double(parameters_(0));
        double delta = double(parameters_(1));
        double res = boost::math::log1p(-std::pow(1-u,theta));
        res = delta * theta *std::pow((-1)*res,delta-1);
        return res*std::pow(1-u,theta-1)/(std::pow(1-u,theta)-1);
    }

    double Bb6Bicop::generator_derivative2(const double& u)
    {
        double theta = double(parameters_(0));
        double delta = double(parameters_(1));
        double tmp = std::pow(1-u,theta);
        double tmp2 = boost::math::log1p(-tmp);
        double res = std::pow((-1)*tmp2,delta-2);
        res *= ((delta-1)*theta*tmp-(tmp+theta-1)*tmp2);
        return res*delta*theta*std::pow(1-u, theta-2)/std::pow(tmp - 1,2);
    }

    double Bb6Bicop::parameters_to_tau(const Eigen::MatrixXd& parameters)
    {
        double theta = parameters(0);
        double delta = parameters(1);
        auto f = [&theta, &delta](const double& v) {
            double res = -4*(1-v-std::pow(1-v,-theta)+std::pow(1-v,-theta)*v);
            return 1/(delta*theta)*boost::math::log1p(-std::pow(1-v,theta))*res;
        };
        return 1 + tools_integration::integrate_zero_to_one(f);
    }

    Eigen::MatrixXd Bb6Bicop::tau_to_parameters(const double& tau)
    {
        return vinecopulib::no_tau_to_parameters(tau);
    }
}
