// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/misc/tools_integration.hpp>

namespace vinecopulib {
inline Bb7Bicop::Bb7Bicop()
{
    family_ = BicopFamily::bb7;
    parameters_ = Eigen::VectorXd(2);
    parameters_lower_bounds_ = Eigen::VectorXd(2);
    parameters_upper_bounds_ = Eigen::VectorXd(2);
    parameters_ << 1, 1;
    parameters_lower_bounds_ << 1, 0;
    parameters_upper_bounds_ << 6, 30;
}

inline double Bb7Bicop::generator(const double &u)
{
    double theta = double(parameters_(0));
    double delta = double(parameters_(1));
    return std::pow(1 - std::pow(1 - u, theta), -delta) - 1;
}

inline double Bb7Bicop::generator_inv(const double &u)
{
    double theta = double(parameters_(0));
    double delta = double(parameters_(1));
    return 1 - std::pow(1 - std::pow(1 + u, -1 / delta), 1 / theta);
}

inline double Bb7Bicop::generator_derivative(const double &u)
{
    double theta = double(parameters_(0));
    double delta = double(parameters_(1));
    double res =
        delta * theta * std::pow(1 - std::pow(1 - u, theta), -1 - delta);
    return -res * std::pow(1 - u, theta - 1);
}

inline double Bb7Bicop::generator_derivative2(const double &u)
{
    double theta = double(parameters_(0));
    double delta = double(parameters_(1));
    double tmp = std::pow(1 - u, theta);
    double res = delta * theta * std::pow(1 - tmp, -2 - delta) *
                 std::pow(1 - u, theta - 2);
    return res * (theta - 1 + (1 + delta * theta) * tmp);
}

inline double Bb7Bicop::parameters_to_tau(const Eigen::MatrixXd &parameters)
{
    double theta = parameters(0);
    double delta = parameters(1);
    auto f = [&theta, &delta](const double &v) {
        double tmp = std::pow(1 - v, theta);
        double res = -4 * (std::pow(1 - tmp, -delta) - 1) / (theta * delta);
        return res /
               (std::pow(1 - v, theta - 1) * std::pow(1 - tmp, -delta - 1));
    };
    return 1 + tools_integration::integrate_zero_to_one(f);
}

inline Eigen::MatrixXd Bb7Bicop::tau_to_parameters(const double &tau)
{
    return vinecopulib::no_tau_to_parameters(tau);
}
}
