// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/misc/tools_stl.hpp>

namespace vinecopulib {
inline Eigen::VectorXd ArchimedeanBicop::pdf(
    const Eigen::Matrix<double, Eigen::Dynamic, 2> &u
)
{
    auto f = [this](const double &u, const double &v) {
        double temp = generator_inv(generator(u) + generator(v));
        temp = log(std::abs(generator_derivative2(temp))) - 
            3.0 * log(std::abs(generator_derivative(temp)));
        temp += std::log(std::abs(generator_derivative(u)));
        temp += std::log(std::abs(generator_derivative(v)));
        return std::exp(temp);
    };
    return tools_eigen::binaryExpr_or_nan(u, f);
}

inline Eigen::VectorXd ArchimedeanBicop::cdf(
    const Eigen::Matrix<double, Eigen::Dynamic, 2> &u
)
{
    auto f = [this](const double &u, const double &v) {
        return generator_inv(generator(u) + generator(v));
    };
    return tools_eigen::binaryExpr_or_nan(u, f);
}

inline Eigen::VectorXd ArchimedeanBicop::hfunc1(
    const Eigen::Matrix<double, Eigen::Dynamic, 2> &u
)
{
    auto f = [this](const double &u, const double &v) {
        double temp = generator_inv(generator(u) + generator(v));
        return generator_derivative(u) / generator_derivative(temp);
    };
    return tools_eigen::binaryExpr_or_nan(u, f);
}

inline Eigen::VectorXd ArchimedeanBicop::hfunc2(
    const Eigen::Matrix<double, Eigen::Dynamic, 2> &u
)
{
    return hfunc1(tools_eigen::swap_cols(u));
}

inline Eigen::VectorXd ArchimedeanBicop::hinv1(
    const Eigen::Matrix<double, Eigen::Dynamic, 2> &u
)
{
    Eigen::VectorXd hinv = hinv1_num(u);
    return hinv;
}

inline Eigen::VectorXd ArchimedeanBicop::hinv2(
    const Eigen::Matrix<double, Eigen::Dynamic, 2> &u
)
{
    return hinv1(tools_eigen::swap_cols(u));
}

inline Eigen::VectorXd ArchimedeanBicop::get_start_parameters(const double)
{
    Eigen::MatrixXd lb = this->get_parameters_lower_bounds();
    Eigen::VectorXd parameters = lb + Eigen::VectorXd::Constant(2, 0.1);
    return parameters;
}
}
