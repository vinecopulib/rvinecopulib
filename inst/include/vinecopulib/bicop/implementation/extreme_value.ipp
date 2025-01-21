// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/misc/tools_eigen.hpp>
#include <vinecopulib/misc/tools_integration.hpp>

namespace vinecopulib {
inline Eigen::VectorXd ExtremeValueBicop::cdf(
    const Eigen::MatrixXd &u
)
{
    auto f = [this](const double &u1, const double &u2) {
        double t = std::log(u2) / std::log(u1 * u2);
        t = pickands(t);
        t = (std::log(u1) + std::log(u2)) * t;
        return std::exp(t);
    };
    return tools_eigen::binaryExpr_or_nan(u, f);
}

inline Eigen::VectorXd ExtremeValueBicop::pdf_raw(
    const Eigen::MatrixXd &u
)
{
    auto f = [this](const double &u1, const double &u2) {
        double t = std::log(u2) / std::log(u1 * u2);
        double t2 = pickands(t);
        double t3 = pickands_derivative(t);
        double t4 = pickands_derivative2(t);

        t3 = std::pow(t2, 2) + (1 - 2 * t) * t3 * t2 - 
                (1 - t) * t * (std::pow(t3, 2) + t4 / std::log(u1 * u2));
        t2 = (std::log(u1) + std::log(u2)) * t2;

        return std::exp(t2) * t3 / (u1 * u2);
    };
    return tools_eigen::binaryExpr_or_nan(u, f);
}

inline Eigen::VectorXd ExtremeValueBicop::hfunc1_raw(
    const Eigen::MatrixXd &u
)
{
    auto f = [this](const double &u1, const double &u2) {
        double t = std::log(u2) / std::log(u1 * u2);
        double t2 = pickands(t);
        double t3 = pickands_derivative(t);
        t3 = t2 - t * t3;
        t2 = (std::log(u1) + std::log(u2)) * t2;

        return std::exp(t2) * t3 / u1;
    };
    return tools_eigen::binaryExpr_or_nan(u, f);
}

inline Eigen::VectorXd ExtremeValueBicop::hfunc2_raw(
    const Eigen::MatrixXd &u
)
{
    auto f = [this](const double &u1, const double &u2) {
        double t = std::log(u2) / std::log(u1 * u2);
        double t2 = pickands(t);
        double t3 = pickands_derivative(t);
        t3 = t2 + (1 - t) * t3;
        t2 = (std::log(u1) + std::log(u2)) * t2;

        return std::exp(t2) * t3 / u2;
    };
    return tools_eigen::binaryExpr_or_nan(u, f);
}

inline Eigen::VectorXd
ExtremeValueBicop::hinv1_raw(const Eigen::MatrixXd& u)
{
  Eigen::VectorXd hinv = hinv1_num(u);
  return hinv;
}

inline Eigen::VectorXd
ExtremeValueBicop::hinv2_raw(const Eigen::MatrixXd& u)
{
  Eigen::VectorXd hinv = hinv2_num(u);
  return hinv;
}

inline double
ExtremeValueBicop::parameters_to_tau(const Eigen::MatrixXd& par)
{
  auto old_par = this->parameters_;
  this->set_parameters(par);
  auto f = [this](const double t) {
    double A = pickands(t);
    double A2 = pickands_derivative2(t);
    return t * (1 - t) * A2 / A;
  };
  double tau = tools_integration::integrate_zero_to_one(f);
  this->parameters_ = old_par;
  return tau;
}
}
