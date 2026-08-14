// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/misc/tools_eigen.hpp>

namespace vinecopulib {
// inline Eigen::VectorXd ArchimedeanBicop::pdf(
//    const Eigen::MatrixXd &u
//)
//{
//    auto f = [this](const double &u1, const double &u2) {
//        double temp = generator_inv(generator(u1) + generator(u2));
//        temp = log(std::abs(generator_derivative2(temp))) -
//            3.0 * log(std::abs(generator_derivative(temp)));
//        temp += std::log(std::abs(generator_derivative(u1)));
//        temp += std::log(std::abs(generator_derivative(u2)));
//        return std::exp(temp);
//    };
//    return tools_eigen::binaryExpr_or_nan(u, f);
//}

inline Eigen::VectorXd
ArchimedeanBicop::cdf(const Eigen::MatrixXd& u,
                      const Eigen::MatrixXd& parameters)
{
  auto f = [this](const double& u1,
                  const double& u2,
                  const Eigen::Ref<const Eigen::VectorXd>& par) {
    return generator_inv(generator(u1, par) + generator(u2, par), par);
  };
  return tools_eigen::binaryExpr_or_nan(u, parameters, f);
}

inline Eigen::VectorXd
ArchimedeanBicop::hfunc1_raw(const Eigen::MatrixXd& u,
                             const Eigen::MatrixXd& parameters)
{
  auto f = [this](const double& u1,
                  const double& u2,
                  const Eigen::Ref<const Eigen::VectorXd>& par) {
    double temp = generator_inv(generator(u1, par) + generator(u2, par), par);
    temp = generator_derivative(u1, par) / generator_derivative(temp, par);
    return std::isnan(temp) ? u2 : std::min(temp, 1.0);
  };
  return tools_eigen::binaryExpr_or_nan(u, parameters, f);
}

inline Eigen::VectorXd
ArchimedeanBicop::hfunc2_raw(const Eigen::MatrixXd& u,
                             const Eigen::MatrixXd& parameters)
{
  return hfunc1_raw(tools_eigen::swap_cols(u), parameters);
}

inline Eigen::VectorXd
ArchimedeanBicop::hfunc2_deriv_raw(const Eigen::MatrixXd& u,
                                   const Eigen::MatrixXd& parameters,
                                   const std::string& deriv)
{
  return hfunc1_deriv_raw(
    tools_eigen::swap_cols(u), parameters, tools_deriv::swap_args(deriv));
}

inline Eigen::VectorXd
ArchimedeanBicop::hfunc2_deriv2_raw(const Eigen::MatrixXd& u,
                                    const Eigen::MatrixXd& parameters,
                                    const std::string& deriv)
{
  return hfunc1_deriv2_raw(
    tools_eigen::swap_cols(u), parameters, tools_deriv::swap_args(deriv));
}

inline Eigen::VectorXd
ArchimedeanBicop::hinv1_raw(const Eigen::MatrixXd& u,
                            const Eigen::MatrixXd& parameters)
{
  return hinv1_num_raw(u, parameters);
}

inline Eigen::VectorXd
ArchimedeanBicop::hinv2_raw(const Eigen::MatrixXd& u,
                            const Eigen::MatrixXd& parameters)
{
  return hinv1_raw(tools_eigen::swap_cols(u), parameters);
}

inline Eigen::VectorXd
ArchimedeanBicop::get_start_parameters(const double)
{
  Eigen::MatrixXd lb = this->get_parameters_lower_bounds();
  Eigen::VectorXd parameters = lb + Eigen::VectorXd::Constant(2, 0.1);
  return parameters;
}
}
