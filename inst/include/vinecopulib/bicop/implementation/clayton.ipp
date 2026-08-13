// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/misc/tools_eigen.hpp>

namespace vinecopulib {
inline ClaytonBicop::ClaytonBicop()
{
  family_ = BicopFamily::clayton;
  parameters_ = Eigen::VectorXd(1);
  parameters_lower_bounds_ = Eigen::VectorXd(1);
  parameters_upper_bounds_ = Eigen::VectorXd(1);
  parameters_ << 1e-10;
  parameters_lower_bounds_ << 1e-10;
  parameters_upper_bounds_ << 28;
}

inline double
ClaytonBicop::generator(const double& u,
                        const Eigen::Ref<const Eigen::VectorXd>& parameters)
{
  double theta = parameters(0);
  return (std::pow(u, -theta) - 1) / theta;
}

inline double
ClaytonBicop::generator_inv(const double& u,
                            const Eigen::Ref<const Eigen::VectorXd>& parameters)
{
  double theta = parameters(0);
  return std::pow(1 + theta * u, -1 / theta);
}

inline double
ClaytonBicop::generator_derivative(
  const double& u,
  const Eigen::Ref<const Eigen::VectorXd>& parameters)
{
  return (-1) * std::pow(u, -1 - parameters(0));
}

inline Eigen::VectorXd
ClaytonBicop::pdf_raw(const Eigen::MatrixXd& u,
                      const Eigen::MatrixXd& parameters)
{
  auto f = [](const double& u1,
              const double& u2,
              const Eigen::Ref<const Eigen::VectorXd>& par) {
    double theta = par(0);
    // avoid numerical issues when copula is too close to independence
    if (theta < 1e-10) {
      return 1.0;
    }
    double temp = std::log1p(theta) - (1.0 + theta) * std::log(u1 * u2);
    temp = temp - (2.0 + 1.0 / (theta)) *
                    std::log(std::pow(u1, -theta) + std::pow(u2, -theta) - 1.0);
    return std::exp(temp);
  };
  return tools_eigen::binaryExpr_or_nan(u, parameters, f);
}

inline Eigen::VectorXd
ClaytonBicop::hinv1_raw(const Eigen::MatrixXd& u,
                        const Eigen::MatrixXd& parameters)
{
  // theta is bounded above by 28 < 75, so the closed form is always used
  const Eigen::Index n = u.rows();
  Eigen::ArrayXd theta =
    tools_eigen::parameter_as_vector(parameters, 0, n).array();
  Eigen::ArrayXd u0 = u.col(0).array();
  Eigen::ArrayXd u1 = u.col(1).array();
  Eigen::ArrayXd hinv = u0.pow(theta + 1.0);
  hinv = u1 * hinv;
  hinv = hinv.pow(-theta / (theta + 1.0));
  Eigen::ArrayXd x = u0.pow(-theta);
  hinv = hinv - x + 1.0;
  hinv = hinv.pow(-1.0 / theta);
  return hinv.matrix();
}

inline Eigen::MatrixXd
ClaytonBicop::tau_to_parameters(const double& tau)
{
  Eigen::VectorXd parameters(1);
  parameters(0) = 2 * std::fabs(tau) / (1 - std::fabs(tau));
  return parameters.cwiseMax(parameters_lower_bounds_)
    .cwiseMin(parameters_upper_bounds_);
}

inline double
ClaytonBicop::parameters_to_tau(const Eigen::MatrixXd& parameters)
{
  return parameters(0) / (2 + std::fabs(parameters(0)));
}

inline Eigen::MatrixXd
ClaytonBicop::parameters_to_taildep(const Eigen::MatrixXd& parameters)
{
  Eigen::MatrixXd taildep = Eigen::MatrixXd::Zero(2, 2);
  taildep(0, 0) = std::pow(2.0, -1.0 / parameters(0)); // lower tail dependence
  return taildep;
}

inline Eigen::VectorXd
ClaytonBicop::get_start_parameters(const double tau)
{
  Eigen::VectorXd par = tau_to_parameters(tau);
  par = par.cwiseMax(parameters_lower_bounds_);
  par = par.cwiseMin(parameters_upper_bounds_);
  return par;
}
}
