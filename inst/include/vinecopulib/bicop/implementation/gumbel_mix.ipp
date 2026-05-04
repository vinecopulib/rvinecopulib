// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <cmath>
#include <vinecopulib/misc/tools_eigen.hpp>
#include <vinecopulib/misc/tools_stats.hpp>

namespace vinecopulib {
namespace {
inline double
xtd_param_to_theta(const double par)
{
  return std::abs(par) + 1.0;
}

inline Eigen::VectorXd
gumbel_cdf_component(const Eigen::MatrixXd& u, const double theta)
{
  auto f = [theta](const double& u1, const double& u2) {
    const double t = std::pow(-std::log(u1), theta) +
                     std::pow(-std::log(u2), theta);
    return std::exp(-std::pow(t, 1.0 / theta));
  };
  return tools_eigen::binaryExpr_or_nan(u, f);
}

inline Eigen::VectorXd
gumbel_pdf_component(const Eigen::MatrixXd& u, const double theta)
{
  const double theta_inv = 1.0 / theta;
  auto f = [theta, theta_inv](const double& u1, const double& u2) {
    const double t = std::pow(-std::log(u1), theta) +
                     std::pow(-std::log(u2), theta);
    const double temp = -std::pow(t, theta_inv) +
                        (2.0 * theta_inv - 2.0) * std::log(t) +
                        (theta - 1.0) * std::log(std::log(u1) * std::log(u2)) -
                        std::log(u1 * u2) +
                        std::log1p((theta - 1.0) * std::pow(t, -theta_inv));
    return std::exp(temp);
  };
  return tools_eigen::binaryExpr_or_nan(u, f);
}

inline Eigen::VectorXd
gumbel_hfunc1_component(const Eigen::MatrixXd& u, const double theta)
{
  Eigen::VectorXd c = gumbel_cdf_component(u, theta);
  auto f = [theta](const double& u1, const double& cdf_val) {
    const double gd_u1 = std::pow(std::log(1.0 / u1), theta - 1.0) *
                         (-theta / u1);
    const double gd_c = std::pow(std::log(1.0 / cdf_val), theta - 1.0) *
                        (-theta / cdf_val);
    const double h = gd_u1 / gd_c;
    return std::isnan(h) ? 0.5 : std::min(h, 1.0);
  };
  return tools_eigen::binaryExpr_or_nan(
    (Eigen::MatrixXd(u.rows(), 2) << u.col(0), c).finished(), f);
}

inline Eigen::VectorXd
gumbel_hfunc2_component(const Eigen::MatrixXd& u, const double theta)
{
  return gumbel_hfunc1_component(tools_eigen::swap_cols(u), theta);
}

inline Eigen::VectorXd
gumbel90_cdf_component(const Eigen::MatrixXd& u, const double theta)
{
  Eigen::MatrixXd u_rot = tools_stats::rotate_data(u, 90);
  return (u.col(1).array() - gumbel_cdf_component(u_rot, theta).array())
    .matrix();
}

inline Eigen::VectorXd
gumbel90_pdf_component(const Eigen::MatrixXd& u, const double theta)
{
  return gumbel_pdf_component(tools_stats::rotate_data(u, 90), theta);
}

inline Eigen::VectorXd
gumbel90_hfunc1_component(const Eigen::MatrixXd& u, const double theta)
{
  Eigen::MatrixXd u_rot = tools_stats::rotate_data(u, 90);
  return gumbel_hfunc2_component(u_rot, theta);
}

inline Eigen::VectorXd
gumbel90_hfunc2_component(const Eigen::MatrixXd& u, const double theta)
{
  Eigen::MatrixXd u_rot = tools_stats::rotate_data(u, 90);
  return 1.0 - gumbel_hfunc1_component(u_rot, theta).array();
}
} // namespace

inline GumbelMixBicop::GumbelMixBicop()
{
  family_ = BicopFamily::gumbel_mix;
  parameters_ = Eigen::VectorXd(2);
  parameters_lower_bounds_ = Eigen::VectorXd(2);
  parameters_upper_bounds_ = Eigen::VectorXd(2);
  parameters_ << 0, 0;
  parameters_lower_bounds_ << -3, -3;
  parameters_upper_bounds_ << 3, 3;
}

inline void
GumbelMixBicop::set_parameters(const Eigen::MatrixXd& par)
{
  ParBicop::set_parameters(par);
}

inline Eigen::VectorXd
GumbelMixBicop::pdf_raw(const Eigen::MatrixXd& u)
{
  const double theta0 = xtd_param_to_theta(parameters_(0));
  const double theta90 = xtd_param_to_theta(parameters_(1));
  return 0.5 * gumbel_pdf_component(u, theta0) +
         0.5 * gumbel90_pdf_component(u, theta90);
}

inline Eigen::VectorXd
GumbelMixBicop::cdf(const Eigen::MatrixXd& u)
{
  const double theta0 = xtd_param_to_theta(parameters_(0));
  const double theta90 = xtd_param_to_theta(parameters_(1));
  return 0.5 * gumbel_cdf_component(u, theta0) +
         0.5 * gumbel90_cdf_component(u, theta90);
}

inline Eigen::VectorXd
GumbelMixBicop::hfunc1_raw(const Eigen::MatrixXd& u)
{
  const double theta0 = xtd_param_to_theta(parameters_(0));
  const double theta90 = xtd_param_to_theta(parameters_(1));
  return 0.5 * gumbel_hfunc1_component(u, theta0) +
         0.5 * gumbel90_hfunc1_component(u, theta90);
}

inline Eigen::VectorXd
GumbelMixBicop::hfunc2_raw(const Eigen::MatrixXd& u)
{
  const double theta0 = xtd_param_to_theta(parameters_(0));
  const double theta90 = xtd_param_to_theta(parameters_(1));
  return 0.5 * gumbel_hfunc2_component(u, theta0) +
         0.5 * gumbel90_hfunc2_component(u, theta90);
}

inline Eigen::VectorXd
GumbelMixBicop::hinv1_raw(const Eigen::MatrixXd& u)
{
  return hinv1_num(u);
}

inline Eigen::VectorXd
GumbelMixBicop::hinv2_raw(const Eigen::MatrixXd& u)
{
  return hinv2_num(u);
}

inline double
GumbelMixBicop::parameters_to_tau(const Eigen::MatrixXd& parameters)
{
  const double theta0 = xtd_param_to_theta(parameters(0));
  const double theta90 = xtd_param_to_theta(parameters(1));
  const double tau0 = (theta0 - 1.0) / theta0;
  const double tau90 = (theta90 - 1.0) / theta90;
  return 0.5 * tau0 - 0.5 * tau90;
}

inline Eigen::MatrixXd
GumbelMixBicop::tau_to_parameters(const double& tau)
{
  return no_tau_to_parameters(tau);
}

inline Eigen::VectorXd
GumbelMixBicop::get_start_parameters(const double)
{
  Eigen::VectorXd parameters(2);
  parameters << 0.5, 0.5;
  return parameters;
}
}
