// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <cmath>
#include <vinecopulib/misc/tools_stats.hpp>

namespace vinecopulib {
namespace {
inline Eigen::VectorXd
gaussian_pdf_component(const Eigen::MatrixXd& u, const double rho)
{
  Eigen::Matrix2d L;
  L(0, 0) = 1;
  L(1, 1) = 1 / std::sqrt(1.0 - std::pow(rho, 2.0));
  L(0, 1) = -rho * L(1, 1);
  L(1, 0) = 0;

  Eigen::VectorXd f = Eigen::VectorXd::Ones(u.rows());
  Eigen::MatrixXd tmp = tools_stats::qnorm(u);
  f = f.cwiseQuotient(tools_stats::dnorm(tmp).rowwise().prod());
  tmp = tmp * L;
  f = f.cwiseProduct(tools_stats::dnorm(tmp).rowwise().prod());
  return f / std::sqrt(1.0 - std::pow(rho, 2.0));
}

inline Eigen::VectorXd
gaussian_hfunc1_component(const Eigen::MatrixXd& u, const double rho)
{
  Eigen::MatrixXd tmp = tools_stats::qnorm(u);
  Eigen::VectorXd h = (tmp.col(1) - rho * tmp.col(0)) /
                      std::sqrt(1.0 - std::pow(rho, 2.0));
  return tools_stats::pnorm(h);
}
} // namespace

inline GaussianMixBicop::GaussianMixBicop()
{
  family_ = BicopFamily::gaussian_mix;
  parameters_ = Eigen::VectorXd(2);
  parameters_lower_bounds_ = Eigen::VectorXd(2);
  parameters_upper_bounds_ = Eigen::VectorXd(2);
  parameters_ << 0, 0;
  parameters_lower_bounds_ << -1, -1;
  parameters_upper_bounds_ << 1, 1;
}

inline Eigen::VectorXd
GaussianMixBicop::pdf_raw(const Eigen::MatrixXd& u)
{
  const double rho1 = parameters_(0);
  const double rho2 = parameters_(1);
  return 0.3 * gaussian_pdf_component(u, rho1) +
         0.7 * gaussian_pdf_component(u, rho2);
}

inline Eigen::VectorXd
GaussianMixBicop::cdf(const Eigen::MatrixXd& u)
{
  const double rho1 = parameters_(0);
  const double rho2 = parameters_(1);
  Eigen::MatrixXd z = tools_stats::qnorm(u);
  return 0.3 * tools_stats::pbvnorm(z, rho1) +
         0.7 * tools_stats::pbvnorm(z, rho2);
}

inline Eigen::VectorXd
GaussianMixBicop::hfunc1_raw(const Eigen::MatrixXd& u)
{
  const double rho1 = parameters_(0);
  const double rho2 = parameters_(1);
  return 0.3 * gaussian_hfunc1_component(u, rho1) +
         0.7 * gaussian_hfunc1_component(u, rho2);
}

inline Eigen::VectorXd
GaussianMixBicop::hinv1_raw(const Eigen::MatrixXd& u)
{
  return hinv1_num(u);
}

inline double
GaussianMixBicop::parameters_to_tau(const Eigen::MatrixXd& parameters)
{
  const double tau1 = (2 / constant::pi) * std::asin(parameters(0));
  const double tau2 = (2 / constant::pi) * std::asin(parameters(1));
  return 0.3 * tau1 + 0.7 * tau2;
}

inline Eigen::MatrixXd
GaussianMixBicop::tau_to_parameters(const double& tau)
{
  return no_tau_to_parameters(tau);
}

inline Eigen::VectorXd
GaussianMixBicop::get_start_parameters(const double tau)
{
  Eigen::VectorXd parameters = this->parameters_;
  const double rho = std::sin(tau * constant::pi / 2);
  parameters << 0, 0;
  return parameters;
}
}
