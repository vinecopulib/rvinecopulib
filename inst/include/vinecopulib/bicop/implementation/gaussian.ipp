// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/misc/tools_eigen.hpp>
#include <vinecopulib/misc/tools_stats.hpp>

namespace vinecopulib {
inline GaussianBicop::GaussianBicop()
{
  family_ = BicopFamily::gaussian;
  parameters_ = Eigen::VectorXd(1);
  parameters_lower_bounds_ = Eigen::VectorXd(1);
  parameters_upper_bounds_ = Eigen::VectorXd(1);
  parameters_ << 0;
  parameters_lower_bounds_ << -1;
  parameters_upper_bounds_ << 1;
}

inline Eigen::VectorXd
GaussianBicop::pdf_raw(const Eigen::MatrixXd& u,
                       const Eigen::MatrixXd& parameters)
{
  const Eigen::Index n = u.rows();
  if (parameters.rows() == 1) {
    // broadcast: scalar inverse-Cholesky form (bit-identical to the historical
    // state-based path)
    double rho = parameters(0, 0);
    Eigen::Matrix2d L;
    L(0, 0) = 1;
    L(1, 1) = 1 / sqrt(1.0 - pow(rho, 2.0));
    L(0, 1) = -rho * L(1, 1);
    L(1, 0) = 0;

    Eigen::VectorXd f = Eigen::VectorXd::Ones(n);
    Eigen::MatrixXd tmp = tools_stats::qnorm(u);
    f = f.cwiseQuotient(tools_stats::dnorm(tmp).rowwise().prod());
    tmp = tmp * L;
    f = f.cwiseProduct(tools_stats::dnorm(tmp).rowwise().prod());
    return f / sqrt(1.0 - pow(rho, 2.0));
  }

  Eigen::ArrayXd rho =
    tools_eigen::parameter_as_vector(parameters, 0, n).array();
  Eigen::ArrayXd s = (1.0 - rho.square()).sqrt(); // sqrt(1 - rho^2)

  Eigen::MatrixXd z = tools_stats::qnorm(u);
  Eigen::MatrixXd zL(n, 2);
  zL.col(0) = z.col(0);
  zL.col(1) = (z.col(1).array() - rho * z.col(0).array()) / s;

  Eigen::VectorXd f = tools_stats::dnorm(zL).rowwise().prod();
  f = f.cwiseQuotient(tools_stats::dnorm(z).rowwise().prod());
  f = (f.array() / s).matrix();
  return f;
}

inline Eigen::VectorXd
GaussianBicop::cdf(const Eigen::MatrixXd& u, const Eigen::MatrixXd& parameters)
{
  Eigen::MatrixXd z = tools_stats::qnorm(u);
  if (parameters.rows() == 1) {
    return tools_stats::pbvnorm(z, parameters(0, 0));
  }
  const Eigen::Index n = u.rows();
  Eigen::VectorXd p(n);
  for (Eigen::Index i = 0; i < n; ++i) {
    p(i) = tools_stats::pbvnorm(z.row(i), parameters(i, 0))(0);
  }
  return p;
}

inline Eigen::VectorXd
GaussianBicop::hfunc1_raw(const Eigen::MatrixXd& u,
                          const Eigen::MatrixXd& parameters)
{
  Eigen::MatrixXd z = tools_stats::qnorm(u);
  if (parameters.rows() == 1) {
    // broadcast: scalar form (computes sqrt(1 - rho^2) once; bit-identical to
    // the historical state-based path)
    double rho = parameters(0, 0);
    Eigen::VectorXd h = (z.col(1) - rho * z.col(0)) / sqrt(1.0 - pow(rho, 2.0));
    return tools_stats::pnorm(h);
  }
  const Eigen::Index n = u.rows();
  Eigen::ArrayXd rho =
    tools_eigen::parameter_as_vector(parameters, 0, n).array();
  Eigen::VectorXd h =
    ((z.col(1).array() - rho * z.col(0).array()) / (1.0 - rho.square()).sqrt())
      .matrix();
  return tools_stats::pnorm(h);
}

inline Eigen::VectorXd
GaussianBicop::hinv1_raw(const Eigen::MatrixXd& u,
                         const Eigen::MatrixXd& parameters)
{
  Eigen::MatrixXd z = tools_stats::qnorm(u);
  if (parameters.rows() == 1) {
    double rho = parameters(0, 0);
    Eigen::VectorXd hinv =
      z.col(1) * sqrt(1.0 - pow(rho, 2.0)) + rho * z.col(0);
    return tools_stats::pnorm(hinv);
  }
  const Eigen::Index n = u.rows();
  Eigen::ArrayXd rho =
    tools_eigen::parameter_as_vector(parameters, 0, n).array();
  Eigen::VectorXd hinv =
    (z.col(1).array() * (1.0 - rho.square()).sqrt() + rho * z.col(0).array())
      .matrix();
  return tools_stats::pnorm(hinv);
}

inline Eigen::VectorXd
GaussianBicop::get_start_parameters(const double tau)
{
  return tau_to_parameters(tau);
}

inline Eigen::MatrixXd
GaussianBicop::tau_to_parameters(const double& tau)
{
  Eigen::VectorXd parameters = this->parameters_;
  parameters(0) = std::sin(tau * constant::pi / 2);
  return parameters;
}

inline Eigen::MatrixXd
GaussianBicop::parameters_to_taildep(const Eigen::MatrixXd&)
{
  return Eigen::MatrixXd::Zero(2, 2);
}
}
