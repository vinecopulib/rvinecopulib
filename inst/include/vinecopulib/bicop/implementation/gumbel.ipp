// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <cmath>
#include <vinecopulib/misc/tools_eigen.hpp>

namespace vinecopulib {
inline GumbelBicop::GumbelBicop()
{
  family_ = BicopFamily::gumbel;
  parameters_ = Eigen::VectorXd(1);
  parameters_lower_bounds_ = Eigen::VectorXd(1);
  parameters_upper_bounds_ = Eigen::VectorXd(1);
  parameters_ << 1;
  parameters_lower_bounds_ << 1;
  parameters_upper_bounds_ << 50;
}

inline double
GumbelBicop::generator(const double& u)
{
  return std::pow(std::log(1 / u), this->parameters_(0));
}

inline double
GumbelBicop::generator_inv(const double& u)
{
  return std::exp(-std::pow(u, 1 / this->parameters_(0)));
}

inline double
GumbelBicop::generator_derivative(const double& u)
{
  double theta = double(this->parameters_(0));
  return std::pow(std::log(1 / u), theta - 1) * (-theta / u);
}

// inline double GumbelBicop::generator_derivative2(const double &u)
//{
//    double theta = double(this->parameters_(0));
//    return (theta - 1 - std::log(u)) * std::pow(std::log(1 / u), theta - 2) *
//           (theta / std::pow(u, 2));
//}

inline Eigen::VectorXd
GumbelBicop::pdf_raw(const Eigen::MatrixXd& u)
{
  double theta = static_cast<double>(parameters_(0));
  double thetha1 = 1.0 / theta;
  auto f = [theta, thetha1](const double& u1, const double& u2) {
    double t1 = std::pow(-std::log(u1), theta) + std::pow(-std::log(u2), theta);
    double temp = -std::pow(t1, thetha1) + (2 * thetha1 - 2.0) * std::log(t1) +
                  (theta - 1.0) * std::log(std::log(u1) * std::log(u2)) -
                  std::log(u1 * u2) +
                  std::log1p((theta - 1.0) * std::pow(t1, -thetha1));
    return std::exp(temp);
  };
  return tools_eigen::binaryExpr_or_nan(u, f);
}

inline Eigen::VectorXd
GumbelBicop::hinv1_raw(const Eigen::MatrixXd& u)
{
  double theta = double(this->parameters_(0));

  // Define the lambda function for qcondgum
  auto qcondgum_func = [&theta](const double& u1, const double& u2) -> double {
    return qcondgum(u2, u1, theta);
  };

  // Use binaryExpr_or_nan to compute hinv
  return tools_eigen::binaryExpr_or_nan(u, qcondgum_func);
}

inline Eigen::MatrixXd
GumbelBicop::tau_to_parameters(const double& tau)
{
  auto par = Eigen::VectorXd::Constant(1, 1.0 / (1 - std::fabs(tau)));
  return par.cwiseMax(parameters_lower_bounds_)
    .cwiseMin(parameters_upper_bounds_);
}

inline double
GumbelBicop::parameters_to_tau(const Eigen::MatrixXd& parameters)
{
  return (parameters(0) - 1) / parameters(0);
}

inline Eigen::VectorXd
GumbelBicop::get_start_parameters(const double tau)
{
  Eigen::VectorXd par = tau_to_parameters(tau);
  par = par.cwiseMax(parameters_lower_bounds_);
  par = par.cwiseMin(parameters_upper_bounds_);
  return par;
}
}

// This is copy&paste from the VineCopula package
inline double
qcondgum(const double& q, const double& u, const double& de)
{
  double a, p, z1, z2, con, de1, dif;
  double mxdif;
  int iter;

  p = 1 - q;
  z1 = -std::log(u);
  con = std::log(1. - p) - z1 + (1. - de) * std::log(z1);
  de1 = de - 1.;
  a = std::pow(2. * std::pow(z1, de), 1. / (de));
  mxdif = 1;
  iter = 0;
  dif = .1; // needed in case first step leads to NaN
  while ((mxdif > 1.e-6) && (iter < 20)) {
    double g = a + de1 * std::log(a) + con;
    double gp = 1. + de1 / a;
    if ((std::isnan)(g) || (std::isnan)(gp) || (std::isnan)(g / gp)) {
      // added for de>50
      dif /= -2.;
    } else {
      dif = g / gp;
    }
    a -= dif;
    iter++;
    int it = 0;
    while ((a <= z1) && (it < 20)) {
      dif /= 2.;
      a += dif;
      ++it;
    }
    mxdif = fabs(dif);
  }
  z2 = std::pow(std::pow(a, de) - std::pow(z1, de), 1. / (de));
  return (std::exp(-z2));
}
