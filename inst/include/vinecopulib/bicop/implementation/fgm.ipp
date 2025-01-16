// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/bicop/fgm.hpp>
#include <vinecopulib/misc/tools_eigen.hpp>

namespace vinecopulib {

// constructor
inline FGMBicop::FGMBicop()
{
  family_ = BicopFamily::fgm;
  parameters_ = Eigen::VectorXd::Constant(1, 0.0);
  parameters_lower_bounds_ = Eigen::VectorXd::Constant(1, -1.0);
  parameters_upper_bounds_ = Eigen::VectorXd::Constant(1, 1.0);
}

// PDF
inline Eigen::VectorXd
FGMBicop::pdf_raw(const Eigen::MatrixXd& u)
{
  const Eigen::VectorXd& u1 = u.col(0);
  const Eigen::VectorXd& u2 = u.col(1);
  const double theta = parameters_(0);
  return (1 + theta * (1 - 2 * u1.array()) * (1 - 2 * u2.array())).matrix();
}

// CDF
inline Eigen::VectorXd
FGMBicop::cdf(const Eigen::MatrixXd& u)
{
  const Eigen::VectorXd& u1 = u.col(0);
  const Eigen::VectorXd& u2 = u.col(1);
  const double theta = parameters_(0);
  return (u1.array() * u2.array() +
          theta * u1.array() * (1 - u1.array()) * u2.array() * (1 - u2.array()))
    .matrix();
}

// h-function (conditional distribution)
inline Eigen::VectorXd
FGMBicop::hfunc1_raw(const Eigen::MatrixXd& u)
{
  const Eigen::VectorXd& u1 = u.col(0);
  const Eigen::VectorXd& u2 = u.col(1);
  const double theta = parameters_(0);
  return (u2.array() +
          theta * (1 - 2 * u1.array()) * u2.array() * (1 - u2.array()))
    .matrix();
}

inline Eigen::VectorXd
FGMBicop::hfunc2_raw(const Eigen::MatrixXd& u)
{
  const Eigen::VectorXd& u1 = u.col(0);
  const Eigen::VectorXd& u2 = u.col(1);
  const double theta = parameters_(0);
  return (u1.array() +
          theta * u1.array() * (1 - u1.array()) * (1 - 2 * u2.array()))
    .matrix();
}

// inverse h-functions (not closed-form for FGM, using an approximation or
// iterative methods here if needed)
inline Eigen::VectorXd
FGMBicop::hinv1_raw(const Eigen::MatrixXd& u)
{
  // (Nelsen book Ex 3.23)
  auto a = 1 + parameters_(0) * (1 - 2 * u.col(0).array());
  auto b = (a.pow(2) - 4 * (a - 1) * u.col(1).array()).sqrt();
  return 2 * u.col(1).array() / (b + a);
}


inline Eigen::VectorXd
FGMBicop::hinv2_raw(const Eigen::MatrixXd& u)
{
  // (Nelsen book Ex 3.23)
  auto a = 1 + parameters_(0) * (1 - 2 * u.col(1).array());
  auto b = (a.pow(2) - 4 * (a - 1) * u.col(0).array()).sqrt();
  return 2 * u.col(0).array() / (b + a);
}

// parameter conversion
inline Eigen::MatrixXd
FGMBicop::tau_to_parameters(const double& tau)
{
  // FGM's Kendall's tau = 2 * theta / 9
  if (std::abs(tau) > 1.0 / 4.5) {
    throw std::runtime_error("|tau| too large for FGM copula.");
  }
  return Eigen::VectorXd::Constant(1, 4.5 * tau);
}

inline double
FGMBicop::parameters_to_tau(const Eigen::MatrixXd& params)
{
  return params(0, 0) / 4.5;
}

// starting parameters
inline Eigen::VectorXd
FGMBicop::get_start_parameters(const double tau)
{
  if (std::abs(tau) <= 1.0 / 4.5) {
    return tau_to_parameters(tau);
  } else {
    auto theta = std::pow(-1.0, static_cast<double>(tau < 0)) * 0.9;
    return Eigen::VectorXd::Constant(1, theta);
  }
}

} // namespace vinecopulib