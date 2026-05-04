// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/misc/tools_stats.hpp>

namespace vinecopulib {
inline XtdGumbelBicop::XtdGumbelBicop()
{
  family_ = BicopFamily::xtd_gumbel;;
  parameters_ = Eigen::VectorXd(1);
  parameters_lower_bounds_ = Eigen::VectorXd(1);
  parameters_upper_bounds_ = Eigen::VectorXd(1);
  parameters_ << 0;
  parameters_lower_bounds_ << -20;
  parameters_upper_bounds_ << 20;
  gumbel_bicop_ = std::make_shared<GumbelBicop>();
}

inline void XtdGumbelBicop::set_parameters(const Eigen::MatrixXd& par) 
{
  parameters_ = par;
  gumbel_bicop_->set_parameters(par.array().abs() + 1.0);
}

inline Eigen::VectorXd
XtdGumbelBicop::pdf_raw(const Eigen::MatrixXd& u)
{
  if (parameters_(0) >= 0) {
      return gumbel_bicop_->pdf_raw(u);
  } else {
      return gumbel_bicop_->pdf_raw(tools_stats::rotate_data(u, 90));
  }
}

inline Eigen::VectorXd
XtdGumbelBicop::cdf(const Eigen::MatrixXd& u)
{
  throw std::runtime_error(
    "XtdGumbelBicop does not support the cdf function.");
}

inline Eigen::VectorXd
XtdGumbelBicop::hfunc1_raw(const Eigen::MatrixXd& u)
{
  if (parameters_(0) >= 0) {
    return gumbel_bicop_->hfunc1_raw(u);
  } else {
    return gumbel_bicop_->hfunc1_raw(tools_eigen::swap_cols(tools_stats::rotate_data(u, 90)));
  }
}


inline Eigen::VectorXd
XtdGumbelBicop::hfunc2_raw(const Eigen::MatrixXd& u)
{
  if (parameters_(0) >= 0) {
    return gumbel_bicop_->hfunc1_raw(tools_eigen::swap_cols(u));
  } else {
    return gumbel_bicop_->hfunc1_raw(tools_stats::rotate_data(u, 90));
  }
}

inline Eigen::VectorXd
XtdGumbelBicop::hinv1_raw(const Eigen::MatrixXd& u)
{
  if (parameters_(0) >= 0) {
    return gumbel_bicop_->hinv1_raw(u);
  } else {
    return gumbel_bicop_->hinv1_raw(tools_eigen::swap_cols(tools_stats::rotate_data(u, 90)));
  }
}

inline Eigen::VectorXd
XtdGumbelBicop::hinv2_raw(const Eigen::MatrixXd& u)
{
  if (parameters_(0) >= 0) {
    return gumbel_bicop_->hinv1_raw(tools_eigen::swap_cols(u));
  } else {
    return 1 - gumbel_bicop_->hinv1_raw(tools_stats::rotate_data(u, 90)).array();
  }
}


inline Eigen::VectorXd
XtdGumbelBicop::get_start_parameters(const double tau)
{
  return tau_to_parameters(tau);
}

inline Eigen::MatrixXd
XtdGumbelBicop::tau_to_parameters(const double& tau)
{
  auto par = gumbel_bicop_->tau_to_parameters(std::abs(tau));
  if (tau >= 0) {
    par(0) =  par(0) - 1.0;
  } else {
    par(0) = -par(0) + 1.0;
  }
  return par;
}

inline double
XtdGumbelBicop::parameters_to_tau(const Eigen::MatrixXd& par)
{
  if (par(0) >= 0) {
    return gumbel_bicop_->parameters_to_tau(par.array() + 1.0);
  } else {
    return -gumbel_bicop_->parameters_to_tau(-par.array() + 1.0);
  }
}

}
