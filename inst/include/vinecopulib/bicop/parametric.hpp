// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <functional>
#include <vinecopulib/bicop/abstract.hpp>

namespace vinecopulib {

//! @brief An abstract class for parametric copula families.
//!
//! This class is used in the implementation underlying the Bicop class.
//! Users should not use AbstractBicop or derived classes directly, but
//! always work with the Bicop interface.
class ParBicop : public AbstractBicop
{
protected:
  // Getters and setters
  Eigen::MatrixXd get_parameters() const;

  Eigen::MatrixXd get_parameters_lower_bounds() const;

  Eigen::MatrixXd get_parameters_upper_bounds() const;

  void set_parameters(const Eigen::MatrixXd& parameters);

  void flip();

  // Data members
  Eigen::MatrixXd parameters_;
  Eigen::MatrixXd parameters_lower_bounds_;
  Eigen::MatrixXd parameters_upper_bounds_;

  void fit(const Eigen::MatrixXd& data,
           std::string method,
           double,
           size_t,
           const Eigen::VectorXd& weights);

  double get_npars() const;

  void set_npars(const double& npars);

  virtual Eigen::VectorXd get_start_parameters(const double tau) = 0;

  // fallback derivative leaves: central finite differences of the value
  // leaves, so that every parametric family supports the derivative
  // interface; families with closed forms override these.
  Eigen::VectorXd pdf_deriv_raw(const Eigen::MatrixXd& u,
                                const Eigen::MatrixXd& parameters,
                                const std::string& deriv);

  Eigen::VectorXd pdf_deriv2_raw(const Eigen::MatrixXd& u,
                                 const Eigen::MatrixXd& parameters,
                                 const std::string& deriv);

  Eigen::VectorXd hfunc1_deriv_raw(const Eigen::MatrixXd& u,
                                   const Eigen::MatrixXd& parameters,
                                   const std::string& deriv);

  Eigen::VectorXd hfunc1_deriv2_raw(const Eigen::MatrixXd& u,
                                    const Eigen::MatrixXd& parameters,
                                    const std::string& deriv);

  Eigen::VectorXd hfunc2_deriv_raw(const Eigen::MatrixXd& u,
                                   const Eigen::MatrixXd& parameters,
                                   const std::string& deriv);

  Eigen::VectorXd hfunc2_deriv2_raw(const Eigen::MatrixXd& u,
                                    const Eigen::MatrixXd& parameters,
                                    const std::string& deriv);

private:
  Eigen::VectorXd fd_deriv(
    const std::function<Eigen::VectorXd(const Eigen::MatrixXd&,
                                        const Eigen::MatrixXd&)>& f,
    const Eigen::MatrixXd& u,
    const Eigen::MatrixXd& parameters,
    int comp);

  // central finite differences of a scalar objective `f` w.r.t. each
  // optimization variable, used as the gradient fallback when no analytic
  // score is available (e.g. discrete data). Steps are clipped to `[lb, ub]`.
  Eigen::VectorXd fd_grad(
    const std::function<double(const Eigen::VectorXd&)>& f,
    const Eigen::VectorXd& x,
    const Eigen::VectorXd& lb,
    const Eigen::VectorXd& ub);

  double winsorize_tau(double tau) const;

  void adjust_parameters_bounds(Eigen::MatrixXd& lb,
                                Eigen::MatrixXd& ub,
                                const double& tau,
                                const std::string& method);

  void check_parameters(const Eigen::MatrixXd& parameters);

  void check_parameters_size(const Eigen::MatrixXd& parameters);

  void check_parameters_upper(const Eigen::MatrixXd& parameters);

  void check_parameters_lower(const Eigen::MatrixXd& parameters);

  void check_fit_method(const std::string& method);
};
}

#include <vinecopulib/bicop/implementation/parametric.ipp>
