// Copyright © 2016-2021 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <memory>

#include <Eigen/Dense>
#include <vinecopulib/bicop/family.hpp>

namespace vinecopulib {
//! @brief An abstract class for bivariate copula families.
//!
//! This class is used in the implementation underlying the Bicop class.
//! Users should not use AbstractBicop or derived classes directly, but
//! always work with the Bicop interface.
class AbstractBicop
{
  friend class Bicop;

public:
  virtual ~AbstractBicop() = 0;

protected:
  // Factories
  static std::shared_ptr<AbstractBicop> create(
    BicopFamily family = BicopFamily::indep,
    const Eigen::MatrixXd& parameters = Eigen::MatrixXd());

  // Getters and setters
  BicopFamily get_family() const;

  std::string get_family_name() const;

  double get_loglik() const;

  void set_loglik(const double loglik = NAN);

  void set_var_types(const std::vector<std::string>& var_types);

  virtual Eigen::MatrixXd get_parameters() const = 0;

  virtual Eigen::MatrixXd get_parameters_lower_bounds() const = 0;

  virtual Eigen::MatrixXd get_parameters_upper_bounds() const = 0;

  virtual void set_parameters(const Eigen::MatrixXd& parameters) = 0;

  // Virtual methods
  virtual void fit(const Eigen::MatrixXd& data,
                   std::string method,
                   double mult,
                   const Eigen::VectorXd& weights) = 0;

  virtual double get_npars() const = 0;

  virtual void set_npars(const double& npars) = 0;

  virtual double parameters_to_tau(const Eigen::MatrixXd& parameters) = 0;

  virtual void flip() = 0;

  // following are virtual so they can be overriden by KernelBicop
  virtual Eigen::VectorXd pdf(const Eigen::MatrixXd& u);

  virtual Eigen::VectorXd cdf(const Eigen::MatrixXd& u) = 0;

  virtual Eigen::VectorXd hfunc1(const Eigen::MatrixXd& u);

  virtual Eigen::VectorXd hfunc2(const Eigen::MatrixXd& u);

  Eigen::VectorXd hinv1(const Eigen::MatrixXd& u);

  Eigen::VectorXd hinv2(const Eigen::MatrixXd& u);

  virtual Eigen::VectorXd pdf_raw(const Eigen::MatrixXd& u) = 0;

  virtual Eigen::VectorXd hfunc1_raw(const Eigen::MatrixXd& u) = 0;

  virtual Eigen::VectorXd hfunc2_raw(const Eigen::MatrixXd& u) = 0;

  virtual Eigen::VectorXd hinv1_raw(const Eigen::MatrixXd& u) = 0;

  virtual Eigen::VectorXd hinv2_raw(const Eigen::MatrixXd& u) = 0;

  virtual Eigen::MatrixXd tau_to_parameters(const double& tau) = 0;
  Eigen::MatrixXd no_tau_to_parameters(const double&);

  // Misc methods
  Eigen::VectorXd hinv1_num(const Eigen::MatrixXd& u);

  Eigen::VectorXd hinv2_num(const Eigen::MatrixXd& u);

  Eigen::VectorXd pdf_c_d(const Eigen::MatrixXd& u);

  Eigen::VectorXd pdf_d_d(const Eigen::MatrixXd& u);

  double loglik(const Eigen::MatrixXd& u,
                const Eigen::VectorXd weights = Eigen::VectorXd());

  // Data members
  BicopFamily family_;
  double loglik_{ NAN };
  std::vector<std::string> var_types_{ "c", "c" };
};

//! A shared pointer to an object of class AbstracBicop.
typedef std::shared_ptr<AbstractBicop> BicopPtr;
}

#include <vinecopulib/bicop/implementation/abstract.ipp>
