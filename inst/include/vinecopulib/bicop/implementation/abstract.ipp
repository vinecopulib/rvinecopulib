// Copyright © 2016-2021 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <stdexcept>

#include <vinecopulib/bicop/bb1.hpp>
#include <vinecopulib/bicop/bb6.hpp>
#include <vinecopulib/bicop/bb7.hpp>
#include <vinecopulib/bicop/bb8.hpp>
#include <vinecopulib/bicop/clayton.hpp>
#include <vinecopulib/bicop/frank.hpp>
#include <vinecopulib/bicop/gaussian.hpp>
#include <vinecopulib/bicop/gumbel.hpp>
#include <vinecopulib/bicop/indep.hpp>
#include <vinecopulib/bicop/joe.hpp>
#include <vinecopulib/bicop/student.hpp>
#include <vinecopulib/bicop/tll.hpp>
#include <vinecopulib/misc/tools_eigen.hpp>

namespace vinecopulib {

//! virtual destructor
inline AbstractBicop::~AbstractBicop() {}

//! Instantiates a bivariate copula using the default contructor
//!
//! @param family The copula family.
//! @param parameters The copula parameters (optional, must be compatible
//!     with family).
//! @return A pointer to an object that inherits from AbstractBicop.
//! @{
inline BicopPtr
AbstractBicop::create(BicopFamily family, const Eigen::MatrixXd& parameters)
{
  BicopPtr new_bicop;
  switch (family) {
    case BicopFamily::indep:
      new_bicop = BicopPtr(new IndepBicop());
      break;
    case BicopFamily::gaussian:
      new_bicop = BicopPtr(new GaussianBicop());
      break;
    case BicopFamily::student:
      new_bicop = BicopPtr(new StudentBicop());
      break;
    case BicopFamily::clayton:
      new_bicop = BicopPtr(new ClaytonBicop());
      break;
    case BicopFamily::gumbel:
      new_bicop = BicopPtr(new GumbelBicop());
      break;
    case BicopFamily::frank:
      new_bicop = BicopPtr(new FrankBicop());
      break;
    case BicopFamily::joe:
      new_bicop = BicopPtr(new JoeBicop());
      break;
    case BicopFamily::bb1:
      new_bicop = BicopPtr(new Bb1Bicop());
      break;
    case BicopFamily::bb6:
      new_bicop = BicopPtr(new Bb6Bicop());
      break;
    case BicopFamily::bb7:
      new_bicop = BicopPtr(new Bb7Bicop());
      break;
    case BicopFamily::bb8:
      new_bicop = BicopPtr(new Bb8Bicop());
      break;
    case BicopFamily::tll:
      new_bicop = BicopPtr(new TllBicop());
      break;

    default:
      throw std::runtime_error(std::string("Family not implemented"));
  }

  if (parameters.size() > 0) {
    new_bicop->set_parameters(parameters);
  }

  return new_bicop;
}

//!@}

inline Eigen::MatrixXd
AbstractBicop::no_tau_to_parameters(const double&)
{
  throw std::runtime_error("Method not implemented for this family");
}

//! Getters and setters.
//! @{
inline BicopFamily
AbstractBicop::get_family() const
{
  return family_;
}

inline std::string
AbstractBicop::get_family_name() const
{
  return vinecopulib::get_family_name(family_);
}

inline double
AbstractBicop::get_loglik() const
{
  return loglik_;
}

inline void
AbstractBicop::set_loglik(const double loglik)
{
  loglik_ = loglik;
}

inline void
AbstractBicop::set_var_types(const std::vector<std::string>& var_types)
{
  if (var_types.size() != 2) {
    throw std::runtime_error("var_types must have size two.");
  }
  var_types_ = var_types;
}
//! @}

//! evaluates the pdf, but truncates it's value by DBL_MIN and DBL_MAX.
//! @param u Matrix of evaluation points.
inline Eigen::VectorXd
AbstractBicop::pdf(const Eigen::MatrixXd& u)
{

  Eigen::VectorXd pdf(u.rows());
  if (var_types_ == std::vector<std::string>{ "c", "c" }) {
    pdf = pdf_raw(u.leftCols(2));
  } else if (var_types_ == std::vector<std::string>{ "d", "d" }) {
    pdf = pdf_d_d(u);
  } else {
    pdf = pdf_c_d(u);
  }
  tools_eigen::trim(pdf, DBL_MIN, DBL_MAX);
  return pdf;
}

inline Eigen::VectorXd
AbstractBicop::pdf_c_d(const Eigen::MatrixXd& u)
{
  if (var_types_[0] != "c") {
    return (hfunc2_raw(u.leftCols(2)) - hfunc2_raw(u.rightCols(2)))
      .cwiseQuotient(u.col(0) - u.col(2))
      .cwiseAbs();
  } else {
    return (hfunc1_raw(u.leftCols(2)) - hfunc1_raw(u.rightCols(2)))
      .cwiseQuotient(u.col(1) - u.col(3))
      .cwiseAbs();
  }
}

inline Eigen::VectorXd
AbstractBicop::pdf_d_d(const Eigen::MatrixXd& u)
{
  Eigen::MatrixXd umax = u.leftCols(2);
  Eigen::MatrixXd umin = u.rightCols(2);
  Eigen::VectorXd pdf = cdf(umax) + cdf(umin);
  umax.col(0).swap(umin.col(0));
  pdf -= cdf(umax) + cdf(umin);
  pdf = pdf.array() / (u.col(0) - u.col(2)).array();
  pdf = pdf.array() / (u.col(1) - u.col(3)).array();
  return pdf;
}

inline Eigen::VectorXd
AbstractBicop::hfunc1(const Eigen::MatrixXd& u)
{
  if (var_types_[0] == "d") {
    auto uu = u;
    uu.col(3) = uu.col(1);
    return ((cdf(uu.leftCols(2)) - cdf(uu.rightCols(2))).array() /
            (uu.col(0) - uu.col(2)).array())
      .abs();
  } else {
    return hfunc1_raw(u.leftCols(2));
  }
}

inline Eigen::VectorXd
AbstractBicop::hfunc2(const Eigen::MatrixXd& u)
{
  if (var_types_[1] == "d") {
    auto uu = u;
    uu.col(2) = uu.col(0);
    return ((cdf(uu.leftCols(2)) - cdf(uu.rightCols(2))).array() /
            (uu.col(1) - uu.col(3)).array())
      .abs();
  } else {
    return hfunc2_raw(u.leftCols(2));
  }
}

inline Eigen::VectorXd
AbstractBicop::hinv1(const Eigen::MatrixXd& u)
{
  if (var_types_[0] == "c") {
    return hinv1_raw(u.leftCols(2));
  } else {
    return hinv1_num(u);
  }
}

inline Eigen::VectorXd
AbstractBicop::hinv2(const Eigen::MatrixXd& u)
{
  if (var_types_[1] == "c") {
    return hinv2_raw(u.leftCols(2));
  } else {
    return hinv2_num(u);
  }
}

//! evaluates the log-likelihood.
//! @param u Data matrix.
//! @param weights Optional weights for each observation.
inline double
AbstractBicop::loglik(const Eigen::MatrixXd& u, const Eigen::VectorXd weights)
{
  Eigen::MatrixXd log_pdf = this->pdf(u).array().log();
  if (weights.size() > 0) {
    log_pdf = log_pdf.cwiseProduct(weights);
  }
  tools_eigen::remove_nans(log_pdf);
  return log_pdf.sum();
}

//! Numerical inversion of h-functions
//!
//! These are generic functions to invert the hfunctions numerically.
//! They can be used in derived classes to define \c hinv1 and \c hinv2.
//!
//! @param u \f$m \times 2\f$ matrix of evaluation points.
//! @return The numerical inverse of h-functions.
//! @{
inline Eigen::VectorXd
AbstractBicop::hinv1_num(const Eigen::MatrixXd& u)
{
  Eigen::MatrixXd u_new = u;
  auto h1 = [&](const Eigen::VectorXd& v) {
    u_new.col(1) = v;
    return hfunc1(u_new);
  };

  return tools_eigen::invert_f(u.col(1), h1);
}

inline Eigen::VectorXd
AbstractBicop::hinv2_num(const Eigen::MatrixXd& u)
{
  Eigen::MatrixXd u_new = u;
  auto h1 = [&](const Eigen::VectorXd& x) {
    u_new.col(0) = x;
    return hfunc2(u_new);
  };

  return tools_eigen::invert_f(u.col(0), h1);
}
//! @}
}
