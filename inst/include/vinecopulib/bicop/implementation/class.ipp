// Copyright © 2016-2021 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <mutex>
#include <vinecopulib/bicop/abstract.hpp>
#include <vinecopulib/bicop/tools_select.hpp>
#include <vinecopulib/misc/tools_interface.hpp>
#include <vinecopulib/misc/tools_serialization.hpp>
#include <vinecopulib/misc/tools_stats.hpp>
#include <vinecopulib/misc/tools_stl.hpp>

//! Tools for bivariate and vine copula modeling
namespace vinecopulib {

//! @brief Instantiates a specific bivariate copula model.
//! @param family The copula family.
//! @param rotation The rotation of the copula; one of 0, 90, 180, or 270
//!     (for Independence, Gaussian, Student, Frank, and nonparametric
//!     families, only 0 is allowed).
//! @param parameters The copula parameters.
//! @param var_types Two strings specifying the types of the variables,
//!   e.g., `("c", "d")` means first variable continuous, second discrete.
inline Bicop::Bicop(const BicopFamily family,
                    const int rotation,
                    const Eigen::MatrixXd& parameters,
                    const std::vector<std::string>& var_types)
{
  bicop_ = AbstractBicop::create(family, parameters);
  // family must be set before checking the rotation
  set_rotation(rotation);
  if (bicop_->get_family() != BicopFamily::indep) {
    bicop_->set_loglik();
  } else {
    bicop_->set_loglik(0.0);
  }
  set_var_types(var_types);
}

//! @brief Instantiates from data.
//!
//! Equivalent to creating a default `Bicop()` and then selecting the model
//! using `select()`.
//!
//! @param data See `select()`.
//! @param controls See `select()`.
//! @param var_types Two strings specifying the types of the variables,
//!   e.g., `("c", "d")` means first variable continuous, second discrete.
inline Bicop::Bicop(const Eigen::MatrixXd& data,
                    const FitControlsBicop& controls,
                    const std::vector<std::string>& var_types)
{
  set_var_types(var_types);
  select(data, controls);
}

//! @brief Copy constructor (deep copy)
//!
//! @param other Bicop object to copy.
inline Bicop::Bicop(const Bicop& other)
  : Bicop(other.get_family(),
          other.get_rotation(),
          other.get_parameters(),
          other.get_var_types())
{
  nobs_ = other.nobs_;
  bicop_->set_loglik(other.bicop_->get_loglik());
  bicop_->set_npars(other.bicop_->get_npars());
}

//! @brief Copy assignment operator (deep copy)
//!
//! @param other Bicop object to copy.
inline Bicop&
Bicop::operator=(Bicop other)
{
  // copy/swap idiom
  std::swap(bicop_, other.bicop_);
  std::swap(rotation_, other.rotation_);
  std::swap(nobs_, other.nobs_);
  std::swap(var_types_, other.var_types_);
  return *this;
}

//! @brief Instantiates from a nlohmann::json object.
//! @param input The nlohmann::json object to convert from
//! (see `to_json()` for the structure of the input).
inline Bicop::Bicop(const nlohmann::json& input)
  : Bicop(get_family_enum(input["family"]),
          static_cast<int>(input["rotation"]),
          tools_serialization::json_to_matrix<double>(input["parameters"]))
{
  // try block for backwards compatibility
  try {
    var_types_ =
      tools_serialization::json_to_vector<std::string>(input["var_types"]);
    nobs_ = static_cast<size_t>(input["nobs_"]);
    bicop_->set_loglik(input["loglik"]);
    bicop_->set_npars(input["npars"]);
  } catch (...) {
  }
}

//! @brief Instantiates from a JSON file.
//!
//! The input file contains four attributes:
//! `"family"`, `"rotation"`, `"parameters"`, `"var_types"` respectively a
//! string for the family name, an integer for the rotation, and a numeric
//! matrix for the parameters, and a list of two strings for the variable
//! types.
//!
//! @param filename The name of the JSON file to read.
inline Bicop::Bicop(const std::string& filename)
  : Bicop(tools_serialization::file_to_json(filename))
{}

//! @brief Convert the copula into a nlohmann::json object.
//!
//! The nlohmann::json is contains of three values named
//! `"family"`, `"rotation"`, `"parameters"`, `"var_types"`,
//! respectively a string for the family name, an integer for the rotation,
//! a numeric matrix for the parameters and a list of two strings for the
//! variables types.
//!
//! @return the nlohmann::json object containing the copula.
inline nlohmann::json
Bicop::to_json() const
{
  nlohmann::json output;
  output["family"] = get_family_name();
  output["rotation"] = rotation_;
  output["parameters"] = tools_serialization::matrix_to_json(get_parameters());
  output["var_types"] = tools_serialization::vector_to_json(var_types_);
  output["nobs_"] = nobs_;
  output["loglik"] = bicop_->get_loglik();
  output["npars"] = bicop_->get_npars();

  return output;
}

//! @brief Write the copula object into a JSON file.
//!
//! The written file contains four attributes:
//! `"family"`, `"rotation"`, `"parameters"`, `"var_types"` respectively a
//! string for the family name, an integer for the rotation, and a numeric
//! matrix for the parameters, and a list of two strings for the variable
//! types.
//!
//! @param filename The name of the file to write.
inline void
Bicop::to_file(const std::string& filename) const
{
  tools_serialization::json_to_file(filename, to_json());
}

//! @brief Evaluates the copula density.
//!
//! The copula density is defined as joint density divided by marginal
//! densities, irrespective of variable types.
//!
//! @param u An \f$ n \times (2 + k) \f$ matrix of observations contained in
//!   \f$(0, 1) \f$, where \f$ k \f$ is the number of discrete variables.
//! @return The copula density evaluated at \c u.
inline Eigen::VectorXd
Bicop::pdf(const Eigen::MatrixXd& u) const
{
  check_data(u);
  return bicop_->pdf(prep_for_abstract(u));
}

//! @brief Evaluates the copula distribution.
//!
//! @param u An \f$ n \times (2 + k) \f$ matrix of observations contained in
//!   \f$(0, 1) \f$, where \f$ k \f$ is the number of discrete variables.
//! @return The copula distribution evaluated at \c u.
inline Eigen::VectorXd
Bicop::cdf(const Eigen::MatrixXd& u) const
{
  check_data(u);
  Eigen::VectorXd p = bicop_->cdf(prep_for_abstract(u).leftCols(2));
  switch (rotation_) {
    default:
      return p;

    case 90:
      return u.col(1) - p;

    case 180:
      return p.array() - 1 + u.leftCols(2).rowwise().sum().array();

    case 270:
      return u.col(0) - p;
  }
}

//! @brief Evaluates the first h-function.
//!
//! The first h-function is
//! \f$ h_1(u_1, u_2) = P(U_2 \le u_2 | U_1 = u_1) \f$.
//! @param u An \f$ n \times (2 + k) \f$ matrix of observations contained in
//!   \f$(0, 1) \f$, where \f$ k \f$ is the number of discrete variables.
inline Eigen::VectorXd
Bicop::hfunc1(const Eigen::MatrixXd& u) const
{
  check_data(u);
  Eigen::VectorXd h(u.rows());
  switch (rotation_) {
    default:
      h = bicop_->hfunc1(prep_for_abstract(u));
      break;

    case 90:
      h = bicop_->hfunc2(prep_for_abstract(u));
      break;

    case 180:
      h = 1.0 - bicop_->hfunc1(prep_for_abstract(u)).array();
      break;

    case 270:
      h = 1.0 - bicop_->hfunc2(prep_for_abstract(u)).array();
      break;
  }
  tools_eigen::trim(h, 0.0, 1.0);
  return h;
}

//! @brief Evaluates the second h-function.
//!
//! The second h-function is
//! \f$ h_2(u_1, u_2) = P(U_1 \le u_1 | U_2 = u_2)  \f$.
//! @param u An \f$ n \times (2 + k) \f$ matrix of observations contained in
//!   \f$(0, 1) \f$, where \f$ k \f$ is the number of discrete variables.
inline Eigen::VectorXd
Bicop::hfunc2(const Eigen::MatrixXd& u) const
{
  check_data(u);
  Eigen::VectorXd h(u.rows());
  switch (rotation_) {
    default:
      h = bicop_->hfunc2(prep_for_abstract(u));
      break;

    case 90:
      h = 1.0 - bicop_->hfunc1(prep_for_abstract(u)).array();
      break;

    case 180:
      h = 1.0 - bicop_->hfunc2(prep_for_abstract(u)).array();
      break;

    case 270:
      h = bicop_->hfunc1(prep_for_abstract(u)).array();
      break;
  }
  tools_eigen::trim(h, 0.0, 1.0);
  return h;
}

//! @brief Evaluates the inverse of the first h-function.
//!
//! The first h-function is
//! \f$ h_1(u_1, u_2) = P(U_2 \le u_2 | U_1 = u_1) \f$.
//! The inverse is calulated w.r.t. the second argument.
//! @param u An \f$ n \times (2 + k) \f$ matrix of observations contained in
//!   \f$(0, 1) \f$, where \f$ k \f$ is the number of discrete variables.
inline Eigen::VectorXd
Bicop::hinv1(const Eigen::MatrixXd& u) const
{
  check_data(u);
  Eigen::VectorXd hi(u.rows());
  switch (rotation_) {
    default:
      hi = bicop_->hinv1(prep_for_abstract(u));
      break;

    case 90:
      hi = bicop_->hinv2(prep_for_abstract(u));
      break;

    case 180:
      hi = 1.0 - bicop_->hinv1(prep_for_abstract(u)).array();
      break;

    case 270:
      hi = 1.0 - bicop_->hinv2(prep_for_abstract(u)).array();
      break;
  }
  tools_eigen::trim(hi, 0.0, 1.0);
  return hi;
}

//! @brief Evaluates the inverse of the second h-function.
//!
//! The second h-function is
//! \f$ h_2(u_1, u_2) = P(U_1 \le u_1 | U_2 = u_2)  \f$.
//! The inverse is calculated w.r.t. the first argument.
//! @param u An \f$ n \times (2 + k) \f$ matrix of observations contained in
//!   \f$(0, 1) \f$, where \f$ k \f$ is the number of discrete variables.
inline Eigen::VectorXd
Bicop::hinv2(const Eigen::MatrixXd& u) const
{
  check_data(u);
  Eigen::VectorXd hi(u.rows());
  switch (rotation_) {
    default:
      hi = bicop_->hinv2(prep_for_abstract(u));
      break;

    case 90:
      hi = 1.0 - bicop_->hinv1(prep_for_abstract(u)).array();
      break;

    case 180:
      hi = 1.0 - bicop_->hinv2(prep_for_abstract(u)).array();
      break;

    case 270:
      hi = bicop_->hinv1(prep_for_abstract(u));
      break;
  }
  tools_eigen::trim(hi, 0.0, 1.0);
  return hi;
}
//! @}

//! @brief Simulates from a bivariate copula.
//!
//! If `qrng = TRUE`, generalized Halton sequences are used.
//! For more information on Generalized Halton sequences, see
//! Faure, H., Lemieux, C. (2009). Generalized Halton Sequences in 2008:
//! A Comparative Study. ACM-TOMACS 19(4), Article 15.
//!
//! @param n Number of observations.
//! @param qrng Set to true for quasi-random numbers.
//! @param seeds Seeds of the (quasi-)random number generator; if empty
//! (default), the (quasi-)random number generator is seeded randomly.
//! @return An \f$ n \times 2 \f$ matrix of samples from the copula model.
inline Eigen::MatrixXd
Bicop::simulate(const size_t& n,
                const bool qrng,
                const std::vector<int>& seeds) const
{
  auto u = tools_stats::simulate_uniform(n, 2, qrng, seeds);
  // use inverse Rosenblatt transform to generate a sample from the copula
  // (always simulate continuous data)
  u.col(1) = this->as_continuous().hinv1(u);
  return u;
}

//! @brief Evaluates the log-likelihood.
//!
//! The log-likelihood is defined as
//! \f[ \mathrm{loglik} = \sum_{i = 1}^n \log c(U_{1, i}, U_{2, i}), \f]
//! where \f$ c \f$ is the copula density, see `Bicop::pdf()`.
//!
//! @param u An \f$ n \times (2 + k) \f$ matrix of observations contained in
//!   \f$(0, 1) \f$, where \f$ k \f$ is the number of discrete variables.
inline double
Bicop::loglik(const Eigen::MatrixXd& u) const
{
  if (u.rows() < 1) {
    return get_loglik();
  } else {
    tools_eigen::check_if_in_unit_cube(u);
    return bicop_->loglik(prep_for_abstract(u));
  }
}

//! @brief Evaluates the Akaike information criterion (AIC).
//!
//! The AIC is defined as
//! \f[ \mathrm{AIC} = -2\, \mathrm{loglik} + 2 p, \f]
//! where \f$ \mathrm{loglik} \f$ is the log-liklihood (see `loglik()`)
//! and \f$ p \f$ is the (effective) number of parameters of the model.
//! The AIC is a consistent model selection criterion even
//! for nonparametric models.
//!
//! @param u An \f$ n \times (2 + k) \f$ matrix of observations contained in
//!   \f$(0, 1) \f$, where \f$ k \f$ is the number of discrete variables.
inline double
Bicop::aic(const Eigen::MatrixXd& u) const
{
  return -2 * loglik(u) + 2 * get_npars();
}

//! @brief Evaluates the Bayesian information criterion (BIC).
//!
//! The BIC is defined as
//! \f[ \mathrm{BIC} = -2\, \mathrm{loglik} +  \log(n) p, \f]
//! where \f$ \mathrm{loglik} \f$ is the log-liklihood (see `loglik()`)
//! and \f$ p \f$ is the (effective) number of parameters of the model.
//! The BIC is a consistent model selection criterion
//! for parametric models.
//!
//! @param u An \f$ n \times (2 + k) \f$ matrix of observations contained in
//!   \f$(0, 1) \f$, where \f$ k \f$ is the number of discrete variables.
inline double
Bicop::bic(const Eigen::MatrixXd& u) const
{
  Eigen::MatrixXd u_no_nan = u;
  double n = static_cast<double>(nobs_);
  if (u.rows() > 0) {
    tools_eigen::remove_nans(u_no_nan);
    n = static_cast<double>(u_no_nan.rows());
  }
  return -2 * loglik(u_no_nan) + get_npars() * log(n);
}

// clang-format off
//! @brief Evaluates the modified Bayesian information criterion (mBIC).
//!
//! The mBIC is defined as
//! \f[ \mathrm{BIC} = -2\, \mathrm{loglik} +  p \log(n) - 2 (I \log(\psi_0) + (1 - I) \log(1 - \psi_0), \f]
//! where \f$ \mathrm{loglik} \f$ is the \log-liklihood
//! (see `loglik()`), \f$ p \f$ is the (effective) number of parameters of the
//! model, and \f$ \psi_0 \f$ is the prior probability of having a
//! non-independence copula and \f$ I \f$ is an indicator for the family being
//! non-independence.
//!
//! @param u An \f$ n \times (2 + k) \f$ matrix of observations contained in
//!   \f$(0, 1) \f$, where \f$ k \f$ is the number of discrete variables.
//! @param psi0 Prior probability of a non-independence copula.
// clang-format on
inline double
Bicop::mbic(const Eigen::MatrixXd& u, const double psi0) const
{
  Eigen::MatrixXd u_no_nan = u;
  bool is_indep = (this->get_family() == BicopFamily::indep);
  double npars = this->get_npars();
  double log_prior = static_cast<double>(!is_indep) * std::log(psi0) +
                     static_cast<double>(is_indep) * std::log(1.0 - psi0);
  double n = static_cast<double>(nobs_);
  if (u.rows() > 0) {
    n = static_cast<double>(u_no_nan.rows());
  }
  return -2 * this->loglik(u_no_nan) + std::log(n) * npars - 2 * log_prior;
}

//! @brief Returns the actual number of parameters for parameteric families.
//!
//! For nonparametric families, there is a conceptually similar definition in
//! the sense that it can be used in the calculation of fit statistics.
inline double
Bicop::get_npars() const
{
  return bicop_->get_npars();
}

//! @brief Converts a Kendall's \f$ \tau \f$ into copula parameters.
//!
//! It only works for one-parameter families.
//! @param tau A value in \f$ (-1, 1) \f$.
inline Eigen::MatrixXd
Bicop::tau_to_parameters(const double& tau) const
{
  return bicop_->tau_to_parameters(tau);
}

//! @brief Converts the copula parameters to Kendall's \f$ tau \f$.
//!
//! @param parameters The parameters (must be a valid parametrization of
//!     the current family).
inline double
Bicop::parameters_to_tau(const Eigen::MatrixXd& parameters) const
{
  double tau = bicop_->parameters_to_tau(parameters);
  if (tools_stl::is_member(rotation_, { 90, 270 })) {
    tau *= -1;
  }
  return tau;
}

//! @name Getters and setters
//!
//! @{

//! @brief Gets the copula family.
inline BicopFamily
Bicop::get_family() const
{
  return bicop_->get_family();
}

//! @brief Gets the copula family as a string.
inline std::string
Bicop::get_family_name() const
{
  return bicop_->get_family_name();
}

//! @brief Gets the rotation.
inline int
Bicop::get_rotation() const
{
  return rotation_;
}

//! @brief Gets the parameters.
inline Eigen::MatrixXd
Bicop::get_parameters() const
{
  return bicop_->get_parameters();
}

//! @brief Gets the log-likelihood (only for fitted objects).
inline double
Bicop::get_loglik() const
{
  check_fitted();
  return bicop_->get_loglik();
}

//! @brief Gets the number of observations (only for fitted objects).
inline size_t
Bicop::get_nobs() const
{
  check_fitted();
  return nobs_;
}

//! @brief Gets the aic (only for fitted objects).
inline double
Bicop::get_aic() const
{
  check_fitted();
  return -2 * bicop_->get_loglik() + 2 * bicop_->get_npars();
}

//! @brief Gets the bic (only for fitted objects).
inline double
Bicop::get_bic() const
{
  check_fitted();
  double npars = bicop_->get_npars();
  return -2 * bicop_->get_loglik() + std::log(nobs_) * npars;
}

//! @brief Gets the modified bic (only for fitted objects).
inline double
Bicop::get_mbic(const double psi0) const
{
  check_fitted();
  return -2 * bicop_->get_loglik() + compute_mbic_penalty(nobs_, psi0);
}

inline double
Bicop::compute_mbic_penalty(const size_t nobs, const double psi0) const
{
  double npars = bicop_->get_npars();
  bool is_indep = (this->get_family() == BicopFamily::indep);
  double log_prior = static_cast<double>(!is_indep) * std::log(psi0) +
                     static_cast<double>(is_indep) * std::log(1.0 - psi0);
  return std::log(nobs) * npars - 2 * log_prior;
}

//! @brief Gets the Kendall's tau.
inline double
Bicop::get_tau() const
{
  return parameters_to_tau(bicop_->get_parameters());
}

//! @brief Sets the rotation.
inline void
Bicop::set_rotation(const int rotation)
{
  check_rotation(rotation);
  if ((rotation_ - rotation % 180) != 0) {
    flip_abstract_var_types();
  }
  rotation_ = rotation;
  bicop_->set_loglik();
}

inline void
Bicop::check_data(const Eigen::MatrixXd& u) const
{
  check_data_dim(u);
  tools_eigen::check_if_in_unit_cube(u);
}

inline void
Bicop::check_data_dim(const Eigen::MatrixXd& u) const
{
  size_t n_cols = u.cols();
  int n_disc = get_n_discrete();
  unsigned short n_cols_exp = static_cast<unsigned short>(2 + n_disc);
  if ((n_cols != n_cols_exp) & (n_cols != 4)) {
    std::stringstream msg;
    msg << "data has wrong number of columns; "
        << "expected: " << n_cols_exp << " or 4, actual: " << n_cols
        << " (model contains ";
    if (n_disc == 0) {
      msg << "no discrete variables)." << std::endl;
    } else if (n_disc == 1) {
      msg << "1 discrete variable)." << std::endl;
    } else {
      msg << get_n_discrete() << " discrete variables)." << std::endl;
    }
    throw std::runtime_error(msg.str());
  }
}

inline void
Bicop::flip_abstract_var_types()
{
  std::swap(bicop_->var_types_[0], bicop_->var_types_[1]);
}

inline void
Bicop::set_parameters(const Eigen::MatrixXd& parameters)
{
  bicop_->set_parameters(parameters);
  bicop_->set_loglik();
}

//! @brief Sets variable types.
//! @param var_types A vector of size two specifying the types of the variables,
//!   e.g., `{"c", "d"}` means first variable continuous, second discrete.
inline void
Bicop::set_var_types(const std::vector<std::string>& var_types)
{
  check_var_types(var_types);
  var_types_ = var_types;
  if (bicop_) {
    bicop_->set_var_types(var_types);
    if (tools_stl::is_member(static_cast<size_t>(rotation_), { 90, 270 })) {
      flip_abstract_var_types();
    }
  }
}

//! @brief Gets variable types.
inline std::vector<std::string>
Bicop::get_var_types() const
{
  return var_types_;
}
//! @}

//! @name Utilities
//! @{
//! useful functions for bivariate copulas

//! Adjusts the copula model to a change in the variable order.
inline void
Bicop::flip()
{
  BicopFamily family = bicop_->get_family();
  // change internal representation
  if (tools_stl::is_member(family, bicop_families::flip_by_rotation)) {
    double loglik = bicop_->get_loglik();
    if (rotation_ == 90) {
      set_rotation(270);
    } else if (rotation_ == 270) {
      set_rotation(90);
    }
    bicop_->set_loglik(loglik);
  } else {
    flip_abstract_var_types();
    bicop_->flip();
  }
  // change Bicop-level var_types
  std::swap(var_types_[0], var_types_[1]);
}

//! @brief Summarizes the model into a string (can be used for printing).
inline std::string
Bicop::str() const
{
  std::stringstream bicop_str;
  bicop_str << get_family_name();
  if (get_rotation() != 0) {
    bicop_str << " " << get_rotation() << "°";
  }
  if (get_family() == BicopFamily::tll) {
    bicop_str << ", parameters = [30x30 grid]";
  } else if (get_family() != BicopFamily::indep) {
    bicop_str << ", parameters = " << get_parameters();
  }
  return bicop_str.str().c_str();
}

//! @brief Gets lower bounds for copula parameters.
inline Eigen::MatrixXd
Bicop::get_parameters_lower_bounds() const
{
  return bicop_->get_parameters_lower_bounds();
}

//! @brief Gets upper bounds for copula parameters.
inline Eigen::MatrixXd
Bicop::get_parameters_upper_bounds() const
{
  return bicop_->get_parameters_upper_bounds();
}

//! @}

inline BicopPtr
Bicop::get_bicop() const
{
  return bicop_;
}

inline Bicop
Bicop::as_continuous() const
{
  std::vector<std::string> cc = { "c", "c" };
  if (var_types_ == cc)
    return *this;
  auto bc_new = *this;
  bc_new.set_var_types(cc);
  return bc_new;
}

//! Fits a bivariate copula (with fixed family) to data.
//!
//! For parametric models, two different methods are available. `"mle"` fits
//! the parameters by maximum-likelihood. `"itau"` uses inversion of
//! Kendall's \f$ \tau \f$, but is only available for one-parameter families
//! and the Student t copula. For the latter, there is a one-to-one
//! transformation for the first parameter, the second is found by profile
//! likelihood optimization (with accuracy of at least 0.5). Nonparametric
//! families have specialized methods, no specification is required.
//!
//! @details When at least one variable is discrete, two types of "observations"
//! are required: the first \f$ n \times 2 \f$ block contains realizations of
//! \f$ F_{X_1}(X_1), F_{X_2}(X_2) \f$. Let \f$ k \f$ denote the number of
//! discrete variables (either one or two). Then the second \f$ n \times k \f$
//! block contains realizations of \f$ F_{X_k}(X_k^-) \f$. The minus indicates a
//! left-sided limit of the cdf. For continuous variables the left limit and the
//! cdf itself coincide. For, e.g., an integer-valued variable, it holds \f$
//! F_{X_k}(X_k^-) = F_{X_k}(X_k - 1) \f$.
//!
//! @param data An \f$ n \times (2 + k) \f$ matrix of observations contained in
//!   \f$(0, 1) \f$, where \f$ k \f$ is the number of discrete variables.
//! @param controls The controls (see FitControlsBicop).
inline void
Bicop::fit(const Eigen::MatrixXd& data, const FitControlsBicop& controls)
{
  std::string method;
  if (tools_stl::is_member(bicop_->get_family(), bicop_families::parametric)) {
    method = controls.get_parametric_method();
  } else {
    method = controls.get_nonparametric_method();
  }
  tools_eigen::check_if_in_unit_cube(data);

  auto w = controls.get_weights();
  Eigen::MatrixXd data_no_nan = data;
  check_weights_size(w, data);
  tools_eigen::remove_nans(data_no_nan, w);

  bicop_->fit(prep_for_abstract(data_no_nan),
              method,
              controls.get_nonparametric_mult(),
              w);
  nobs_ = data_no_nan.rows();
}

//

//! @brief Selects the best fitting model.
//!
//! The function calls `fit()` for all families in
//! `family_set` and selecting the best fitting model by either BIC or AIC,
//! see `bic()` and `aic()`.
//!
//! @details When at least one variable is discrete, two types of "observations"
//! are required: the first \f$ n \times 2 \f$ block contains realizations of
//! \f$ F_{X_1}(X_1), F_{X_2}(X_2) \f$. Let \f$ k \f$ denote the number of
//! discrete variables (either one or two). Then the second \f$ n \times k \f$
//! block contains realizations of \f$ F_{X_k}(X_k^-) \f$. The minus indicates a
//! left-sided limit of the cdf. For continuous variables the left limit and the
//! cdf itself coincide. For, e.g., an integer-valued variable, it holds \f$
//! F_{X_k}(X_k^-) = F_{X_k}(X_k - 1) \f$.
//!
//! @param data An \f$ n \times (2 + k) \f$ matrix of observations contained in
//!   \f$(0, 1) \f$, where \f$ k \f$ is the number of discrete variables.
//! @param controls The controls (see FitControlsBicop).
inline void
Bicop::select(const Eigen::MatrixXd& data, FitControlsBicop controls)
{
  using namespace tools_select;
  check_weights_size(controls.get_weights(), data);
  Eigen::MatrixXd data_no_nan = data;
  {
    auto w = controls.get_weights();
    tools_eigen::remove_nans(data_no_nan, w);
    controls.set_weights(w);
  }
  check_data(data_no_nan);
  nobs_ = data_no_nan.rows();

  bicop_ = AbstractBicop::create();
  bicop_->set_var_types(var_types_);
  rotation_ = 0;
  bicop_->set_loglik(0.0);
  if (data_no_nan.rows() >= 10) {
    tools_eigen::trim(data_no_nan);
    std::vector<Bicop> bicops = create_candidate_bicops(data_no_nan, controls);
    for (auto& bc : bicops) {
      bc.set_var_types(var_types_);
    }

    // Estimate all models and select the best one using the
    // selection_criterion
    double fitted_criterion = std::numeric_limits<double>::max();
    std::mutex m;
    auto fit_and_compare = [&](Bicop cop) {
      tools_interface::check_user_interrupt();
      // Estimate the model
      cop.fit(data_no_nan, controls);

      // Compute the selection criterion
      double new_criterion;
      double ll = cop.get_loglik();
      if (controls.get_selection_criterion() == "loglik") {
        new_criterion = -ll;
      } else if (controls.get_selection_criterion() == "aic") {
        new_criterion = -2 * ll + 2 * cop.get_npars();
      } else {
        double n_eff = static_cast<double>(data_no_nan.rows());
        if (controls.get_weights().size() > 0) {
          n_eff = std::pow(controls.get_weights().sum(), 2);
          n_eff /= controls.get_weights().array().pow(2).sum();
        }
        double npars = cop.get_npars();

        new_criterion = -2 * ll + log(n_eff) * npars; // BIC
        if (controls.get_selection_criterion() == "mbic") {
          // correction for mBIC
          bool is_indep = (cop.get_family() == BicopFamily::indep);
          double psi0 = controls.get_psi0();
          double log_prior = static_cast<double>(!is_indep) * log(psi0) +
                             static_cast<double>(is_indep) * log(1.0 - psi0);
          new_criterion -= 2 * log_prior;
        }
      }

      // the following block modifies thread-external variables
      // and is thus shielded by a mutex
      {
        std::lock_guard<std::mutex> lk(m);
        // If the new model is better than the current one,
        // then replace the current model by the new one
        if (new_criterion < fitted_criterion) {
          fitted_criterion = new_criterion;
          bicop_ = cop.get_bicop();
          rotation_ = cop.get_rotation();
        }
      }
    };

    tools_thread::ThreadPool pool(controls.get_num_threads());
    pool.map(fit_and_compare, bicops);
  }
}

//! @brief Adds an additional column if there's only one discrete variable;
//! removes superfluous columns for continuous variables.
//! (continuous models only require two columns, discrete models always four)
inline Eigen::MatrixXd
Bicop::format_data(const Eigen::MatrixXd& u) const
{
  auto n_disc = get_n_discrete();
  if (n_disc == 0) {
    return u.leftCols(2);
  } else if (n_disc == 2) {
    return u;
  }
  // n_disc = 1:
  Eigen::MatrixXd u_new(u.rows(), 4);
  u_new.leftCols(2) = u.leftCols(2);
  int disc_col = (var_types_[1] == "d");
  int cont_col = 1 - disc_col;
  // We already know that there is one discrete and one continuous variable. Now
  // there are two cases:
  // 1. `u.cols() == 3`: then the F(x^-) values for the discrete variable is
  // always in the last column, i.e. `u.col(2)`.
  // 2. `u.cols() == 4`: Then the F(x^-) values for the discrete variable is in
  // the third column if variable 1 is discrete, and in the fourth column if
  // variable 2 is discrete. Thus, `u.col(2 + disc_col)`.
  int old_disc_col = 2 + (u.cols() == 4) * disc_col;
  u_new.col(2 + disc_col) = u.col(old_disc_col);
  u_new.col(2 + cont_col) = u.col(cont_col);
  return u_new;
}

//! @brief Rotates the data corresponding to the models rotation.
//! @param u An `n x 2` matrix.
inline void
Bicop::rotate_data(Eigen::MatrixXd& u) const
{
  // counter-clockwise rotations
  switch (rotation_) {
    case 0:
      break;

    case 90:
      u.col(0).swap(u.col(1));
      u.col(1) = 1 - u.col(1).array();
      if (u.cols() == 4) {
        u.col(2).swap(u.col(3));
        u.col(3) = 1 - u.col(3).array();
      }
      break;

    case 180:
      u = 1 - u.array();
      break;

    case 270:
      u.col(0).swap(u.col(1));
      u.col(0) = 1 - u.col(0).array();
      if (u.cols() == 4) {
        u.col(2).swap(u.col(3));
        u.col(2) = 1 - u.col(2).array();
      }
      break;
  }
}

//! @brief Prepares data for use with the `AbstractBicop` class:
//! - add an additional column if there's only one discrete variable.
//! - trim the data to the interval [1e-10, 1 - 1e-10] for numerical stability.
//! - rotate the data appropriately (`AbstractBicop` is always 0deg-rotation).
inline Eigen::MatrixXd
Bicop::prep_for_abstract(const Eigen::MatrixXd& u) const
{
  auto u_new = format_data(u);
  tools_eigen::trim(u_new);
  rotate_data(u_new);
  return u_new;
}

//! @brief Checks whether the supplied rotation is valid (only 0, 90, 180, 270
//! allowd).
inline void
Bicop::check_rotation(int rotation) const
{
  using namespace tools_stl;
  std::vector<int> allowed_rotations = { 0, 90, 180, 270 };
  if (!is_member(rotation, allowed_rotations)) {
    throw std::runtime_error("rotation must be one of {0, 90, 180, 270}");
  }
  if (is_member(bicop_->get_family(), bicop_families::rotationless)) {
    if (rotation != 0) {
      throw std::runtime_error("rotation must be 0 for the " +
                               bicop_->get_family_name() + " copula");
    }
  }
}

//! @brief Checks whether weights and data have matching sizes.
inline void
Bicop::check_weights_size(const Eigen::VectorXd& weights,
                          const Eigen::MatrixXd& data) const
{
  if ((weights.size() > 0) & (weights.size() != data.rows())) {
    throw std::runtime_error("sizes of weights and data don't match.");
  }
}

//! @brief Checks whether the Bicop object was fitted to data.
inline void
Bicop::check_fitted() const
{
  if ((boost::math::isnan)(bicop_->get_loglik())) {
    throw std::runtime_error("copula has not been fitted from data or its "
                             "parameters have been modified manually");
  }
}

//! @brief Checks whether var_types have the correct length and are either "c"
//! or "d".
inline void
Bicop::check_var_types(const std::vector<std::string>& var_types) const
{
  if (var_types.size() != 2) {
    throw std::runtime_error("var_types must have size two.");
  }
  for (auto t : var_types) {
    if (!tools_stl::is_member(t, { "c", "d" })) {
      throw std::runtime_error("var type must be either 'c' or 'd'.");
    }
  }
}

//! @brief Returns the number of discrete variables.
inline unsigned short
Bicop::get_n_discrete() const
{
  int n_discrete = 0;
  for (auto t : var_types_) {
    n_discrete += (t == "d");
  }
  return static_cast<unsigned short>(n_discrete);
}
}
