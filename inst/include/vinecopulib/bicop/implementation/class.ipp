// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <functional>
#include <mutex>
#include <vinecopulib/bicop/abstract.hpp>
#include <vinecopulib/bicop/tools_select.hpp>
#include <vinecopulib/misc/tools_batch.hpp>
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
//! @details Equivalent to creating a default `Bicop()` and then selecting
//!  the model using `Bicop::select()`.
//!
//! @param data See `Bicop::select()`.
//! @param controls See `Bicop::select()`.
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

//! @brief Instantiates from a `nlohmann::json` object.
//! @param input The `nlohmann::json` object to convert from
//! (see `to_json()` for the structure of the input).
inline Bicop::Bicop(const nlohmann::json& input)
  : Bicop(get_family_enum(input["fam"]),
          static_cast<int>(input["rot"]),
          tools_serialization::json_to_matrix<double>(input["par"]))
{
  // try block for backwards compatibility
  try {
    var_types_ = tools_serialization::json_to_vector<std::string>(input["vt"]);
    nobs_ = static_cast<size_t>(input["nobs"]);
    bicop_->set_loglik(input["ll"]);
    bicop_->set_npars(input["npars"]);
  } catch (...) {
  }
}

//! @brief Instantiates from a JSON or CBOR file.
//!
//! @details Files ending in `.cbor` are read as CBOR. All other filenames are
//! read as JSON for backwards compatibility. The input contains four
//! attributes:
//! `"fam"`, `"rot"`, `"par"`, `"vt"` respectively a
//! string for the family name, an integer for the rotation, and a numeric
//! matrix for the parameters, and a list of two strings for the variable
//! types.
//!
//! @param filename The name of the file to read.
inline Bicop::Bicop(const std::string& filename)
  : Bicop(tools_serialization::file_to_json(filename))
{
}

//! @brief Convert the copula into a nlohmann::json object.
//!
//! @details The `nlohmann::json` is contains of three values named
//! `"fam"`, `"rot"`, `"par"`, `"vt"`,
//! respectively a string for the family name, an integer for the rotation,
//! a numeric matrix for the parameters and a list of two strings for the
//! variables types.
//!
//! @return The `nlohmann::json` object containing the copula.
inline nlohmann::json
Bicop::to_json() const
{
  nlohmann::json output;
  output["fam"] = get_family_name();
  output["rot"] = rotation_;
  output["par"] = tools_serialization::matrix_to_json(get_parameters());
  output["vt"] = tools_serialization::vector_to_json(var_types_);
  output["nobs"] = nobs_;
  output["ll"] = bicop_->get_loglik();
  output["npars"] = bicop_->get_npars();

  return output;
}

//! @brief Writes the copula object into a JSON or CBOR file.
//!
//! @details Filenames ending in `.cbor` are written as CBOR. All other
//! filenames are written as JSON for backwards compatibility. The written
//! representation contains four attributes:
//! `"fam"`, `"rot"`, `"par"`, `"vt"`, `"nobs"`, `"ll"`, `"npars"`
//! respectively a string for the family name, an integer for the rotation, and
//! a numeric matrix for the parameters, a list of two strings for the
//! variable types, an integer for the number of observations (if fitted),
//! a double for the log-likelihood (if fitted), and a double
//! for the number of parameters (can be non-integer in nonparametric
//! models).
//!
//! @param filename The name of the file to write.
inline void
Bicop::to_file(const std::string& filename) const
{
  tools_serialization::json_to_file(filename, to_json());
}

//! @brief Evaluates the copula density.
//!
//! @details The copula density is defined as joint density divided by marginal
//! densities, irrespective of variable types.
//!
//! When at least one variable is discrete, more than two columns are required
//! for `u`: the first \f$ n \times 2 \f$ block contains realizations of
//! \f$ (F_{X_1}(x_1), F_{X_2}(x_2)) \f$. The second \f$ n \times 2 \f$ block
//! contains realizations of \f$ (F_{X_1}(x_1^-), F_{X_2}(x_2^-)) \f$. The minus
//! indicates a left-sided limit of the cdf. For, e.g., an integer-valued
//! variable, it holds \f$ F_{X_1}(x_1^-) = F_{X_1}(x_1 - 1) \f$. For continuous
//! variables the left limit and the cdf itself coincide. Respective columns can
//! be omitted in the second block.
//!
//! @param u An \f$ n \times (2 + k) \f$ matrix of observations contained in
//!   \f$(0, 1) \f$, where \f$ k \f$ is the number of discrete variables.
//! @return A length n vector of copula densities evaluated at \c u.
inline Eigen::VectorXd
Bicop::pdf(const Eigen::MatrixXd& u) const
{
  check_data(u);
  return bicop_->pdf(prep_for_abstract(u));
}

//! @brief Evaluates the copula distribution.
//!
//! @details When at least one variable is discrete, more than two
//! columns are required
//! for `u`: the first \f$ n \times 2 \f$ block contains realizations of
//! \f$ (F_{X_1}(x_1), F_{X_2}(x_2)) \f$. The second \f$ n \times 2 \f$ block
//! contains realizations of \f$ (F_{X_1}(x_1^-), F_{X_2}(x_2^-)) \f$. The minus
//! indicates a left-sided limit of the cdf. For, e.g., an integer-valued
//! variable, it holds \f$ F_{X_1}(x_1^-) = F_{X_1}(x_1 - 1) \f$. For continuous
//! variables the left limit and the cdf itself coincide. Respective columns can
//! be omitted in the second block.
//!
//! @param u An \f$ n \times (2 + k) \f$ matrix of observations contained in
//!   \f$(0, 1) \f$, where \f$ k \f$ is the number of discrete variables.
//! @return A length n vector of copula probabilities evaluated at \c u.
inline Eigen::VectorXd
Bicop::cdf(const Eigen::MatrixXd& u) const
{
  check_data(u);
  Eigen::VectorXd p = bicop_->cdf(prep_for_abstract(u).leftCols(2),
                                  bicop_->get_parameters().transpose());
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
//! @details The first h-function is
//! \f$ h_1(u_1, u_2) = P(U_2 \le u_2 | U_1 = u_1) \f$.
//!
//! When at least one variable is discrete, more than two columns are required
//! for `u`: the first \f$ n \times 2 \f$ block contains realizations of
//! \f$ (F_{X_1}(x_1), F_{X_2}(x_2)) \f$. The second \f$ n \times 2 \f$ block
//! contains realizations of \f$ (F_{X_1}(x_1^-), F_{X_2}(x_2^-)) \f$. The minus
//! indicates a left-sided limit of the cdf. For, e.g., an integer-valued
//! variable, it holds \f$ F_{X_1}(x_1^-) = F_{X_1}(x_1 - 1) \f$. For continuous
//! variables the left limit and the cdf itself coincide. Respective columns can
//! be omitted in the second block.
//!
//! @param u An \f$ n \times (2 + k) \f$ matrix of observations contained in
//!   \f$(0, 1) \f$, where \f$ k \f$ is the number of discrete variables.
//! @return A length n vector of the first h-function evaluated at \c u.
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
//! @details The second h-function is
//! \f$ h_2(u_1, u_2) = P(U_1 \le u_1 | U_2 = u_2)  \f$.
//!
//! When at least one variable is discrete, more than two columns are required
//! for `u`: the first \f$ n \times 2 \f$ block contains realizations of
//! \f$ (F_{X_1}(x_1), F_{X_2}(x_2)) \f$. The second \f$ n \times 2 \f$ block
//! contains realizations of \f$ (F_{X_1}(x_1^-), F_{X_2}(x_2^-)) \f$. The minus
//! indicates a left-sided limit of the cdf. For, e.g., an integer-valued
//! variable, it holds \f$ F_{X_1}(x_1^-) = F_{X_1}(x_1 - 1) \f$. For continuous
//! variables the left limit and the cdf itself coincide. Respective columns can
//! be omitted in the second block.
//!
//! @param u An \f$ n \times (2 + k) \f$ matrix of observations contained in
//!   \f$(0, 1) \f$, where \f$ k \f$ is the number of discrete variables.
//! @return A length n vector of the second h-function evaluated at \c u.
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
//! @details The first h-function is
//! \f$ h_1(u_1, u_2) = P(U_2 \le u_2 | U_1 = u_1) \f$.
//! The inverse is calulated w.r.t. the second argument.
//!
//! When at least one variable is discrete, more than two columns are required
//! for `u`: the first \f$ n \times 2 \f$ block contains realizations of
//! \f$ (F_{X_1}(x_1), F_{X_2}(x_2)) \f$. The second \f$ n \times 2 \f$ block
//! contains realizations of \f$ (F_{X_1}(x_1^-), F_{X_2}(x_2^-)) \f$. The minus
//! indicates a left-sided limit of the cdf. For, e.g., an integer-valued
//! variable, it holds \f$ F_{X_1}(x_1^-) = F_{X_1}(x_1 - 1) \f$. For continuous
//! variables the left limit and the cdf itself coincide. Respective columns can
//! be omitted in the second block.
//!
//! @param u An \f$ n \times (2 + k) \f$ matrix of observations contained in
//!   \f$(0, 1) \f$, where \f$ k \f$ is the number of discrete variables.
//! @return A length n vector of the inverse of the first h-function evaluated
//! at \c u.
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
//! @details The second h-function is
//! \f$ h_2(u_1, u_2) = P(U_1 \le u_1 | U_2 = u_2)  \f$.
//! The inverse is calculated w.r.t. the first argument.
//!
//! When at least one variable is discrete, more than two columns are required
//! for `u`: the first \f$ n \times 2 \f$ block contains realizations of
//! \f$ (F_{X_1}(x_1), F_{X_2}(x_2)) \f$. The second \f$ n \times 2 \f$ block
//! contains realizations of \f$ (F_{X_1}(x_1^-), F_{X_2}(x_2^-)) \f$. The minus
//! indicates a left-sided limit of the cdf. For, e.g., an integer-valued
//! variable, it holds \f$ F_{X_1}(x_1^-) = F_{X_1}(x_1 - 1) \f$. For continuous
//! variables the left limit and the cdf itself coincide. Respective columns can
//! be omitted in the second block.
//!
//! @param u An \f$ n \times (2 + k) \f$ matrix of observations contained in
//!   \f$(0, 1) \f$, where \f$ k \f$ is the number of discrete variables.
//! @return A length n vector of the inverse of the second h-function evaluated
//! at \c u.
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

//! @name Stats methods with per-row parameters
//!
//! @details These overloads evaluate the copula at a *different* parameter set
//! per row of `u`, in a single call, without mutating the object's stored
//! parameters. The family, rotation, and variable types are taken from the
//! object; only the parameter values vary by row. They are available for
//! parametric families only (nonparametric families store an interpolation
//! grid rather than a per-observation parameter vector).
//!
//! @param u An \f$ n \times (2 + k) \f$ matrix of observations (see the
//!   single-argument overloads for the layout with discrete variables).
//! @param parameters An \f$ n \times p \f$ matrix of parameters, where `p` is
//!   the number of family parameters (`get_parameters().size()`) and row `i`
//!   holds the parameter set used for row `i` of `u`. Parameters are given in
//!   the family's natural (unrotated) parameterization, as for
//!   `get_parameters()`.
//! @param num_threads The number of threads to parallelize the evaluation over
//!   rows.
//! @{

//! @brief Evaluates the copula density with per-row parameters.
inline Eigen::VectorXd
Bicop::pdf(const Eigen::MatrixXd& u,
           const Eigen::MatrixXd& parameters,
           const size_t num_threads) const
{
  Eigen::MatrixXd par_t = format_parameters(u, parameters);
  return eval_in_batches(u,
                         par_t,
                         num_threads,
                         [this](const Eigen::MatrixXd& ub,
                                const Eigen::MatrixXd& pb) -> Eigen::VectorXd {
                           return bicop_->pdf(prep_for_abstract(ub), pb);
                         });
}

//! @brief Evaluates the copula distribution with per-row parameters.
inline Eigen::VectorXd
Bicop::cdf(const Eigen::MatrixXd& u,
           const Eigen::MatrixXd& parameters,
           const size_t num_threads) const
{
  Eigen::MatrixXd par_t = format_parameters(u, parameters);
  return eval_in_batches(
    u,
    par_t,
    num_threads,
    [this](const Eigen::MatrixXd& ub,
           const Eigen::MatrixXd& pb) -> Eigen::VectorXd {
      Eigen::VectorXd p = bicop_->cdf(prep_for_abstract(ub).leftCols(2), pb);
      switch (rotation_) {
        case 90:
          return ub.col(1) - p;
        case 180:
          return (p.array() - 1 + ub.leftCols(2).rowwise().sum().array())
            .matrix();
        case 270:
          return ub.col(0) - p;
        default:
          return p;
      }
    });
}

//! @brief Evaluates the first h-function with per-row parameters.
inline Eigen::VectorXd
Bicop::hfunc1(const Eigen::MatrixXd& u,
              const Eigen::MatrixXd& parameters,
              const size_t num_threads) const
{
  Eigen::MatrixXd par_t = format_parameters(u, parameters);
  return eval_in_batches(
    u,
    par_t,
    num_threads,
    [this](const Eigen::MatrixXd& ub,
           const Eigen::MatrixXd& pb) -> Eigen::VectorXd {
      Eigen::VectorXd h(ub.rows());
      switch (rotation_) {
        case 90:
          h = bicop_->hfunc2(prep_for_abstract(ub), pb);
          break;
        case 180:
          h = 1.0 - bicop_->hfunc1(prep_for_abstract(ub), pb).array();
          break;
        case 270:
          h = 1.0 - bicop_->hfunc2(prep_for_abstract(ub), pb).array();
          break;
        default:
          h = bicop_->hfunc1(prep_for_abstract(ub), pb);
          break;
      }
      tools_eigen::trim(h, 0.0, 1.0);
      return h;
    });
}

//! @brief Evaluates the second h-function with per-row parameters.
inline Eigen::VectorXd
Bicop::hfunc2(const Eigen::MatrixXd& u,
              const Eigen::MatrixXd& parameters,
              const size_t num_threads) const
{
  Eigen::MatrixXd par_t = format_parameters(u, parameters);
  return eval_in_batches(
    u,
    par_t,
    num_threads,
    [this](const Eigen::MatrixXd& ub,
           const Eigen::MatrixXd& pb) -> Eigen::VectorXd {
      Eigen::VectorXd h(ub.rows());
      switch (rotation_) {
        case 90:
          h = 1.0 - bicop_->hfunc1(prep_for_abstract(ub), pb).array();
          break;
        case 180:
          h = 1.0 - bicop_->hfunc2(prep_for_abstract(ub), pb).array();
          break;
        case 270:
          h = bicop_->hfunc1(prep_for_abstract(ub), pb);
          break;
        default:
          h = bicop_->hfunc2(prep_for_abstract(ub), pb);
          break;
      }
      tools_eigen::trim(h, 0.0, 1.0);
      return h;
    });
}

//! @brief Evaluates the inverse of the first h-function with per-row
//! parameters.
inline Eigen::VectorXd
Bicop::hinv1(const Eigen::MatrixXd& u,
             const Eigen::MatrixXd& parameters,
             const size_t num_threads) const
{
  Eigen::MatrixXd par_t = format_parameters(u, parameters);
  return eval_in_batches(
    u,
    par_t,
    num_threads,
    [this](const Eigen::MatrixXd& ub,
           const Eigen::MatrixXd& pb) -> Eigen::VectorXd {
      Eigen::VectorXd hi(ub.rows());
      switch (rotation_) {
        case 90:
          hi = bicop_->hinv2(prep_for_abstract(ub), pb);
          break;
        case 180:
          hi = 1.0 - bicop_->hinv1(prep_for_abstract(ub), pb).array();
          break;
        case 270:
          hi = 1.0 - bicop_->hinv2(prep_for_abstract(ub), pb).array();
          break;
        default:
          hi = bicop_->hinv1(prep_for_abstract(ub), pb);
          break;
      }
      tools_eigen::trim(hi, 0.0, 1.0);
      return hi;
    });
}

//! @brief Evaluates the inverse of the second h-function with per-row
//! parameters.
inline Eigen::VectorXd
Bicop::hinv2(const Eigen::MatrixXd& u,
             const Eigen::MatrixXd& parameters,
             const size_t num_threads) const
{
  Eigen::MatrixXd par_t = format_parameters(u, parameters);
  return eval_in_batches(
    u,
    par_t,
    num_threads,
    [this](const Eigen::MatrixXd& ub,
           const Eigen::MatrixXd& pb) -> Eigen::VectorXd {
      Eigen::VectorXd hi(ub.rows());
      switch (rotation_) {
        case 90:
          hi = 1.0 - bicop_->hinv1(prep_for_abstract(ub), pb).array();
          break;
        case 180:
          hi = 1.0 - bicop_->hinv2(prep_for_abstract(ub), pb).array();
          break;
        case 270:
          hi = bicop_->hinv1(prep_for_abstract(ub), pb);
          break;
        default:
          hi = bicop_->hinv2(prep_for_abstract(ub), pb);
          break;
      }
      tools_eigen::trim(hi, 0.0, 1.0);
      return hi;
    });
}

//! @brief Evaluates the log-likelihood contributions with per-row parameters
//! and returns their sum (NaNs are ignored).
inline double
Bicop::loglik(const Eigen::MatrixXd& u,
              const Eigen::MatrixXd& parameters,
              const size_t num_threads) const
{
  Eigen::VectorXd lpdf = pdf(u, parameters, num_threads).array().log();
  double ll = 0.0;
  for (Eigen::Index i = 0; i < lpdf.size(); ++i) {
    if (!(std::isnan)(lpdf(i))) {
      ll += lpdf(i);
    }
  }
  return ll;
}
//! @}

//! @name Derivatives of the density and h-functions
//!
//! @details These methods evaluate partial derivatives of the copula density
//! \f$ c(u_1, u_2; \theta) \f$, its logarithm, and the h-functions with
//! respect to the parameters and/or the arguments. The derivative is chosen
//! by the selector `deriv`, a concatenation of the components `"par1"`,
//! `"par2"`, ... (the k-th parameter, in the family's natural parameterization
//! as returned by `get_parameters()`), `"u1"`, and `"u2"`. First-order
//! methods take one component (`"par"` is short for `"par1"`); second-order
//! methods take two in any order (a single component means differentiating
//! twice, so `"par"` is short for `"par1par1"`).
//!
//! Rotations are handled internally via the chain rule, so derivatives are
//! always w.r.t. the arguments and (positive, natural) parameters of the
//! rotated copula. Closed-form expressions are used for the families in
//! `bicop_families::analytic_derivs`; other parametric families fall back to
//! central finite differences of the corresponding function. Derivatives
//! require continuous variable types; nonparametric families throw.
//!
//! The per-row-parameter overloads evaluate the derivative at a different
//! parameter set per row of `u` (see the corresponding `pdf()` overload for
//! the layout and validation rules).
//!
//! @param u An \f$ n \times 2 \f$ matrix of observations contained in
//!   \f$ (0, 1)^2 \f$.
//! @param deriv The derivative selector.
//! @return A length n vector of derivatives evaluated at `u`.
//! @{

//! @brief Evaluates a first derivative of the copula density.
//!
//! @details `deriv` is one of `"par1"`, `"par2"`, ... (short form `"par"`),
//! `"u1"`, or `"u2"`.
inline Eigen::VectorXd
Bicop::pdf_deriv(const Eigen::MatrixXd& u, const std::string& deriv) const
{
  check_data(u);
  check_deriv_preconditions();
  auto spec = map_pdf_deriv(tools_deriv::canonicalize(deriv, 1, deriv_npars()));
  return spec.sign * bicop_->pdf_deriv_raw(prep_for_abstract(u).leftCols(2),
                                           bicop_->get_parameters().transpose(),
                                           spec.deriv);
}

//! @brief Evaluates a second derivative of the copula density.
//!
//! @details `deriv` combines two first-order components, e.g. `"par1par1"`,
//! `"par1u1"`, `"u1u2"`; a single component means differentiating twice.
inline Eigen::VectorXd
Bicop::pdf_deriv2(const Eigen::MatrixXd& u, const std::string& deriv) const
{
  check_data(u);
  check_deriv_preconditions();
  auto spec = map_pdf_deriv(tools_deriv::canonicalize(deriv, 2, deriv_npars()));
  return spec.sign *
         bicop_->pdf_deriv2_raw(prep_for_abstract(u).leftCols(2),
                                bicop_->get_parameters().transpose(),
                                spec.deriv);
}

//! @brief Evaluates a first derivative of the first h-function
//! \f$ h_1(u_1, u_2) = P(U_2 \le u_2 | U_1 = u_1) \f$.
//!
//! @details `deriv` is one of `"par1"`, `"par2"`, ... (short form `"par"`),
//! `"u1"` (the conditioning argument), or `"u2"` (which equals the copula
//! density).
inline Eigen::VectorXd
Bicop::hfunc1_deriv(const Eigen::MatrixXd& u, const std::string& deriv) const
{
  check_data(u);
  check_deriv_preconditions();
  auto canonical = tools_deriv::canonicalize(deriv, 1, deriv_npars());
  if (canonical == "u2") {
    return pdf(u); // dh1/du2 = c, at every rotation
  }
  auto spec = map_hfunc_deriv(canonical, true);
  Eigen::MatrixXd u_r = prep_for_abstract(u).leftCols(2);
  Eigen::MatrixXd pars = bicop_->get_parameters().transpose();
  if (spec.swap_hfunc) {
    return spec.sign * bicop_->hfunc2_deriv_raw(u_r, pars, spec.deriv);
  }
  return spec.sign * bicop_->hfunc1_deriv_raw(u_r, pars, spec.deriv);
}

//! @brief Evaluates a second derivative of the first h-function.
//!
//! @details `deriv` combines two first-order components; any selector
//! containing `"u2"` reduces to a density derivative (e.g. `"par1u2"` equals
//! `pdf_deriv(u, "par1")`).
inline Eigen::VectorXd
Bicop::hfunc1_deriv2(const Eigen::MatrixXd& u, const std::string& deriv) const
{
  check_data(u);
  check_deriv_preconditions();
  auto canonical = tools_deriv::canonicalize(deriv, 2, deriv_npars());
  auto comps = tools_deriv::parse_components(canonical);
  if (comps[1] == -2) { // sorted, so "u2" can only be the last component
    // dh1/du2 = c: reduce to a first derivative of the density
    return pdf_deriv(u, tools_deriv::comp_to_string(comps[0]));
  }
  auto spec = map_hfunc_deriv(canonical, true);
  Eigen::MatrixXd u_r = prep_for_abstract(u).leftCols(2);
  Eigen::MatrixXd pars = bicop_->get_parameters().transpose();
  if (spec.swap_hfunc) {
    return spec.sign * bicop_->hfunc2_deriv2_raw(u_r, pars, spec.deriv);
  }
  return spec.sign * bicop_->hfunc1_deriv2_raw(u_r, pars, spec.deriv);
}

//! @brief Evaluates a first derivative of the second h-function
//! \f$ h_2(u_1, u_2) = P(U_1 \le u_1 | U_2 = u_2) \f$.
//!
//! @details `deriv` is one of `"par1"`, `"par2"`, ... (short form `"par"`),
//! `"u2"` (the conditioning argument), or `"u1"` (which equals the copula
//! density).
inline Eigen::VectorXd
Bicop::hfunc2_deriv(const Eigen::MatrixXd& u, const std::string& deriv) const
{
  check_data(u);
  check_deriv_preconditions();
  auto canonical = tools_deriv::canonicalize(deriv, 1, deriv_npars());
  if (canonical == "u1") {
    return pdf(u); // dh2/du1 = c, at every rotation
  }
  auto spec = map_hfunc_deriv(canonical, false);
  Eigen::MatrixXd u_r = prep_for_abstract(u).leftCols(2);
  Eigen::MatrixXd pars = bicop_->get_parameters().transpose();
  if (spec.swap_hfunc) {
    return spec.sign * bicop_->hfunc1_deriv_raw(u_r, pars, spec.deriv);
  }
  return spec.sign * bicop_->hfunc2_deriv_raw(u_r, pars, spec.deriv);
}

//! @brief Evaluates a second derivative of the second h-function.
//!
//! @details `deriv` combines two first-order components; any selector
//! containing `"u1"` reduces to a density derivative (e.g. `"par1u1"` equals
//! `pdf_deriv(u, "par1")`).
inline Eigen::VectorXd
Bicop::hfunc2_deriv2(const Eigen::MatrixXd& u, const std::string& deriv) const
{
  check_data(u);
  check_deriv_preconditions();
  auto canonical = tools_deriv::canonicalize(deriv, 2, deriv_npars());
  auto comps = tools_deriv::parse_components(canonical);
  if ((comps[0] == -1) || (comps[1] == -1)) {
    // dh2/du1 = c: reduce to a first derivative of the density
    int other = (comps[0] == -1) ? comps[1] : comps[0];
    return pdf_deriv(u, tools_deriv::comp_to_string(other));
  }
  auto spec = map_hfunc_deriv(canonical, false);
  Eigen::MatrixXd u_r = prep_for_abstract(u).leftCols(2);
  Eigen::MatrixXd pars = bicop_->get_parameters().transpose();
  if (spec.swap_hfunc) {
    return spec.sign * bicop_->hfunc1_deriv2_raw(u_r, pars, spec.deriv);
  }
  return spec.sign * bicop_->hfunc2_deriv2_raw(u_r, pars, spec.deriv);
}

//! @brief Evaluates a first derivative of the log-density
//! \f$ \partial \log c / \partial \cdot \f$.
//!
//! @details For parameter selectors, dedicated closed forms are used where
//! available (numerically stabler than `pdf_deriv() / pdf()` when the
//! density is small); argument selectors are composed by the quotient rule.
inline Eigen::VectorXd
Bicop::logpdf_deriv(const Eigen::MatrixXd& u, const std::string& deriv) const
{
  check_data(u);
  check_deriv_preconditions();
  auto canonical = tools_deriv::canonicalize(deriv, 1, deriv_npars());
  auto comps = tools_deriv::parse_components(canonical);
  if (comps[0] >= 0) {
    // parameter selectors are rotation-invariant
    return bicop_->logpdf_deriv_raw(prep_for_abstract(u).leftCols(2),
                                    bicop_->get_parameters().transpose(),
                                    canonical);
  }
  // pdf() is trimmed to [DBL_MIN, DBL_MAX], so the denominator is positive
  Eigen::ArrayXd c = pdf(u).array();
  return (pdf_deriv(u, canonical).array() / c).matrix();
}

//! @brief Evaluates a second derivative of the log-density.
//!
//! @details See `logpdf_deriv()`; selectors involving the arguments are
//! composed from density derivatives by the quotient rule.
inline Eigen::VectorXd
Bicop::logpdf_deriv2(const Eigen::MatrixXd& u, const std::string& deriv) const
{
  check_data(u);
  check_deriv_preconditions();
  auto canonical = tools_deriv::canonicalize(deriv, 2, deriv_npars());
  auto comps = tools_deriv::parse_components(canonical);
  if (comps[1] >= 0) { // sorted, so both components are parameters
    return bicop_->logpdf_deriv2_raw(prep_for_abstract(u).leftCols(2),
                                     bicop_->get_parameters().transpose(),
                                     canonical);
  }
  // pdf() is trimmed to [DBL_MIN, DBL_MAX], so the denominator is positive
  Eigen::ArrayXd c = pdf(u).array();
  Eigen::ArrayXd c_xy = pdf_deriv2(u, canonical).array();
  Eigen::ArrayXd c_x =
    pdf_deriv(u, tools_deriv::comp_to_string(comps[0])).array();
  Eigen::ArrayXd c_y =
    (comps[0] == comps[1])
      ? c_x
      : pdf_deriv(u, tools_deriv::comp_to_string(comps[1])).array();
  return (c_xy / c - (c_x / c) * (c_y / c)).matrix();
}

//! @brief Evaluates a first derivative of the copula density with per-row
//! parameters.
inline Eigen::VectorXd
Bicop::pdf_deriv(const Eigen::MatrixXd& u,
                 const std::string& deriv,
                 const Eigen::MatrixXd& parameters,
                 const size_t num_threads) const
{
  check_deriv_preconditions();
  auto spec = map_pdf_deriv(tools_deriv::canonicalize(deriv, 1, deriv_npars()));
  Eigen::MatrixXd par_t = format_parameters(u, parameters);
  return eval_in_batches(
    u,
    par_t,
    num_threads,
    [this, spec](const Eigen::MatrixXd& ub,
                 const Eigen::MatrixXd& pb) -> Eigen::VectorXd {
      return spec.sign * bicop_->pdf_deriv_raw(
                           prep_for_abstract(ub).leftCols(2), pb, spec.deriv);
    });
}

//! @brief Evaluates a second derivative of the copula density with per-row
//! parameters.
inline Eigen::VectorXd
Bicop::pdf_deriv2(const Eigen::MatrixXd& u,
                  const std::string& deriv,
                  const Eigen::MatrixXd& parameters,
                  const size_t num_threads) const
{
  check_deriv_preconditions();
  auto spec = map_pdf_deriv(tools_deriv::canonicalize(deriv, 2, deriv_npars()));
  Eigen::MatrixXd par_t = format_parameters(u, parameters);
  return eval_in_batches(
    u,
    par_t,
    num_threads,
    [this, spec](const Eigen::MatrixXd& ub,
                 const Eigen::MatrixXd& pb) -> Eigen::VectorXd {
      return spec.sign * bicop_->pdf_deriv2_raw(
                           prep_for_abstract(ub).leftCols(2), pb, spec.deriv);
    });
}

//! @brief Evaluates a first derivative of the first h-function with per-row
//! parameters.
inline Eigen::VectorXd
Bicop::hfunc1_deriv(const Eigen::MatrixXd& u,
                    const std::string& deriv,
                    const Eigen::MatrixXd& parameters,
                    const size_t num_threads) const
{
  check_deriv_preconditions();
  auto canonical = tools_deriv::canonicalize(deriv, 1, deriv_npars());
  if (canonical == "u2") {
    return pdf(u, parameters, num_threads);
  }
  auto spec = map_hfunc_deriv(canonical, true);
  Eigen::MatrixXd par_t = format_parameters(u, parameters);
  return eval_in_batches(
    u,
    par_t,
    num_threads,
    [this, spec](const Eigen::MatrixXd& ub,
                 const Eigen::MatrixXd& pb) -> Eigen::VectorXd {
      Eigen::MatrixXd u_r = prep_for_abstract(ub).leftCols(2);
      if (spec.swap_hfunc) {
        return spec.sign * bicop_->hfunc2_deriv_raw(u_r, pb, spec.deriv);
      }
      return spec.sign * bicop_->hfunc1_deriv_raw(u_r, pb, spec.deriv);
    });
}

//! @brief Evaluates a second derivative of the first h-function with per-row
//! parameters.
inline Eigen::VectorXd
Bicop::hfunc1_deriv2(const Eigen::MatrixXd& u,
                     const std::string& deriv,
                     const Eigen::MatrixXd& parameters,
                     const size_t num_threads) const
{
  check_deriv_preconditions();
  auto canonical = tools_deriv::canonicalize(deriv, 2, deriv_npars());
  auto comps = tools_deriv::parse_components(canonical);
  if (comps[1] == -2) {
    return pdf_deriv(
      u, tools_deriv::comp_to_string(comps[0]), parameters, num_threads);
  }
  auto spec = map_hfunc_deriv(canonical, true);
  Eigen::MatrixXd par_t = format_parameters(u, parameters);
  return eval_in_batches(
    u,
    par_t,
    num_threads,
    [this, spec](const Eigen::MatrixXd& ub,
                 const Eigen::MatrixXd& pb) -> Eigen::VectorXd {
      Eigen::MatrixXd u_r = prep_for_abstract(ub).leftCols(2);
      if (spec.swap_hfunc) {
        return spec.sign * bicop_->hfunc2_deriv2_raw(u_r, pb, spec.deriv);
      }
      return spec.sign * bicop_->hfunc1_deriv2_raw(u_r, pb, spec.deriv);
    });
}

//! @brief Evaluates a first derivative of the second h-function with per-row
//! parameters.
inline Eigen::VectorXd
Bicop::hfunc2_deriv(const Eigen::MatrixXd& u,
                    const std::string& deriv,
                    const Eigen::MatrixXd& parameters,
                    const size_t num_threads) const
{
  check_deriv_preconditions();
  auto canonical = tools_deriv::canonicalize(deriv, 1, deriv_npars());
  if (canonical == "u1") {
    return pdf(u, parameters, num_threads);
  }
  auto spec = map_hfunc_deriv(canonical, false);
  Eigen::MatrixXd par_t = format_parameters(u, parameters);
  return eval_in_batches(
    u,
    par_t,
    num_threads,
    [this, spec](const Eigen::MatrixXd& ub,
                 const Eigen::MatrixXd& pb) -> Eigen::VectorXd {
      Eigen::MatrixXd u_r = prep_for_abstract(ub).leftCols(2);
      if (spec.swap_hfunc) {
        return spec.sign * bicop_->hfunc1_deriv_raw(u_r, pb, spec.deriv);
      }
      return spec.sign * bicop_->hfunc2_deriv_raw(u_r, pb, spec.deriv);
    });
}

//! @brief Evaluates a second derivative of the second h-function with
//! per-row parameters.
inline Eigen::VectorXd
Bicop::hfunc2_deriv2(const Eigen::MatrixXd& u,
                     const std::string& deriv,
                     const Eigen::MatrixXd& parameters,
                     const size_t num_threads) const
{
  check_deriv_preconditions();
  auto canonical = tools_deriv::canonicalize(deriv, 2, deriv_npars());
  auto comps = tools_deriv::parse_components(canonical);
  if ((comps[0] == -1) || (comps[1] == -1)) {
    int other = (comps[0] == -1) ? comps[1] : comps[0];
    return pdf_deriv(
      u, tools_deriv::comp_to_string(other), parameters, num_threads);
  }
  auto spec = map_hfunc_deriv(canonical, false);
  Eigen::MatrixXd par_t = format_parameters(u, parameters);
  return eval_in_batches(
    u,
    par_t,
    num_threads,
    [this, spec](const Eigen::MatrixXd& ub,
                 const Eigen::MatrixXd& pb) -> Eigen::VectorXd {
      Eigen::MatrixXd u_r = prep_for_abstract(ub).leftCols(2);
      if (spec.swap_hfunc) {
        return spec.sign * bicop_->hfunc1_deriv2_raw(u_r, pb, spec.deriv);
      }
      return spec.sign * bicop_->hfunc2_deriv2_raw(u_r, pb, spec.deriv);
    });
}

//! @brief Evaluates a first derivative of the log-density with per-row
//! parameters.
inline Eigen::VectorXd
Bicop::logpdf_deriv(const Eigen::MatrixXd& u,
                    const std::string& deriv,
                    const Eigen::MatrixXd& parameters,
                    const size_t num_threads) const
{
  check_deriv_preconditions();
  auto canonical = tools_deriv::canonicalize(deriv, 1, deriv_npars());
  auto comps = tools_deriv::parse_components(canonical);
  if (comps[0] >= 0) {
    Eigen::MatrixXd par_t = format_parameters(u, parameters);
    return eval_in_batches(
      u,
      par_t,
      num_threads,
      [this, canonical](const Eigen::MatrixXd& ub,
                        const Eigen::MatrixXd& pb) -> Eigen::VectorXd {
        return bicop_->logpdf_deriv_raw(
          prep_for_abstract(ub).leftCols(2), pb, canonical);
      });
  }
  Eigen::ArrayXd c = pdf(u, parameters, num_threads).array();
  return (pdf_deriv(u, canonical, parameters, num_threads).array() / c)
    .matrix();
}

//! @brief Evaluates a second derivative of the log-density with per-row
//! parameters.
inline Eigen::VectorXd
Bicop::logpdf_deriv2(const Eigen::MatrixXd& u,
                     const std::string& deriv,
                     const Eigen::MatrixXd& parameters,
                     const size_t num_threads) const
{
  check_deriv_preconditions();
  auto canonical = tools_deriv::canonicalize(deriv, 2, deriv_npars());
  auto comps = tools_deriv::parse_components(canonical);
  if (comps[1] >= 0) {
    Eigen::MatrixXd par_t = format_parameters(u, parameters);
    return eval_in_batches(
      u,
      par_t,
      num_threads,
      [this, canonical](const Eigen::MatrixXd& ub,
                        const Eigen::MatrixXd& pb) -> Eigen::VectorXd {
        return bicop_->logpdf_deriv2_raw(
          prep_for_abstract(ub).leftCols(2), pb, canonical);
      });
  }
  Eigen::ArrayXd c = pdf(u, parameters, num_threads).array();
  Eigen::ArrayXd c_xy =
    pdf_deriv2(u, canonical, parameters, num_threads).array();
  Eigen::ArrayXd c_x =
    pdf_deriv(u, tools_deriv::comp_to_string(comps[0]), parameters, num_threads)
      .array();
  Eigen::ArrayXd c_y =
    (comps[0] == comps[1])
      ? c_x
      : pdf_deriv(
          u, tools_deriv::comp_to_string(comps[1]), parameters, num_threads)
          .array();
  return (c_xy / c - (c_x / c) * (c_y / c)).matrix();
}
//! @}

//! @brief Assembles an \f$ n \times p \f$ score matrix from a per-parameter
//! column evaluator (`col(k)` returns column `k`). Shared by the fixed- and
//! per-row-parameter `scores()` overloads so each keeps its own optimal
//! `logpdf_deriv()` path (broadcast vs. per-row) while the loop lives once.
inline Eigen::MatrixXd
assemble_scores(Eigen::Index n,
                Eigen::Index p,
                const std::function<Eigen::VectorXd(Eigen::Index)>& col)
{
  Eigen::MatrixXd s(n, p);
  for (Eigen::Index k = 0; k < p; ++k) {
    s.col(k) = col(k);
  }
  return s;
}

//! @brief Assembles the averaged, symmetric \f$ p \times p \f$ Hessian from a
//! per-\f$ (a, b) \f$ second-derivative column evaluator (upper triangle only).
inline Eigen::MatrixXd
assemble_hessian(
  Eigen::Index p,
  const std::function<Eigen::VectorXd(Eigen::Index, Eigen::Index)>& col)
{
  Eigen::MatrixXd h(p, p);
  for (Eigen::Index a = 0; a < p; ++a) {
    for (Eigen::Index b = a; b < p; ++b) {
      h(a, b) = col(a, b).mean();
      h(b, a) = h(a, b);
    }
  }
  return h;
}

//! @brief Assembles the per-observation, symmetric \f$ p \times p \f$ Hessians
//! (one per row of `u`) from the same per-\f$ (a, b) \f$ column evaluator.
inline std::vector<Eigen::MatrixXd>
assemble_hessian_full(
  Eigen::Index n,
  Eigen::Index p,
  const std::function<Eigen::VectorXd(Eigen::Index, Eigen::Index)>& col)
{
  std::vector<Eigen::MatrixXd> hess(static_cast<size_t>(n),
                                    Eigen::MatrixXd(p, p));
  for (Eigen::Index a = 0; a < p; ++a) {
    for (Eigen::Index b = a; b < p; ++b) {
      Eigen::VectorXd d = col(a, b);
      for (Eigen::Index i = 0; i < n; ++i) {
        hess[static_cast<size_t>(i)](a, b) = d(i);
        hess[static_cast<size_t>(i)](b, a) = d(i);
      }
    }
  }
  return hess;
}

//! @name Scores, gradient, and Hessian of the log-likelihood
//!
//! @details These methods aggregate the log-density parameter derivatives into
//! the score (the gradient of each observation's log-density contribution), its
//! observation-average (`gradient`), and the Hessian. They are thin wrappers
//! around `logpdf_deriv()` / `logpdf_deriv2()`: with `p = get_parameters()`
//! parameters, the score of observation `i` w.r.t. parameter `k` is
//! \f$ \partial \log c(u_i; \theta) / \partial \theta_k \f$, and the (per-obs
//! and averaged) Hessians collect the second log-density derivatives. Like the
//! derivatives they build on, they require parametric families and continuous
//! variable types; nonparametric or discrete models throw.
//!
//! The per-row-parameter overloads evaluate at a different parameter set per
//! row of `u` (see the corresponding `pdf()` overload for the layout and
//! validation rules).
//!
//! @param u An \f$ n \times 2 \f$ matrix of observations contained in
//!   \f$ (0, 1)^2 \f$.
//! @{

//! @brief Evaluates the per-observation scores.
//!
//! @return An \f$ n \times p \f$ matrix whose column `k` is
//! \f$ \partial \log c / \partial \theta_{k+1} \f$.
inline Eigen::MatrixXd
Bicop::scores(const Eigen::MatrixXd& u) const
{
  check_data(u);
  check_deriv_preconditions();
  return assemble_scores(
    u.rows(), static_cast<Eigen::Index>(deriv_npars()), [&](Eigen::Index k) {
      return logpdf_deriv(u, "par" + std::to_string(k + 1));
    });
}

//! @brief Evaluates the gradient of the average log-likelihood.
//!
//! @return The observation-average of `scores()`, a vector of length `p`.
inline Eigen::VectorXd
Bicop::gradient(const Eigen::MatrixXd& u) const
{
  return scores(u).colwise().mean().transpose();
}

//! @brief Evaluates the Hessian of the average log-likelihood.
//!
//! @return A symmetric \f$ p \times p \f$ matrix whose entry \f$ (a, b) \f$ is
//! the observation-average of
//! \f$ \partial^2 \log c / \partial \theta_{a+1} \partial \theta_{b+1} \f$.
inline Eigen::MatrixXd
Bicop::hessian(const Eigen::MatrixXd& u) const
{
  check_data(u);
  check_deriv_preconditions();
  return assemble_hessian(
    static_cast<Eigen::Index>(deriv_npars()),
    [&](Eigen::Index a, Eigen::Index b) {
      return logpdf_deriv2(
        u, "par" + std::to_string(a + 1) + "par" + std::to_string(b + 1));
    });
}

//! @brief Evaluates the per-observation Hessians.
//!
//! @return A vector of `n` symmetric \f$ p \times p \f$ matrices; entry `i`'s
//! \f$ (a, b) \f$ element is
//! \f$ \partial^2 \log c(u_i; \theta) / \partial \theta_{a+1} \partial
//! \theta_{b+1} \f$.
inline std::vector<Eigen::MatrixXd>
Bicop::hessian_full(const Eigen::MatrixXd& u) const
{
  check_data(u);
  check_deriv_preconditions();
  return assemble_hessian_full(
    u.rows(),
    static_cast<Eigen::Index>(deriv_npars()),
    [&](Eigen::Index a, Eigen::Index b) {
      return logpdf_deriv2(
        u, "par" + std::to_string(a + 1) + "par" + std::to_string(b + 1));
    });
}

//! @brief Computes the covariance matrix of the scores.
//!
//! @return The mean-centered, divided-by-`n` covariance of `scores()`, a
//! \f$ p \times p \f$ matrix.
inline Eigen::MatrixXd
Bicop::scores_cov(const Eigen::MatrixXd& u) const
{
  Eigen::MatrixXd s = scores(u);
  // materialize the centered scores; a lazy expression would be evaluated
  // twice by the product below
  Eigen::MatrixXd sc = s.rowwise() - s.colwise().mean();
  return (sc.adjoint() * sc) / static_cast<double>(s.rows());
}

//! @brief Evaluates the scores, bundled in a `ScoresResult`.
//!
//! @details Provided for parity with `Vinecop::scores_full()`; a single pair
//! copula has no cascade caches, so the result only carries the score matrix.
inline Bicop::ScoresResult
Bicop::scores_full(const Eigen::MatrixXd& u) const
{
  ScoresResult result;
  result.scores = scores(u);
  return result;
}

//! @brief Evaluates the per-observation scores with per-row parameters.
inline Eigen::MatrixXd
Bicop::scores(const Eigen::MatrixXd& u,
              const Eigen::MatrixXd& parameters,
              const size_t num_threads) const
{
  check_deriv_preconditions();
  return assemble_scores(
    u.rows(), static_cast<Eigen::Index>(deriv_npars()), [&](Eigen::Index k) {
      return logpdf_deriv(
        u, "par" + std::to_string(k + 1), parameters, num_threads);
    });
}

//! @brief Evaluates the gradient of the average log-likelihood with per-row
//! parameters.
inline Eigen::VectorXd
Bicop::gradient(const Eigen::MatrixXd& u,
                const Eigen::MatrixXd& parameters,
                const size_t num_threads) const
{
  return scores(u, parameters, num_threads).colwise().mean().transpose();
}

//! @brief Evaluates the Hessian of the average log-likelihood with per-row
//! parameters.
inline Eigen::MatrixXd
Bicop::hessian(const Eigen::MatrixXd& u,
               const Eigen::MatrixXd& parameters,
               const size_t num_threads) const
{
  check_deriv_preconditions();
  return assemble_hessian(static_cast<Eigen::Index>(deriv_npars()),
                          [&](Eigen::Index a, Eigen::Index b) {
                            return logpdf_deriv2(u,
                                                 "par" + std::to_string(a + 1) +
                                                   "par" +
                                                   std::to_string(b + 1),
                                                 parameters,
                                                 num_threads);
                          });
}

//! @brief Evaluates the per-observation Hessians with per-row parameters.
inline std::vector<Eigen::MatrixXd>
Bicop::hessian_full(const Eigen::MatrixXd& u,
                    const Eigen::MatrixXd& parameters,
                    const size_t num_threads) const
{
  check_deriv_preconditions();
  return assemble_hessian_full(u.rows(),
                               static_cast<Eigen::Index>(deriv_npars()),
                               [&](Eigen::Index a, Eigen::Index b) {
                                 return logpdf_deriv2(
                                   u,
                                   "par" + std::to_string(a + 1) + "par" +
                                     std::to_string(b + 1),
                                   parameters,
                                   num_threads);
                               });
}

//! @brief Computes the covariance matrix of the scores with per-row parameters.
inline Eigen::MatrixXd
Bicop::scores_cov(const Eigen::MatrixXd& u,
                  const Eigen::MatrixXd& parameters,
                  const size_t num_threads) const
{
  Eigen::MatrixXd s = scores(u, parameters, num_threads);
  Eigen::MatrixXd sc = s.rowwise() - s.colwise().mean();
  return (sc.adjoint() * sc) / static_cast<double>(s.rows());
}

//! @brief Evaluates the scores with per-row parameters, bundled in a
//! `ScoresResult`.
inline Bicop::ScoresResult
Bicop::scores_full(const Eigen::MatrixXd& u,
                   const Eigen::MatrixXd& parameters,
                   const size_t num_threads) const
{
  ScoresResult result;
  result.scores = scores(u, parameters, num_threads);
  return result;
}
//! @}

//! checks that derivatives are available for the model (parametric family,
//! continuous variable types).
inline void
Bicop::check_deriv_preconditions() const
{
  if (!tools_stl::is_member(get_family(), bicop_families::parametric)) {
    throw std::runtime_error("derivatives are not implemented for the " +
                             get_family_name() + " copula");
  }
  if (var_types_ != std::vector<std::string>{ "c", "c" }) {
    throw std::runtime_error(
      "derivatives are only available for continuous variable types");
  }
}

//! the number of parameters used to validate derivative selectors.
inline size_t
Bicop::deriv_npars() const
{
  return static_cast<size_t>(bicop_->get_parameters().size());
}

//! @brief Resolves the rotation for a density-derivative selector.
//!
//! @details The rotated density is the unrotated one at transformed
//! arguments (90: \f$ (u_2, 1 - u_1) \f$, 180: \f$ (1 - u_1, 1 - u_2) \f$,
//! 270: \f$ (1 - u_2, u_1) \f$), so parameter components pass through
//! unchanged while argument components map to the leaf's argument slot and
//! pick up the chain-rule sign of the transform.
inline Bicop::DerivSpec
Bicop::map_pdf_deriv(const std::string& canonical) const
{
  auto comps = tools_deriv::parse_components(canonical);
  double sign = 1.0;
  for (auto& comp : comps) {
    if (comp >= 0) {
      continue;
    }
    switch (rotation_) {
      default:
        break;
      case 90:
        sign *= (comp == -1) ? -1.0 : 1.0;
        comp = (comp == -1) ? -2 : -1;
        break;
      case 180:
        sign *= -1.0;
        break;
      case 270:
        sign *= (comp == -2) ? -1.0 : 1.0;
        comp = (comp == -1) ? -2 : -1;
        break;
    }
  }
  return { tools_deriv::components_to_string(comps), sign, false };
}

//! @brief Resolves the rotation for an h-function-derivative selector.
//!
//! @details Under 90/270 rotations the rotated h-function is built from the
//! *other* h-function's leaf (`swap_hfunc`), mirroring `hfunc1()`/`hfunc2()`.
//! The sign is the product of the `1 - h` output flip (180/270 for the first
//! h-function, 90/180 for the second) and one argument chain-rule factor per
//! conditioning-argument component; parameter components pass through
//! unchanged. Selectors containing the conditioned argument must be reduced
//! to density derivatives before calling this.
inline Bicop::DerivSpec
Bicop::map_hfunc_deriv(const std::string& canonical, bool first_hfunc) const
{
  auto comps = tools_deriv::parse_components(canonical);
  bool swap = (rotation_ == 90) || (rotation_ == 270);
  double sign, chain;
  if (first_hfunc) {
    sign = ((rotation_ == 180) || (rotation_ == 270)) ? -1.0 : 1.0;
    chain = ((rotation_ == 90) || (rotation_ == 180)) ? -1.0 : 1.0;
  } else {
    sign = ((rotation_ == 90) || (rotation_ == 180)) ? -1.0 : 1.0;
    chain = ((rotation_ == 180) || (rotation_ == 270)) ? -1.0 : 1.0;
  }
  for (auto& comp : comps) {
    if (comp >= 0) {
      continue;
    }
    sign *= chain;
    if (swap) {
      comp = (comp == -1) ? -2 : -1;
    }
  }
  return { tools_deriv::components_to_string(comps), sign, swap };
}

//! @brief Validates per-row parameters (the internal leaves use the same
//! `n x p` layout as the public API, one parameter set per row).
inline Eigen::MatrixXd
Bicop::format_parameters(const Eigen::MatrixXd& u,
                         const Eigen::MatrixXd& parameters) const
{
  check_data(u);
  if (!tools_stl::is_member(get_family(), bicop_families::parametric)) {
    throw std::runtime_error(
      "per-row parameters are only supported for parametric families; "
      "nonparametric families store an interpolation grid rather than a "
      "per-observation parameter vector.");
  }
  const Eigen::Index n = u.rows();
  const Eigen::Index p = get_parameters().rows();
  if (parameters.rows() != n) {
    throw std::runtime_error("parameters must have one row per row of u "
                             "(parameters.rows() must equal u.rows()).");
  }
  if (parameters.cols() != p) {
    std::stringstream msg;
    msg << "parameters has wrong number of columns; expected " << p
        << " (the number of family parameters), actual " << parameters.cols()
        << ".";
    throw std::runtime_error(msg.str());
  }
  if (!parameters.allFinite()) {
    throw std::runtime_error("parameters must not contain NaN or Inf.");
  }
  if (p > 0 && n > 0) {
    Eigen::MatrixXd lb = get_parameters_lower_bounds();
    Eigen::MatrixXd ub = get_parameters_upper_bounds();
    for (Eigen::Index k = 0; k < p; ++k) {
      if (parameters.col(k).minCoeff() < lb(k) ||
          parameters.col(k).maxCoeff() > ub(k)) {
        std::stringstream msg;
        msg << "parameter " << k << " is out of bounds [" << lb(k) << ", "
            << ub(k) << "].";
        throw std::runtime_error(msg.str());
      }
    }
  }
  return parameters;
}

//! @brief Evaluates `f` over row-batches of `u`/`parameters`, possibly in
//! parallel, and assembles the results.
inline Eigen::VectorXd
Bicop::eval_in_batches(
  const Eigen::MatrixXd& u,
  const Eigen::MatrixXd& parameters,
  const size_t num_threads,
  const std::function<Eigen::VectorXd(const Eigen::MatrixXd&,
                                      const Eigen::MatrixXd&)>& f) const
{
  const size_t n = static_cast<size_t>(u.rows());
  Eigen::VectorXd out(n);
  if (n == 0) {
    return out;
  }
  auto do_batch = [&](const tools_batch::Batch& b) {
    out.segment(b.begin, b.size) =
      f(u.middleRows(b.begin, b.size), parameters.middleRows(b.begin, b.size));
  };
  if (num_threads <= 1) {
    do_batch(tools_batch::Batch{ 0, n });
  } else {
    tools_thread::ThreadPool pool(num_threads);
    pool.map(do_batch, tools_batch::create_batches(n, num_threads));
    pool.join();
  }
  return out;
}

//! @brief Simulates from a bivariate copula.
//!
//! @details If `qrng = TRUE`, generalized Halton sequences are used.
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
//! @details The log-likelihood is defined as
//! \f[ \mathrm{loglik} = \sum_{i = 1}^n \log c(U_{1, i}, U_{2, i}), \f]
//! where \f$ c \f$ is the copula density, see `Bicop::pdf()`.
//!
//! When at least one variable is discrete, more than two columns are required
//! for `u`: the first \f$ n \times 2 \f$ block contains realizations of
//! \f$ (F_{X_1}(x_1), F_{X_2}(x_2)) \f$. The second \f$ n \times 2 \f$ block
//! contains realizations of \f$ (F_{X_1}(x_1^-), F_{X_2}(x_2^-)) \f$. The minus
//! indicates a left-sided limit of the cdf. For, e.g., an integer-valued
//! variable, it holds \f$ F_{X_1}(x_1^-) = F_{X_1}(x_1 - 1) \f$. For continuous
//! variables the left limit and the cdf itself coincide. Respective columns can
//! be omitted in the second block.
//!
//! @param u An \f$ n \times (2 + k) \f$ matrix of observations contained in
//!   \f$(0, 1) \f$, where \f$ k \f$ is the number of discrete variables.
//! @return The log-likelihood evaluated at \c u.
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
//! @details The AIC is defined as
//! \f[ \mathrm{AIC} = -2\, \mathrm{loglik} + 2 p, \f]
//! where \f$ \mathrm{loglik} \f$ is the log-liklihood (see `Bicop::loglik()`)
//! and \f$ p \f$ is the (effective) number of parameters of the model.
//! The AIC is a consistent model selection criterion even
//! for nonparametric models.
//!
//! @param u An \f$ n \times (2 + k) \f$ matrix of observations contained in
//!   \f$(0, 1) \f$, where \f$ k \f$ is the number of discrete variables.
//! @return The AIC evaluated at \c u.
inline double
Bicop::aic(const Eigen::MatrixXd& u) const
{
  return -2 * loglik(u) + 2 * get_npars();
}

//! @brief Evaluates the Bayesian information criterion (BIC).
//!
//! The BIC is defined as
//! \f[ \mathrm{BIC} = -2\, \mathrm{loglik} +  \log(n) p, \f]
//! where \f$ \mathrm{loglik} \f$ is the log-liklihood (see `Bicop::loglik()`)
//! and \f$ p \f$ is the (effective) number of parameters of the model.
//! The BIC is a consistent model selection criterion
//! for parametric models.
//!
//! @param u An \f$ n \times (2 + k) \f$ matrix of observations contained in
//!   \f$(0, 1) \f$, where \f$ k \f$ is the number of discrete variables.
//! @return The BIC evaluated at \c u.
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
//! @details The mBIC is defined as
//! \f[ \mathrm{BIC} = -2\, \mathrm{loglik} +  p \log(n) - 2 (I \log(\psi_0) + (1 - I) \log(1 - \psi_0), \f]
//! where \f$ \mathrm{loglik} \f$ is the \log-liklihood
//! (see `Bicop::loglik()`), \f$ p \f$ is the (effective) number of parameters of the
//! model, and \f$ \psi_0 \f$ is the prior probability of having a
//! non-independence copula and \f$ I \f$ is an indicator for the family being
//! non-independence.
//!
//! @param u An \f$ n \times (2 + k) \f$ matrix of observations contained in
//!   \f$(0, 1) \f$, where \f$ k \f$ is the number of discrete variables.
//! @param psi0 Prior probability of a non-independence copula.
//! @return The mBIC evaluated at \c u.
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

//! @brief The number of parameters of the copula model.
//!
//! @details Returns the actual number of parameters for parameteric families.
//! For nonparametric families, there is a conceptually similar definition in
//! the sense that it can be used in the calculation of fit statistics.
inline double
Bicop::get_npars() const
{
  return bicop_->get_npars();
}

//! @brief Converts a Kendall's \f$ \tau \f$ into copula parameters
//! for one-parameter families.
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

//! @brief Converts the copula parameters to the tail dependence coefficients.
//!
//! @details The result is a \f$ 2 \times 2 \f$ matrix \f$ M \f$ collecting the
//! tail dependence coefficients in the four corners of the unit square:
//! \f$ M(i, j) \f$ is the coefficient as \f$ U_1 \to i \f$ and
//! \f$ U_2 \to j \f$, with \f$ i, j \in \{0, 1\} \f$ (0 = lower, 1 = upper).
//! Thus \f$ M(0, 0) \f$ is the classical lower and \f$ M(1, 1) \f$ the
//! classical upper tail dependence coefficient.
//!
//! For the unrotated family, the lower \f$ \lambda_L = M(0, 0) \f$ and upper
//! \f$ \lambda_U = M(1, 1) \f$ coefficients are:
//!
//! | Family | Lower \f$ \lambda_L \f$ | Upper \f$ \lambda_U \f$ |
//! | --- | --- | --- |
//! | Independence, Gaussian, Frank | \f$ 0 \f$ | \f$ 0 \f$ |
//! | Student t | \f$ \lambda_t \f$ (see below) | \f$ \lambda_t \f$ |
//! | Clayton | \f$ 2^{-1/\theta} \f$ | \f$ 0 \f$ |
//! | Gumbel | \f$ 0 \f$ | \f$ 2 - 2^{1/\theta} \f$ |
//! | Joe | \f$ 0 \f$ | \f$ 2 - 2^{1/\theta} \f$ |
//! | BB1 | \f$ 2^{-1/(\theta\delta)} \f$ | \f$ 2 - 2^{1/\delta} \f$ |
//! | BB6 | \f$ 0 \f$ | \f$ 2 - 2^{1/(\theta\delta)} \f$ |
//! | BB7 | \f$ 2^{-1/\delta} \f$ | \f$ 2 - 2^{1/\theta} \f$ |
//! | BB8 | \f$ 0 \f$ | \f$ 2 - 2^{1/\theta} \f$ if \f$ \delta=1 \f$, else 0 |
//! | Tawn (extreme-value) | \f$ 0 \f$ | \f$ 2 (1 - A(1/2)) \f$ |
//! | TLL (nonparametric) | \c NaN | \c NaN |
//!
//! Here \f$ \theta \f$ (and \f$ \delta \f$ for two-parameter families) denote
//! the copula parameters and \f$ A \f$ the Pickands dependence function. For
//! the Student t copula, with correlation \f$ \rho \f$, degrees of freedom
//! \f$ \nu \f$, and \f$ t_{\nu+1} \f$ the Student t cdf,
//! \f$ \lambda_t = 2\, t_{\nu+1}(-\sqrt{(\nu+1)(1-\rho)/(1+\rho)}) \f$; it
//! additionally has equal dependence in the two off-diagonal (discordant)
//! corners, obtained by replacing \f$ \rho \f$ with \f$ -\rho \f$. All other
//! parametric families have zero off-diagonal coefficients. Rotations permute
//! the corners: 180 degrees swaps lower/upper, while 90/270 degrees move
//! dependence to the off-diagonal corners.
//!
//! @param parameters The parameters (must be a valid parametrization of
//!     the current family).
inline Eigen::MatrixXd
Bicop::parameters_to_taildep(const Eigen::MatrixXd& parameters) const
{
  Eigen::MatrixXd m = bicop_->parameters_to_taildep(parameters);
  // rotate the 2x2 matrix like a grid, counter-clockwise, to match the
  // (counter-clockwise) rotation applied to the data in `rotate_data`.
  switch (rotation_) {
    case 90:
      return m.transpose().colwise().reverse();
    case 180:
      return m.reverse();
    case 270:
      return m.transpose().rowwise().reverse();
    default:
      return m;
  }
}

//! @brief Converts the copula parameters to Blomqvist's \f$ \beta \f$.
//!
//! @details Blomqvist's beta is computed from the copula cdf as
//! \f$ \beta = 4\, C(0.5, 0.5) - 1 \f$, using the same formula for every
//! family (including the nonparametric `tll`).
//!
//! @param parameters The parameters (must be a valid parametrization of
//!     the current family).
inline double
Bicop::parameters_to_beta(const Eigen::MatrixXd& parameters) const
{
  double beta = bicop_->parameters_to_beta(parameters);
  if (tools_stl::is_member(rotation_, { 90, 270 })) {
    beta *= -1;
  }
  return beta;
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

//! @brief Gets the tail dependence coefficients.
//!
//! @details See Bicop::parameters_to_taildep() for the layout of the returned
//! \f$ 2 \times 2 \f$ matrix.
inline Eigen::MatrixXd
Bicop::get_taildep() const
{
  return parameters_to_taildep(bicop_->get_parameters());
}

//! @brief Gets Blomqvist's beta.
inline double
Bicop::get_beta() const
{
  return parameters_to_beta(bicop_->get_parameters());
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
  if (bicop_->get_family_name() != "Independence") {
    bicop_->set_parameters(parameters);
  }
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
  // change var_types
  // ----------------
  std::swap(var_types_[0], var_types_[1]);
  flip_abstract_var_types();

  // change internal representation
  // ------------------------------
  // For families that are exchangeable (most in our library), flipping 0° or
  // 180° rotated versions does not alter the shape of the copula.
  //
  // For 90° and 270° rotations, their definition itself involves exchanging
  // arguments plus a reflection in one of the arguments (e.g.,
  // c(u, v) -> c(1 - v, u)). To indicate that the the argument which needs to
  // be reflected changes, we always switch 90° <-> 270° . For families that are
  // exchangeable in their 0° variant, this will be the only change this is
  // needed.
  if (rotation_ == 90) {
    rotation_ = 270;
  } else if (rotation_ == 270) {
    rotation_ = 90;
  }
  // The following implements any changes to the shape beyond the change in
  // rotation. Formost of our families, it does nothing.
  bicop_->flip();
}

//! @brief Summarizes the model into a string (can be used for printing).
inline std::string
Bicop::str() const
{
  std::stringstream bicop_str;
  bicop_str << std::setprecision(2); // set precision to 2 decimal places
  bicop_str << "Bivariate copula: \n";
  bicop_str << "  family = " << get_family_name() << "\n";
  bicop_str << "  rotation = " << get_rotation() << "\n";
  bicop_str << "  var_types = " << var_types_[0] << "," << var_types_[1]
            << "\n";
  if (get_family() == BicopFamily::tll) {
    bicop_str << "  parameters = [30x30 grid] with " << get_npars()
              << " d.f.\n";
  } else if (get_family() != BicopFamily::indep) {
    bicop_str << "  parameters = " << get_parameters() << "\n";
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

//! @brief Fits a bivariate copula (with fixed family) to data.
//!
//! @details For parametric models, two different methods are available. `"mle"`
//! fits the parameters by maximum-likelihood. `"itau"` uses inversion of
//! Kendall's \f$ \tau \f$, but is only available for one-parameter families
//! and the Student t copula. For the latter, there is a one-to-one
//! transformation for the first parameter, the second is found by profile
//! likelihood optimization (with accuracy of at least 0.5). Nonparametric
//! families have specialized methods, no specification is required.
//!
//! When at least one variable is discrete, two types of "observations"
//! are required: the first \f$ n \times 2 \f$ block contains realizations of
//! \f$ F_{X_1}(X_1), F_{X_2}(X_2) \f$. Let \f$ k \f$ denote the number of
//! discrete variables (either one or two). Then the second \f$ n \times k \f$
//! block contains realizations of \f$ F_{X_k}(X_k^-) \f$. The minus indicates a
//! left-sided limit of the cdf. For continuous variables the left limit and the
//! cdf itself coincide. For, e.g., an integer-valued variable, it holds \f$
//! F_{X_k}(X_k^-) = F_{X_k}(X_k - 1) \f$.
//!
//! Incomplete observations (i.e., ones with a NaN value) are discarded.
//!
//! @param data An \f$ n \times (2 + k) \f$ matrix of observations contained in
//!   \f$(0, 1) \f$, where \f$ k \f$ is the number of discrete variables.
//! @param controls The controls (see `FitControlsBicop`).
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
              controls.get_nonparametric_grid_size(),
              w);
  nobs_ = data_no_nan.rows();
}

//

//! @brief Selects the best fitting model.
//!
//! @details The function calls `Bicop::fit()` for all families in
//! `family_set` and selecting the best fitting model by either BIC or AIC,
//! see `Bicop::bic()` and `Bicop::aic()`.
//!
//! When at least one variable is discrete, two types of "observations"
//! are required: the first \f$ n \times 2 \f$ block contains realizations of
//! \f$ F_{X_1}(X_1), F_{X_2}(X_2) \f$. Let \f$ k \f$ denote the number of
//! discrete variables (either one or two). Then the second \f$ n \times k \f$
//! block contains realizations of \f$ F_{X_k}(X_k^-) \f$. The minus indicates a
//! left-sided limit of the cdf. For continuous variables the left limit and the
//! cdf itself coincide. For, e.g., an integer-valued variable, it holds \f$
//! F_{X_k}(X_k^-) = F_{X_k}(X_k - 1) \f$.
//!
//! Incomplete observations (i.e., ones with a NaN value) are discarded.
//!
//! @param data An \f$ n \times (2 + k) \f$ matrix of observations contained in
//!   \f$(0, 1) \f$, where \f$ k \f$ is the number of discrete variables.
//! @param controls The controls (see `FitControlsBicop`).
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
        if (controls.get_selection_criterion() == "mbic" ||
            controls.get_selection_criterion() == "mbicv") {
          // correction for mBIC or mBICV
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
    pool.wait();
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
  if ((weights.size() > 0) && (weights.size() != data.rows())) {
    throw std::runtime_error("sizes of weights and data don't match.");
  }
}

//! @brief Checks whether the Bicop object was fitted to data.
inline void
Bicop::check_fitted() const
{
  if ((std::isnan)(bicop_->get_loglik())) {
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
