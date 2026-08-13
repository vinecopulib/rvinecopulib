// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <algorithm>
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
#include <vinecopulib/bicop/tawn.hpp>
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
    case BicopFamily::tawn:
      new_bicop = BicopPtr(new TawnBicop());
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

//! Default tail dependence: not implemented for this family, so all four
//! corners are reported as NaN. Families with a closed form override this
//! (including those that genuinely have zero tail dependence, e.g. `indep`,
//! `gaussian`, `frank`).
inline Eigen::MatrixXd
AbstractBicop::parameters_to_taildep(const Eigen::MatrixXd&)
{
  return Eigen::MatrixXd::Constant(2, 2, NAN);
}

//! Blomqvist's beta computed generically from the copula cdf as
//! \f$ \beta = 4 C(0.5, 0.5) - 1 \f$. Works for all families (including the
//! nonparametric kernel estimator).
inline double
AbstractBicop::parameters_to_beta(const Eigen::MatrixXd& parameters)
{
  Eigen::MatrixXd u(1, 2);
  u << 0.5, 0.5;
  // the cdf leaf expects an (m x p) parameter matrix (one row per evaluation
  // point); parameters is a (p x 1) column, so transpose it to a single row.
  Eigen::MatrixXd par_row =
    (parameters.cols() == 1) ? parameters.transpose() : parameters;
  return 4.0 * this->cdf(u, par_row)(0) - 1.0;
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
  if (var_types_ != std::vector<std::string>{ "c", "c" }) {
    // discrete margins go through the parameter-aware difference quotients
    return pdf(u, get_parameters().transpose());
  }
  Eigen::VectorXd pdf = pdf_raw(u.leftCols(2), get_parameters().transpose());
  tools_eigen::trim(pdf, DBL_MIN, DBL_MAX);
  return pdf;
}

inline Eigen::VectorXd
AbstractBicop::hfunc1(const Eigen::MatrixXd& u)
{
  if (var_types_[0] == "d") {
    return hfunc1(u, get_parameters().transpose());
  }
  return hfunc1_raw(u.leftCols(2), get_parameters().transpose());
}

inline Eigen::VectorXd
AbstractBicop::hfunc2(const Eigen::MatrixXd& u)
{
  if (var_types_[1] == "d") {
    return hfunc2(u, get_parameters().transpose());
  }
  return hfunc2_raw(u.leftCols(2), get_parameters().transpose());
}

inline Eigen::VectorXd
AbstractBicop::hinv1(const Eigen::MatrixXd& u)
{
  if (var_types_[0] == "c") {
    return hinv1_raw(u.leftCols(2), get_parameters().transpose());
  } else {
    return hinv1_num(u);
  }
}

inline Eigen::VectorXd
AbstractBicop::hinv2(const Eigen::MatrixXd& u)
{
  if (var_types_[1] == "c") {
    return hinv2_raw(u.leftCols(2), get_parameters().transpose());
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
  return hinv1_num(u, get_parameters().transpose());
}

inline Eigen::VectorXd
AbstractBicop::hinv2_num(const Eigen::MatrixXd& u)
{
  return hinv2_num(u, get_parameters().transpose());
}
//! @}

//! @name Parameter-aware leaves and dispatchers
//!
//! `parameters` has shape m x p with m in {1, n}: row i holds the p parameters
//! for observation i; a single row is broadcast to all observations. These are
//! the sole evaluation interface; the state-based methods above call them with
//! the stored parameters as a 1 x p broadcast row.
//! @{

inline Eigen::VectorXd
AbstractBicop::pdf(const Eigen::MatrixXd& u, const Eigen::MatrixXd& parameters)
{
  Eigen::VectorXd pdf(u.rows());
  if (var_types_ == std::vector<std::string>{ "c", "c" }) {
    pdf = pdf_raw(u.leftCols(2), parameters);
  } else if (var_types_ == std::vector<std::string>{ "d", "d" }) {
    pdf = pdf_d_d(u, parameters);
  } else {
    pdf = pdf_c_d(u, parameters);
  }
  tools_eigen::trim(pdf, DBL_MIN, DBL_MAX);
  return pdf;
}

inline Eigen::VectorXd
AbstractBicop::pdf_c_d(const Eigen::MatrixXd& u,
                       const Eigen::MatrixXd& parameters)
{
  Eigen::VectorXd pdf(u.rows());
  Eigen::MatrixXd umax = u.leftCols(2);
  Eigen::MatrixXd umin = u.rightCols(2);
  Eigen::VectorXd udiff(u.rows());

  if (var_types_[0] != "c") {
    udiff = (u.col(0) - u.col(2)).cwiseAbs();
  } else {
    udiff = (u.col(1) - u.col(3)).cwiseAbs();
  }

  const bool bc = (parameters.rows() != u.rows());
  for (Eigen::Index i = 0; i < u.rows(); i++) {
    const Eigen::MatrixXd par_i =
      parameters.rows() == 0 ? parameters : parameters.row(bc ? 0 : i);
    if (udiff(i) > 5e-5) {
      if (var_types_[0] != "c") {
        pdf(i) =
          (hfunc2_raw(umax.row(i), par_i) - hfunc2_raw(umin.row(i), par_i))(0);
      } else {
        pdf(i) =
          (hfunc1_raw(umax.row(i), par_i) - hfunc1_raw(umin.row(i), par_i))(0);
      }
      pdf(i) /= udiff(i);
    } else {
      pdf(i) = pdf_raw((umax.row(i) + umin.row(i)) / 2, par_i)(0);
    }
  }
  return pdf.cwiseAbs();
}

inline Eigen::VectorXd
AbstractBicop::pdf_d_d(const Eigen::MatrixXd& u,
                       const Eigen::MatrixXd& parameters)
{
  Eigen::VectorXd pdf(u.rows());
  Eigen::MatrixXd umax = u.leftCols(2);
  Eigen::MatrixXd umin = u.rightCols(2);
  Eigen::MatrixXd udiff = (umax - umin).cwiseAbs();

  const bool bc = (parameters.rows() != u.rows());
  for (Eigen::Index i = 0; i < u.rows(); i++) {
    const Eigen::MatrixXd par_i =
      parameters.rows() == 0 ? parameters : parameters.row(bc ? 0 : i);
    // the difference quotient can be instable, use derivative if denominator
    // too small
    if (udiff.row(i).maxCoeff() < 5e-5) {
      pdf(i) = pdf_raw((umax.row(i) + umin.row(i)) / 2, par_i)(0);
    } else if (udiff(i, 0) < 5e-5) {
      umax(i, 0) = (umax(i, 0) + umin(i, 0)) / 2;
      umin(i, 0) = (umax(i, 0) + umin(i, 0)) / 2;
      pdf(i) = (hfunc1_raw(umax.row(i), par_i)(0) -
                hfunc1_raw(umin.row(i), par_i)(0)) /
               udiff(i, 1);
    } else if (udiff(i, 1) < 5e-5) {
      umax(i, 1) = (umax(i, 1) + umin(i, 1)) / 2;
      umin(i, 1) = (umax(i, 1) + umin(i, 1)) / 2;
      pdf(i) = (hfunc2_raw(umax.row(i), par_i)(0) -
                hfunc2_raw(umin.row(i), par_i)(0)) /
               udiff(i, 0);
    } else {
      pdf(i) = cdf(umax.row(i), par_i)(0) + cdf(umin.row(i), par_i)(0);
      std::swap(umax(i, 0), umin(i, 0));
      pdf(i) -= cdf(umax.row(i), par_i)(0) + cdf(umin.row(i), par_i)(0);
      pdf(i) /= udiff(i, 0) * udiff(i, 1);
    }
  }

  return pdf.cwiseAbs();
}

inline Eigen::VectorXd
AbstractBicop::hfunc1(const Eigen::MatrixXd& u,
                      const Eigen::MatrixXd& parameters)
{
  if (var_types_[0] == "d") {
    auto uu = u;
    uu.col(3) = uu.col(1);
    auto u1diff = (uu.col(0) - uu.col(2)).cwiseAbs();
    Eigen::VectorXd h(u.rows());

    const bool bc = (parameters.rows() != u.rows());
    for (Eigen::Index i = 0; i < u.rows(); i++) {
      const Eigen::MatrixXd par_i =
        parameters.rows() == 0 ? parameters : parameters.row(bc ? 0 : i);
      if (std::abs(u1diff(i)) > 5e-5) {
        h(i) = cdf(uu.row(i).leftCols(2), par_i)(0) -
               cdf(uu.row(i).rightCols(2), par_i)(0);
        h(i) /= u1diff(i);
      } else {
        uu(i, 0) = (uu(i, 0) + uu(i, 2)) / 2;
        h(i) = hfunc1_raw(uu.row(i).leftCols(2), par_i)(0);
      }
    }
    return h.cwiseAbs();
  } else {
    return hfunc1_raw(u.leftCols(2), parameters);
  }
}

inline Eigen::VectorXd
AbstractBicop::hfunc2(const Eigen::MatrixXd& u,
                      const Eigen::MatrixXd& parameters)
{
  if (var_types_[1] == "d") {
    auto uu = u;
    uu.col(2) = uu.col(0);
    auto u2diff = (uu.col(1) - uu.col(3)).cwiseAbs();
    Eigen::VectorXd h(u.rows());

    const bool bc = (parameters.rows() != u.rows());
    for (Eigen::Index i = 0; i < u.rows(); i++) {
      const Eigen::MatrixXd par_i =
        parameters.rows() == 0 ? parameters : parameters.row(bc ? 0 : i);
      if (u2diff(i) > 5e-5) {
        h(i) = cdf(uu.row(i).leftCols(2), par_i)(0) -
               cdf(uu.row(i).rightCols(2), par_i)(0);
        h(i) /= u2diff(i);
      } else {
        uu(i, 1) = (uu(i, 1) + uu(i, 3)) / 2;
        h(i) = hfunc2_raw(uu.row(i).leftCols(2), par_i)(0);
      }
    }
    return h.cwiseAbs();
  } else {
    return hfunc2_raw(u.leftCols(2), parameters);
  }
}

inline Eigen::VectorXd
AbstractBicop::hinv1(const Eigen::MatrixXd& u,
                     const Eigen::MatrixXd& parameters)
{
  if (var_types_[0] == "c") {
    return hinv1_raw(u.leftCols(2), parameters);
  } else {
    return hinv1_num(u, parameters);
  }
}

inline Eigen::VectorXd
AbstractBicop::hinv2(const Eigen::MatrixXd& u,
                     const Eigen::MatrixXd& parameters)
{
  if (var_types_[1] == "c") {
    return hinv2_raw(u.leftCols(2), parameters);
  } else {
    return hinv2_num(u, parameters);
  }
}

inline Eigen::VectorXd
AbstractBicop::hinv1_num(const Eigen::MatrixXd& u,
                         const Eigen::MatrixXd& parameters)
{
  Eigen::MatrixXd u_new = u;
  auto h1 = [&](const Eigen::VectorXd& v) {
    u_new.col(1) = v;
    return hfunc1(u_new, parameters);
  };

  return tools_eigen::invert_f(u.col(1), h1);
}

inline Eigen::VectorXd
AbstractBicop::hinv2_num(const Eigen::MatrixXd& u,
                         const Eigen::MatrixXd& parameters)
{
  Eigen::MatrixXd u_new = u;
  auto h1 = [&](const Eigen::VectorXd& x) {
    u_new.col(0) = x;
    return hfunc2(u_new, parameters);
  };

  return tools_eigen::invert_f(u.col(0), h1);
}

inline Eigen::VectorXd
AbstractBicop::hinv1_num_raw(const Eigen::MatrixXd& u,
                             const Eigen::MatrixXd& parameters)
{
  Eigen::MatrixXd u_new = u;
  auto h1 = [&](const Eigen::VectorXd& v) {
    u_new.col(1) = v;
    return hfunc1_raw(u_new.leftCols(2), parameters);
  };

  return tools_eigen::invert_f(u.col(1), h1);
}

inline Eigen::VectorXd
AbstractBicop::hinv2_num_raw(const Eigen::MatrixXd& u,
                             const Eigen::MatrixXd& parameters)
{
  Eigen::MatrixXd u_new = u;
  auto h1 = [&](const Eigen::VectorXd& x) {
    u_new.col(0) = x;
    return hfunc2_raw(u_new.leftCols(2), parameters);
  };

  return tools_eigen::invert_f(u.col(0), h1);
}
//! @}

//! @name Derivative leaves (defaults)
//!
//! Nonparametric families have no analytical derivatives; ParBicop overrides
//! these with finite-difference fallbacks and analytic families override
//! those in turn.
//! @{

inline Eigen::VectorXd
AbstractBicop::pdf_deriv_raw(const Eigen::MatrixXd&,
                             const Eigen::MatrixXd&,
                             const std::string&)
{
  throw std::runtime_error("derivatives are not implemented for the " +
                           get_family_name() + " copula");
}

inline Eigen::VectorXd
AbstractBicop::pdf_deriv2_raw(const Eigen::MatrixXd&,
                              const Eigen::MatrixXd&,
                              const std::string&)
{
  throw std::runtime_error("derivatives are not implemented for the " +
                           get_family_name() + " copula");
}

inline Eigen::VectorXd
AbstractBicop::hfunc1_deriv_raw(const Eigen::MatrixXd&,
                                const Eigen::MatrixXd&,
                                const std::string&)
{
  throw std::runtime_error("derivatives are not implemented for the " +
                           get_family_name() + " copula");
}

inline Eigen::VectorXd
AbstractBicop::hfunc1_deriv2_raw(const Eigen::MatrixXd&,
                                 const Eigen::MatrixXd&,
                                 const std::string&)
{
  throw std::runtime_error("derivatives are not implemented for the " +
                           get_family_name() + " copula");
}

inline Eigen::VectorXd
AbstractBicop::hfunc2_deriv_raw(const Eigen::MatrixXd&,
                                const Eigen::MatrixXd&,
                                const std::string&)
{
  throw std::runtime_error("derivatives are not implemented for the " +
                           get_family_name() + " copula");
}

inline Eigen::VectorXd
AbstractBicop::hfunc2_deriv2_raw(const Eigen::MatrixXd&,
                                 const Eigen::MatrixXd&,
                                 const std::string&)
{
  throw std::runtime_error("derivatives are not implemented for the " +
                           get_family_name() + " copula");
}

inline Eigen::VectorXd
AbstractBicop::logpdf_deriv_raw(const Eigen::MatrixXd& u,
                                const Eigen::MatrixXd& parameters,
                                const std::string& deriv)
{
  Eigen::ArrayXd c = pdf_raw(u, parameters).array().max(DBL_MIN);
  return (pdf_deriv_raw(u, parameters, deriv).array() / c).matrix();
}

inline Eigen::VectorXd
AbstractBicop::logpdf_deriv2_raw(const Eigen::MatrixXd& u,
                                 const Eigen::MatrixXd& parameters,
                                 const std::string& deriv)
{
  auto comps = tools_deriv::parse_components(deriv);
  Eigen::ArrayXd c = pdf_raw(u, parameters).array().max(DBL_MIN);
  Eigen::ArrayXd c_xy = pdf_deriv2_raw(u, parameters, deriv).array();
  Eigen::ArrayXd c_x =
    pdf_deriv_raw(u, parameters, tools_deriv::comp_to_string(comps[0])).array();
  Eigen::ArrayXd c_y =
    (comps[0] == comps[1])
      ? c_x
      : pdf_deriv_raw(u, parameters, tools_deriv::comp_to_string(comps[1]))
          .array();
  return (c_xy / c - (c_x / c) * (c_y / c)).matrix();
}
//! @}

namespace tools_deriv {

//! splits a derivative selector into components; the encoding is the 0-based
//! parameter index for `"par<k>"`, `-1` for `"u1"`, and `-2` for `"u2"`.
inline std::vector<int>
parse_components(const std::string& deriv)
{
  std::vector<int> comps;
  size_t pos = 0;
  while (pos < deriv.size()) {
    if (deriv.compare(pos, 1, "u") == 0) {
      if ((pos + 1 >= deriv.size()) ||
          ((deriv[pos + 1] != '1') && (deriv[pos + 1] != '2'))) {
        throw std::runtime_error("invalid derivative selector: '" + deriv +
                                 "'");
      }
      comps.push_back((deriv[pos + 1] == '1') ? -1 : -2);
      pos += 2;
    } else if (deriv.compare(pos, 3, "par") == 0) {
      pos += 3;
      size_t k = 0;
      size_t num_digits = 0;
      while ((pos < deriv.size()) && (deriv[pos] >= '0') &&
             (deriv[pos] <= '9')) {
        k = 10 * k + static_cast<size_t>(deriv[pos] - '0');
        pos++;
        num_digits++;
      }
      if (num_digits == 0) {
        k = 1; // "par" is short for "par1"
      }
      if ((k < 1) || (num_digits > 3)) {
        throw std::runtime_error("invalid derivative selector: '" + deriv +
                                 "'");
      }
      comps.push_back(static_cast<int>(k) - 1);
    } else {
      throw std::runtime_error("invalid derivative selector: '" + deriv + "'");
    }
  }
  if (comps.empty()) {
    throw std::runtime_error("derivative selector cannot be empty");
  }
  return comps;
}

inline std::string
comp_to_string(int comp)
{
  if (comp == -1) {
    return "u1";
  } else if (comp == -2) {
    return "u2";
  }
  return "par" + std::to_string(comp + 1);
}

//! sorts components canonically (parameters by index, then u1, then u2) and
//! concatenates them.
inline std::string
components_to_string(std::vector<int> comps)
{
  auto rank = [](int comp) { return (comp >= 0) ? comp : 1000 - comp; };
  std::sort(comps.begin(), comps.end(), [&](int lhs, int rhs) {
    return rank(lhs) < rank(rhs);
  });
  std::string str;
  for (auto comp : comps) {
    str += comp_to_string(comp);
  }
  return str;
}

//! validates a user-facing selector against the derivative order and the
//! number of parameters and returns its canonical form; a single component
//! in a second-order selector means differentiating twice w.r.t. it.
inline std::string
canonicalize(const std::string& deriv, size_t order, size_t npars)
{
  auto comps = parse_components(deriv);
  if ((comps.size() == 1) && (order == 2)) {
    comps.push_back(comps[0]);
  }
  if (comps.size() != order) {
    throw std::runtime_error("derivative selector '" + deriv + "' has " +
                             std::to_string(comps.size()) +
                             " components; expected " + std::to_string(order));
  }
  for (auto comp : comps) {
    if (comp >= static_cast<int>(npars)) {
      throw std::runtime_error(
        "derivative selector '" + deriv + "' refers to parameter " +
        std::to_string(comp + 1) + ", but the family has " +
        std::to_string(npars) + " parameter(s)");
    }
  }
  return components_to_string(comps);
}

//! swaps `"u1"` and `"u2"` in a selector (for exchangeable families).
inline std::string
swap_args(const std::string& deriv)
{
  auto comps = parse_components(deriv);
  for (auto& comp : comps) {
    if (comp == -1) {
      comp = -2;
    } else if (comp == -2) {
      comp = -1;
    }
  }
  return components_to_string(comps);
}

//! whether a selector involves `"u2"` but not `"u1"` (an exchangeable family
//! can then route it to its `"u1"`-flavored leaf via `swap_args`).
inline bool
is_u2_only(const std::string& deriv)
{
  auto comps = parse_components(deriv);
  bool has_u1 = std::find(comps.begin(), comps.end(), -1) != comps.end();
  bool has_u2 = std::find(comps.begin(), comps.end(), -2) != comps.end();
  return has_u2 && !has_u1;
}
}
}
