// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <memory>

#include <Eigen/Dense>
#include <vinecopulib/bicop/family.hpp>

namespace vinecopulib {

//! @brief Utilities for parsing derivative selectors.
//!
//! A selector is a concatenation of components `"par<k>"` (k-th parameter,
//! 1-based; `"par"` is short for `"par1"`), `"u1"`, and `"u2"`. Components
//! are encoded as integers: the 0-based parameter index for parameters, `-1`
//! for `"u1"`, and `-2` for `"u2"`.
namespace tools_deriv {

std::vector<int>
parse_components(const std::string& deriv);

std::string
comp_to_string(int comp);

std::string
components_to_string(std::vector<int> comps);

std::string
canonicalize(const std::string& deriv, size_t order, size_t npars);

std::string
swap_args(const std::string& deriv);

bool
is_u2_only(const std::string& deriv);
}

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
                   size_t grid_size,
                   const Eigen::VectorXd& weights) = 0;

  virtual double get_npars() const = 0;

  virtual void set_npars(const double& npars) = 0;

  virtual double parameters_to_tau(const Eigen::MatrixXd& parameters) = 0;

  // following two are non-pure: tail dependence defaults to NaN ("not
  // implemented") and Blomqvist's beta has a generic implementation via the
  // cdf.
  virtual Eigen::MatrixXd parameters_to_taildep(
    const Eigen::MatrixXd& parameters);

  virtual double parameters_to_beta(const Eigen::MatrixXd& parameters);

  virtual void flip() = 0;

  // following are virtual so they can be overriden by KernelBicop
  virtual Eigen::VectorXd pdf(const Eigen::MatrixXd& u);

  virtual Eigen::VectorXd hfunc1(const Eigen::MatrixXd& u);

  virtual Eigen::VectorXd hfunc2(const Eigen::MatrixXd& u);

  Eigen::VectorXd hinv1(const Eigen::MatrixXd& u);

  Eigen::VectorXd hinv2(const Eigen::MatrixXd& u);

  // Evaluation leaves. `parameters` has shape m x p with m in {1, n}: row i
  // holds the p parameters for observation i; a single row is broadcast to all
  // observations. The state-based dispatchers above call these with the stored
  // parameters (a 1 x p broadcast row). This is the sole evaluation interface;
  // every family implements the (stateless) math, and nonparametric families
  // ignore `parameters` and read their interpolation grid instead.
  Eigen::VectorXd pdf(const Eigen::MatrixXd& u,
                      const Eigen::MatrixXd& parameters);

  Eigen::VectorXd hfunc1(const Eigen::MatrixXd& u,
                         const Eigen::MatrixXd& parameters);

  Eigen::VectorXd hfunc2(const Eigen::MatrixXd& u,
                         const Eigen::MatrixXd& parameters);

  Eigen::VectorXd hinv1(const Eigen::MatrixXd& u,
                        const Eigen::MatrixXd& parameters);

  Eigen::VectorXd hinv2(const Eigen::MatrixXd& u,
                        const Eigen::MatrixXd& parameters);

  virtual Eigen::VectorXd cdf(const Eigen::MatrixXd& u,
                              const Eigen::MatrixXd& parameters) = 0;

  virtual Eigen::VectorXd pdf_raw(const Eigen::MatrixXd& u,
                                  const Eigen::MatrixXd& parameters) = 0;

  virtual Eigen::VectorXd hfunc1_raw(const Eigen::MatrixXd& u,
                                     const Eigen::MatrixXd& parameters) = 0;

  virtual Eigen::VectorXd hfunc2_raw(const Eigen::MatrixXd& u,
                                     const Eigen::MatrixXd& parameters) = 0;

  virtual Eigen::VectorXd hinv1_raw(const Eigen::MatrixXd& u,
                                    const Eigen::MatrixXd& parameters) = 0;

  virtual Eigen::VectorXd hinv2_raw(const Eigen::MatrixXd& u,
                                    const Eigen::MatrixXd& parameters) = 0;

  // Derivative leaves. `deriv` is a canonical selector (see tools_deriv):
  // first order `"par1"`, `"par2"`, ..., `"u1"`, `"u2"`; second order a
  // sorted concatenation like `"par1u1"` (parameters first, then `"u1"`,
  // then `"u2"`). The facade canonicalizes user input and resolves rotations
  // before calling these, so implementations only see canonical selectors
  // and 0-degree-rotation data. The defaults here throw (nonparametric
  // families); ParBicop overrides them with central finite differences so
  // every parametric family works; families with closed forms override
  // those in turn. The logpdf defaults compose the pdf leaves by the
  // quotient rule; families override the parameter selectors where dedicated
  // closed forms exist.
  virtual Eigen::VectorXd pdf_deriv_raw(const Eigen::MatrixXd& u,
                                        const Eigen::MatrixXd& parameters,
                                        const std::string& deriv);

  virtual Eigen::VectorXd pdf_deriv2_raw(const Eigen::MatrixXd& u,
                                         const Eigen::MatrixXd& parameters,
                                         const std::string& deriv);

  virtual Eigen::VectorXd hfunc1_deriv_raw(const Eigen::MatrixXd& u,
                                           const Eigen::MatrixXd& parameters,
                                           const std::string& deriv);

  virtual Eigen::VectorXd hfunc1_deriv2_raw(const Eigen::MatrixXd& u,
                                            const Eigen::MatrixXd& parameters,
                                            const std::string& deriv);

  virtual Eigen::VectorXd hfunc2_deriv_raw(const Eigen::MatrixXd& u,
                                           const Eigen::MatrixXd& parameters,
                                           const std::string& deriv);

  virtual Eigen::VectorXd hfunc2_deriv2_raw(const Eigen::MatrixXd& u,
                                            const Eigen::MatrixXd& parameters,
                                            const std::string& deriv);

  virtual Eigen::VectorXd logpdf_deriv_raw(const Eigen::MatrixXd& u,
                                           const Eigen::MatrixXd& parameters,
                                           const std::string& deriv);

  virtual Eigen::VectorXd logpdf_deriv2_raw(const Eigen::MatrixXd& u,
                                            const Eigen::MatrixXd& parameters,
                                            const std::string& deriv);

  virtual Eigen::MatrixXd tau_to_parameters(const double& tau) = 0;
  Eigen::MatrixXd no_tau_to_parameters(const double&);

  // Misc methods
  Eigen::VectorXd hinv1_num(const Eigen::MatrixXd& u);

  Eigen::VectorXd hinv2_num(const Eigen::MatrixXd& u);

  Eigen::VectorXd hinv1_num(const Eigen::MatrixXd& u,
                            const Eigen::MatrixXd& parameters);

  Eigen::VectorXd hinv2_num(const Eigen::MatrixXd& u,
                            const Eigen::MatrixXd& parameters);

  // continuous numeric inverses: invert the *_raw leaves (which ignore
  // var_types_), used by the `_raw` primitives so a continuous inverse never
  // routes through the discrete h-function dispatcher
  Eigen::VectorXd hinv1_num_raw(const Eigen::MatrixXd& u,
                                const Eigen::MatrixXd& parameters);

  Eigen::VectorXd hinv2_num_raw(const Eigen::MatrixXd& u,
                                const Eigen::MatrixXd& parameters);

  Eigen::VectorXd pdf_c_d(const Eigen::MatrixXd& u,
                          const Eigen::MatrixXd& parameters);

  Eigen::VectorXd pdf_d_d(const Eigen::MatrixXd& u,
                          const Eigen::MatrixXd& parameters);

  double loglik(const Eigen::MatrixXd& u,
                const Eigen::VectorXd& weights = Eigen::VectorXd());

  // Data members
  BicopFamily family_;
  double loglik_{ NAN };
  std::vector<std::string> var_types_{ "c", "c" };
};

//! A shared pointer to an object of class AbstracBicop.
using BicopPtr = std::shared_ptr<AbstractBicop>;
}

#include <vinecopulib/bicop/implementation/abstract.ipp>
