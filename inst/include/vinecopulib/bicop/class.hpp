// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <vinecopulib/bicop/fit_controls.hpp>
#include <vinecopulib/misc/tools_serialization.hpp>

namespace vinecopulib {

// forward declaration of Abstract class
class AbstractBicop;
class BicopView;
using BicopPtr = std::shared_ptr<AbstractBicop>;

//! @brief A class for bivariate copula models.
//!
//! @details The model is fully characterized by the family,
//! rotation (one of `0`, `90`, `180`, `270`), a matrix of parameters, and
//! variable types (two strings, one for each variable, either `"c"`
//! for continuous or `"d"` for discrete).
//!
//! Implemented families (see `BicopFamily`):
//!
//! ```
//! | type          | full name             | string identifier     |
//! |---------------|-----------------------|-----------------------|
//! | -             | Independence          | "indep"               |
//! | Elliptical    | Gaussian              | "gaussian"            |
//! |               | Student t             | "student"             |
//! | Archimedean   | Clayton               | "clayton"             |
//! |               | Gumbel                | "gumbel"              |
//! |               | Frank                 | "frank"               |
//! |               | Joe                   | "joe"                 |
//! |               | Clayton-Gumbel (BB1)  | "bb1"                 |
//! |               | Joe-Gumbel (BB6)      | "bb6"                 |
//! |               | Joe-Clayton (BB7)     | "bb7"                 |
//! |               | Joe-Frank (BB8)       | "bb8"                 |
//! | Extreme-Value | Tawn                  | "tawn"                |
//! | Nonparametric | Transformation kernel | "tll"                 |
//! ```
//!
class Bicop
{
  friend class BicopView;

public:
  // Constructors
  Bicop(const BicopFamily family = BicopFamily::indep,
        const int rotation = 0,
        const Eigen::MatrixXd& parameters = Eigen::MatrixXd(),
        const std::vector<std::string>& var_types = { "c", "c" });

  explicit Bicop(const Eigen::MatrixXd& data,
                 const FitControlsBicop& controls = FitControlsBicop(),
                 const std::vector<std::string>& var_types = { "c", "c" });

  Bicop(const Bicop& other);

  // Required, not redundant: the copy constructor above suppresses the implicit
  // move, and copying deep-clones the family object.
  Bicop(Bicop&& other) = default;

  explicit Bicop(const std::string& filename);

  explicit Bicop(const nlohmann::json& input);

  Bicop& operator=(Bicop other);

  // Serialize
  nlohmann::json to_json() const;

  void to_file(const std::string& filename) const;

  // Getters and setters

  //! @return the copula family.
  BicopFamily get_family() const;

  //! @return the human-readable name of the copula family.
  std::string get_family_name() const;

  //! @return the copula rotation in degrees (one of `0`, `90`, `180`, `270`).
  int get_rotation() const;

  //! @return the copula parameter(s) as a matrix (shape depends on the family).
  Eigen::MatrixXd get_parameters() const;

  //! @return Kendall's tau implied by the current parameters.
  double get_tau() const;

  //! @return the tail dependence coefficients implied by the current
  //! parameters, as a 2x2 matrix; see parameters_to_taildep().
  Eigen::MatrixXd get_taildep() const;

  //! @return Blomqvist's beta implied by the current parameters.
  double get_beta() const;

  //! @return the number of parameters in the copula model. For nonparametric
  //! families, this is a conceptually similar effective-parameter count.
  double get_npars() const;

  //! @return the log-likelihood of the data the copula was fitted to.
  double get_loglik() const;
  //! @return the number of observations used to fit the copula (`0` until
  //! `fit()` has been called).
  size_t get_nobs() const;
  //! @return Akaike's information criterion at the fitted parameters.
  double get_aic() const;
  //! @return the Bayesian information criterion at the fitted parameters.
  double get_bic() const;
  //! @return the modified Bayesian information criterion at the fitted
  //! parameters; `psi0` is the prior probability of a non-independence
  //! copula.
  double get_mbic(const double psi0 = 0.9) const;

  //! Sets the copula rotation.
  //! @param rotation one of `0`, `90`, `180`, `270`.
  void set_rotation(const int rotation);

  //! Sets the copula parameters.
  //! @param parameters parameter matrix; shape must match the family's
  //! expected parameter layout.
  void set_parameters(const Eigen::MatrixXd& parameters);

  //! Sets the variable types.
  //! @param var_types a length-2 vector with each entry either `"c"`
  //! (continuous, default) or `"d"` (discrete).
  void set_var_types(const std::vector<std::string>& var_types = { "c", "c" });

  //! @return the variable types of the two variables (each `"c"` for
  //! continuous or `"d"` for discrete).
  std::vector<std::string> get_var_types() const;

  // Stats methods
  Eigen::VectorXd pdf(const Eigen::MatrixXd& u) const;

  Eigen::VectorXd cdf(const Eigen::MatrixXd& u) const;

  Eigen::VectorXd hfunc1(const Eigen::MatrixXd& u) const;

  Eigen::VectorXd hfunc2(const Eigen::MatrixXd& u) const;

  Eigen::VectorXd hinv1(const Eigen::MatrixXd& u) const;

  Eigen::VectorXd hinv2(const Eigen::MatrixXd& u) const;

  // Stats methods with per-row parameters (parametric families only)
  Eigen::VectorXd pdf(const Eigen::MatrixXd& u,
                      const Eigen::MatrixXd& parameters,
                      const size_t num_threads = 1) const;

  Eigen::VectorXd cdf(const Eigen::MatrixXd& u,
                      const Eigen::MatrixXd& parameters,
                      const size_t num_threads = 1) const;

  Eigen::VectorXd hfunc1(const Eigen::MatrixXd& u,
                         const Eigen::MatrixXd& parameters,
                         const size_t num_threads = 1) const;

  Eigen::VectorXd hfunc2(const Eigen::MatrixXd& u,
                         const Eigen::MatrixXd& parameters,
                         const size_t num_threads = 1) const;

  Eigen::VectorXd hinv1(const Eigen::MatrixXd& u,
                        const Eigen::MatrixXd& parameters,
                        const size_t num_threads = 1) const;

  Eigen::VectorXd hinv2(const Eigen::MatrixXd& u,
                        const Eigen::MatrixXd& parameters,
                        const size_t num_threads = 1) const;

  double loglik(const Eigen::MatrixXd& u,
                const Eigen::MatrixXd& parameters,
                const size_t num_threads = 1) const;

  // Derivatives of the density and h-functions w.r.t. parameters/arguments
  Eigen::VectorXd pdf_deriv(const Eigen::MatrixXd& u,
                            const std::string& deriv) const;

  Eigen::VectorXd pdf_deriv2(const Eigen::MatrixXd& u,
                             const std::string& deriv) const;

  Eigen::VectorXd hfunc1_deriv(const Eigen::MatrixXd& u,
                               const std::string& deriv) const;

  Eigen::VectorXd hfunc1_deriv2(const Eigen::MatrixXd& u,
                                const std::string& deriv) const;

  Eigen::VectorXd hfunc2_deriv(const Eigen::MatrixXd& u,
                               const std::string& deriv) const;

  Eigen::VectorXd hfunc2_deriv2(const Eigen::MatrixXd& u,
                                const std::string& deriv) const;

  Eigen::VectorXd logpdf_deriv(const Eigen::MatrixXd& u,
                               const std::string& deriv) const;

  Eigen::VectorXd logpdf_deriv2(const Eigen::MatrixXd& u,
                                const std::string& deriv) const;

  // Derivatives with per-row parameters (parametric families only)
  Eigen::VectorXd pdf_deriv(const Eigen::MatrixXd& u,
                            const std::string& deriv,
                            const Eigen::MatrixXd& parameters,
                            const size_t num_threads = 1) const;

  Eigen::VectorXd pdf_deriv2(const Eigen::MatrixXd& u,
                             const std::string& deriv,
                             const Eigen::MatrixXd& parameters,
                             const size_t num_threads = 1) const;

  Eigen::VectorXd hfunc1_deriv(const Eigen::MatrixXd& u,
                               const std::string& deriv,
                               const Eigen::MatrixXd& parameters,
                               const size_t num_threads = 1) const;

  Eigen::VectorXd hfunc1_deriv2(const Eigen::MatrixXd& u,
                                const std::string& deriv,
                                const Eigen::MatrixXd& parameters,
                                const size_t num_threads = 1) const;

  Eigen::VectorXd hfunc2_deriv(const Eigen::MatrixXd& u,
                               const std::string& deriv,
                               const Eigen::MatrixXd& parameters,
                               const size_t num_threads = 1) const;

  Eigen::VectorXd hfunc2_deriv2(const Eigen::MatrixXd& u,
                                const std::string& deriv,
                                const Eigen::MatrixXd& parameters,
                                const size_t num_threads = 1) const;

  Eigen::VectorXd logpdf_deriv(const Eigen::MatrixXd& u,
                               const std::string& deriv,
                               const Eigen::MatrixXd& parameters,
                               const size_t num_threads = 1) const;

  Eigen::VectorXd logpdf_deriv2(const Eigen::MatrixXd& u,
                                const std::string& deriv,
                                const Eigen::MatrixXd& parameters,
                                const size_t num_threads = 1) const;

  // Scores, gradient, and Hessian of the log-likelihood (parametric,
  // continuous families only)

  //! @brief Bundles the per-observation scores.
  //!
  //! @details Mirrors `Vinecop::ScoresResult`; a single pair copula has no
  //! cascade caches to expose, so only the score matrix is carried.
  struct ScoresResult
  {
    Eigen::MatrixXd scores;
  };

  Eigen::MatrixXd scores(const Eigen::MatrixXd& u) const;

  Eigen::VectorXd gradient(const Eigen::MatrixXd& u) const;

  Eigen::MatrixXd hessian(const Eigen::MatrixXd& u) const;

  std::vector<Eigen::MatrixXd> hessian_full(const Eigen::MatrixXd& u) const;

  Eigen::MatrixXd scores_cov(const Eigen::MatrixXd& u) const;

  ScoresResult scores_full(const Eigen::MatrixXd& u) const;

  // Scores, gradient, and Hessian with per-row parameters
  Eigen::MatrixXd scores(const Eigen::MatrixXd& u,
                         const Eigen::MatrixXd& parameters,
                         const size_t num_threads = 1) const;

  Eigen::VectorXd gradient(const Eigen::MatrixXd& u,
                           const Eigen::MatrixXd& parameters,
                           const size_t num_threads = 1) const;

  Eigen::MatrixXd hessian(const Eigen::MatrixXd& u,
                          const Eigen::MatrixXd& parameters,
                          const size_t num_threads = 1) const;

  std::vector<Eigen::MatrixXd> hessian_full(const Eigen::MatrixXd& u,
                                            const Eigen::MatrixXd& parameters,
                                            const size_t num_threads = 1) const;

  Eigen::MatrixXd scores_cov(const Eigen::MatrixXd& u,
                             const Eigen::MatrixXd& parameters,
                             const size_t num_threads = 1) const;

  ScoresResult scores_full(const Eigen::MatrixXd& u,
                           const Eigen::MatrixXd& parameters,
                           const size_t num_threads = 1) const;

  Eigen::MatrixXd simulate(
    const size_t& n,
    const bool qrng = false,
    const std::vector<int>& seeds = std::vector<int>()) const;

  Eigen::MatrixXd simulate(const Eigen::MatrixXd& parameters,
                           const bool qrng = false,
                           const std::vector<int>& seeds = std::vector<int>(),
                           const size_t num_threads = 1) const;

  // Methods modifying the family/rotation/parameters
  void fit(const Eigen::MatrixXd& data,
           const FitControlsBicop& controls = FitControlsBicop());

  void select(const Eigen::MatrixXd& data,
              FitControlsBicop controls = FitControlsBicop());

  // Fit statistics
  double loglik(const Eigen::MatrixXd& u = Eigen::MatrixXd()) const;

  double aic(const Eigen::MatrixXd& u = Eigen::MatrixXd()) const;

  double bic(const Eigen::MatrixXd& u = Eigen::MatrixXd()) const;

  double mbic(const Eigen::MatrixXd& u = Eigen::MatrixXd(),
              const double psi0 = 0.9) const;

  // Misc
  std::string str() const;

  double parameters_to_tau(const Eigen::MatrixXd& parameters) const;

  Eigen::MatrixXd parameters_to_taildep(
    const Eigen::MatrixXd& parameters) const;

  double parameters_to_beta(const Eigen::MatrixXd& parameters) const;

  Eigen::MatrixXd tau_to_parameters(const double& tau) const;

  void flip();

  Eigen::MatrixXd get_parameters_lower_bounds() const;

  Eigen::MatrixXd get_parameters_upper_bounds() const;

  Bicop as_continuous() const;

private:
  // Evaluate the continuous copula represented by this object without
  // consulting or changing its stored variable types. These are used by the
  // inverse Rosenblatt transform, whose input and output always have
  // continuous uniform margins even when the fitted model has discrete
  // variables.
  Eigen::VectorXd hfunc1_continuous(const Eigen::MatrixXd& u,
                                    bool flipped) const;

  Eigen::VectorXd hinv2_continuous(const Eigen::MatrixXd& u,
                                   bool flipped) const;

  Eigen::MatrixXd format_data(const Eigen::MatrixXd& u) const;

  void rotate_data(Eigen::MatrixXd& u) const;

  Eigen::MatrixXd prep_for_abstract(const Eigen::MatrixXd& u) const;

  Eigen::MatrixXd prep_for_abstract_continuous(const Eigen::MatrixXd& u) const;

  struct ConditionalSpec
  {
    bool use_first;
    bool complement;
  };

  ConditionalSpec get_conditional_spec(bool first_function) const;

  static Eigen::VectorXd finalize_conditional(Eigen::VectorXd value,
                                              bool complement);

  Eigen::MatrixXd format_parameters(const Eigen::MatrixXd& u,
                                    const Eigen::MatrixXd& parameters) const;

  Eigen::VectorXd eval_in_batches(
    const Eigen::MatrixXd& u,
    const Eigen::MatrixXd& parameters_t,
    const size_t num_threads,
    const std::function<Eigen::VectorXd(const Eigen::MatrixXd&,
                                        const Eigen::MatrixXd&)>& f) const;

  // rotation-resolved derivative call: canonical selector for the (unrotated)
  // leaf, chain-rule sign, and, for h-functions, whether the rotation swaps
  // which h-function's leaf is used
  struct DerivSpec
  {
    std::string deriv;
    double sign;
    bool swap_hfunc;
  };

  size_t deriv_npars() const;

  void check_deriv_preconditions() const;

  DerivSpec map_pdf_deriv(const std::string& canonical) const;

  DerivSpec map_hfunc_deriv(const std::string& canonical,
                            bool first_hfunc) const;

  void check_rotation(int rotation) const;

  void check_data(const Eigen::MatrixXd& u) const;

  void check_data_dim(const Eigen::MatrixXd& u) const;

  void check_var_types(const std::vector<std::string>& var_types) const;

  void flip_abstract_var_types();

  void check_weights_size(const Eigen::VectorXd& weights,
                          const Eigen::MatrixXd& data) const;

  void check_fitted() const;

  unsigned short get_n_discrete() const;

  double compute_mbic_penalty(const size_t nobs, const double psi0) const;

  BicopPtr get_bicop() const;

  BicopPtr bicop_;
  int rotation_{ 0 };
  size_t nobs_{ 0 };
  mutable std::vector<std::string> var_types_;
};

//! @cond INTERNAL
//! A non-owning, optionally transposed view of a bivariate copula.
class BicopView
{
public:
  explicit BicopView(const Bicop& bicop,
                     bool flipped = false,
                     bool continuous = false);

  BicopView as_continuous() const;
  std::vector<std::string> get_var_types() const;
  Eigen::VectorXd hfunc1(const Eigen::MatrixXd& u) const;
  Eigen::VectorXd hfunc2(const Eigen::MatrixXd& u) const;
  Eigen::VectorXd hinv2(const Eigen::MatrixXd& u) const;

private:
  static Eigen::MatrixXd swap_arguments(const Eigen::MatrixXd& u);

  const Bicop* bicop_;
  bool flipped_;
  bool continuous_;
};
//! @endcond
}

#include <vinecopulib/bicop/implementation/class.ipp>
