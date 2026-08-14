// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <Eigen/Dense>
#include <utility>
#include <vinecopulib/vinecop/fit_controls.hpp>
#include <vinecopulib/vinecop/rvine_structure.hpp>

namespace vinecopulib {

// forward declarations
class Bicop;
namespace tools_select {
class VinecopSelector;
}

//! @brief A class for vine copula models.
//!
//! @details A vine copula model is characterized by its structure (see
//! `RVineStructure` objects) and the pair-copulas (see `Bicop` objects).
class Vinecop
{
public:
  // default constructors
  Vinecop() = default;

  explicit Vinecop(size_t d);

  // Constructors without data
  explicit Vinecop(const RVineStructure& structure,
                   const std::vector<std::vector<Bicop>>& pair_copulas = {},
                   const std::vector<std::string>& var_types = {});

  explicit Vinecop(
    const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& matrix,
    const std::vector<std::vector<Bicop>>& pair_copulas = {},
    const std::vector<std::string>& var_types = {});

  // Constructors from data
  explicit Vinecop(const Eigen::MatrixXd& data,
                   const RVineStructure& structure = RVineStructure(),
                   const std::vector<std::string>& var_types = {},
                   const FitControlsVinecop& controls = FitControlsVinecop());

  // `matrix` must not be defaulted: that makes `Vinecop(data)` ambiguous with
  // the overload above.
  explicit Vinecop(
    const Eigen::MatrixXd& data,
    const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& matrix,
    const std::vector<std::string>& var_types = {},
    const FitControlsVinecop& controls = FitControlsVinecop());

  // Constructors from files/serialized objects
  explicit Vinecop(const std::string& filename, const bool check = true);
  explicit Vinecop(const nlohmann::json& input, const bool check = true);

  // Serialize
  nlohmann::json to_json() const;
  void to_file(const std::string& filename) const;

  // Methods modifying structure and/or families and parameters
  void select(const Eigen::MatrixXd& data,
              const FitControlsVinecop& controls = FitControlsVinecop());

  void fit(const Eigen::MatrixXd& data,
           const FitControlsBicop& controls = FitControlsBicop(),
           const size_t num_threads = 1);

  // Getters for a single pair copula

  //! @return the pair copula at the given `(tree, edge)` position.
  Bicop get_pair_copula(size_t tree, size_t edge) const;

  //! @return the family of the pair copula at the given `(tree, edge)`.
  BicopFamily get_family(size_t tree, size_t edge) const;

  //! @return the rotation of the pair copula at the given `(tree, edge)`,
  //! in degrees (one of `0`, `90`, `180`, `270`).
  int get_rotation(size_t tree, size_t edge) const;

  //! @return the parameter matrix of the pair copula at the given
  //! `(tree, edge)`.
  Eigen::MatrixXd get_parameters(size_t tree, size_t edge) const;

  //! @return Kendall's tau of the pair copula at the given `(tree, edge)`.
  double get_tau(size_t tree, size_t edge) const;

  //! @return the truncation level (number of fitted trees; pair copulas in
  //! trees above this level are forced to independence).
  size_t get_trunc_lvl() const;

  // Getters for all pair copulas

  //! @return all pair copulas, indexed as `[tree][edge]`.
  std::vector<std::vector<Bicop>> get_all_pair_copulas() const;

  //! @return the families of all pair copulas, indexed as `[tree][edge]`.
  std::vector<std::vector<BicopFamily>> get_all_families() const;

  //! @return the rotations of all pair copulas, indexed as `[tree][edge]`.
  std::vector<std::vector<int>> get_all_rotations() const;

  //! @return the parameter matrices of all pair copulas, indexed as
  //! `[tree][edge]`.
  std::vector<std::vector<Eigen::MatrixXd>> get_all_parameters() const;

  //! @return Kendall's tau for every pair copula, indexed as `[tree][edge]`.
  std::vector<std::vector<double>> get_all_taus() const;

  // Getters for the structure

  //! @return the dimension of the vine (number of variables).
  size_t get_dim() const;

  //! @return the natural variable order in the first tree.
  std::vector<size_t> get_order() const;

  //! @return the underlying R-vine structure.
  RVineStructure get_rvine_structure() const;

  //! @return the R-vine structure as a square matrix
  //! (cf. Joe, 2014, Section 6.13).
  Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> get_matrix() const;

  //! @return the R-vine structure as a triangular array.
  //! @param natural_order whether to return the array in natural order or
  //! original variable order.
  TriangularArray<size_t> get_struct_array(bool natural_order = false) const;

  //! @return the vine as a list-of-trees decomposition, with each edge
  //! carrying its fitted pair-copula.
  RVineTrees get_trees() const;

  // getters for fit statistics

  //! @return the absolute-dependence threshold used during structure
  //! selection (`0` if thresholding was disabled).
  double get_threshold() const;
  //! @return the log-likelihood of the data the vine was fitted to.
  double get_loglik() const;
  //! @return the number of observations used to fit the vine (`0` until
  //! `fit()` has been called).
  size_t get_nobs() const;
  //! @return Akaike's information criterion of the fitted vine.
  double get_aic() const;
  //! @return the Bayesian information criterion of the fitted vine.
  double get_bic() const;
  //! @return the modified Bayesian information criterion for vines at the
  //! fitted parameters; `psi0` is the prior probability of a
  //! non-independence pair copula.
  double get_mbicv(const double psi0 = 0.9) const;

  // Stats methods
  Eigen::VectorXd pdf(Eigen::MatrixXd u, const size_t num_threads = 1) const;

  //! @brief The density together with the per-edge quantities computed on the
  //! way, as returned by `pdf_full()`.
  //!
  //! The triangular arrays are indexed `(tree, edge)`. `_sub` holds the
  //! h-functions of the second ("sub") argument needed for discrete variables.
  //! The derivative routines take these as input rather than recomputing the
  //! h-function cascade.
  struct PdfWithHfuncsResult
  {
    Eigen::VectorXd pdf;
    TriangularArray<Eigen::VectorXd> pdf_edges;
    TriangularArray<Eigen::VectorXd> hfunc1;
    TriangularArray<Eigen::VectorXd> hfunc2;
    TriangularArray<Eigen::VectorXd> hfunc1_sub;
    TriangularArray<Eigen::VectorXd> hfunc2_sub;
  };

  PdfWithHfuncsResult pdf_full(Eigen::MatrixXd u,
                               const size_t num_threads = 1,
                               const bool keep_all = true) const;

  // Stats methods with per-observation parameters. `parameters` is an
  // n x npars matrix, one full-vine parameter vector per observation, with
  // columns in the (tree, edge, parameter) order of scores(). Continuous,
  // all-parametric models only.
  Eigen::VectorXd pdf(Eigen::MatrixXd u,
                      const Eigen::MatrixXd& parameters,
                      const size_t num_threads = 1) const;

  PdfWithHfuncsResult pdf_full(Eigen::MatrixXd u,
                               const Eigen::MatrixXd& parameters,
                               const size_t num_threads = 1,
                               const bool keep_all = true) const;

  Eigen::VectorXd cdf(const Eigen::MatrixXd& u,
                      const size_t N = 10000,
                      const size_t num_threads = 1,
                      const std::vector<int>& seeds = std::vector<int>()) const;

  Eigen::MatrixXd simulate(
    const size_t n,
    const bool qrng = false,
    const size_t num_threads = 1,
    const std::vector<int>& seeds = std::vector<int>()) const;

  Eigen::MatrixXd simulate_conditional(
    const Eigen::MatrixXd& u_cond,
    const bool qrng = false,
    const size_t num_threads = 1,
    const std::vector<int>& seeds = std::vector<int>()) const;
  Eigen::MatrixXd simulate_conditional(
    const Eigen::MatrixXd& u_cond,
    const std::vector<size_t>& conditioning_set,
    const bool qrng = false,
    const size_t num_threads = 1,
    const std::vector<int>& seeds = std::vector<int>()) const;

  void reorient(const std::vector<size_t>& conditioning_set);

  Eigen::MatrixXd rosenblatt(Eigen::MatrixXd u,
                             const size_t num_threads = 1,
                             bool randomize_discrete = true,
                             std::vector<int> seeds = {}) const;
  Eigen::MatrixXd rosenblatt(Eigen::MatrixXd u,
                             const std::vector<size_t>& conditioning_set,
                             const size_t num_threads = 1,
                             bool randomize_discrete = true,
                             std::vector<int> seeds = {}) const;
  Eigen::MatrixXd inverse_rosenblatt(const Eigen::MatrixXd& u,
                                     const size_t num_threads = 1) const;
  Eigen::MatrixXd inverse_rosenblatt(
    const Eigen::MatrixXd& u,
    const std::vector<size_t>& conditioning_set,
    const size_t num_threads = 1) const;

  //! Sets every pair copula in one shot.
  //! @param pair_copulas nested list of `Bicop` instances, shaped like
  //! `[tree][edge]` with `dim - 1 - tree` edges in tree `tree`.
  void set_all_pair_copulas(
    const std::vector<std::vector<Bicop>>& pair_copulas);
  //! Sets the variable types.
  //! @param var_types a length-`dim` vector with each entry either `"c"`
  //! (continuous) or `"d"` (discrete).
  void set_var_types(const std::vector<std::string>& var_types);

  //! @return the variable types of each variable (each `"c"` for
  //! continuous or `"d"` for discrete).
  std::vector<std::string> get_var_types() const;

  // Fit statistics
  //! @return the total number of parameters across all pair copulas.
  //! For nonparametric families, this is a conceptually similar
  //! effective-parameter count.
  double get_npars() const;

  double loglik(const Eigen::MatrixXd& u = Eigen::MatrixXd(),
                const size_t num_threads = 1) const;

  //! Log-likelihood with per-observation parameters (see the per-observation
  //! `pdf()` overload for the `parameters` layout and restrictions).
  double loglik(const Eigen::MatrixXd& u,
                const Eigen::MatrixXd& parameters,
                const size_t num_threads = 1) const;

  double aic(const Eigen::MatrixXd& u = Eigen::MatrixXd(),
             const size_t num_threads = 1) const;

  double bic(const Eigen::MatrixXd& u = Eigen::MatrixXd(),
             const size_t num_threads = 1) const;

  double mbicv(const Eigen::MatrixXd& u = Eigen::MatrixXd(),
               const double psi0 = 0.9,
               const size_t num_threads = 1) const;

  // Misc methods
  static std::vector<std::vector<Bicop>> make_pair_copula_store(
    const size_t d,
    const size_t trunc_lvl = std::numeric_limits<size_t>::max());
  void truncate(size_t trunc_lvl);

  std::string str(const std::vector<size_t>& trees = {}) const;

  //! @brief The scores together with the per-edge quantities computed on the
  //! way, as returned by `scores_full()`.
  //!
  //! `scores` is \f$ n \times p \f$ with columns in `(tree, edge, parameter)`
  //! order. The remaining members are indexed `(tree, edge)` and hold the
  //! pair-copula densities and the derivatives of the log-density and
  //! h-functions with respect to the parameters and the arguments.
  struct ScoresResult
  {
    Eigen::MatrixXd scores;
    TriangularArray<Eigen::VectorXd> pdf_edges;
    TriangularArray<std::vector<Eigen::VectorXd>> logpdf_deriv_pars;
    TriangularArray<std::vector<Eigen::VectorXd>> hfunc1_deriv_pars;
    TriangularArray<std::vector<Eigen::VectorXd>> hfunc2_deriv_pars;
    TriangularArray<Eigen::VectorXd> logpdf_deriv_u1;
    TriangularArray<Eigen::VectorXd> logpdf_deriv_u2;
    TriangularArray<Eigen::VectorXd> hfunc1_deriv_u1;
    TriangularArray<Eigen::VectorXd> hfunc2_deriv_u2;
  };

  ScoresResult scores_full(Eigen::MatrixXd u,
                           bool step_wise = true,
                           const size_t num_threads = 1,
                           const bool keep_all = true);
  Eigen::MatrixXd scores(Eigen::MatrixXd u,
                         bool step_wise = true,
                         const size_t num_threads = 1);
  Eigen::VectorXd gradient(Eigen::MatrixXd u,
                           bool step_wise = true,
                           const size_t num_threads = 1);
  Eigen::MatrixXd hessian(Eigen::MatrixXd u,
                          bool step_wise = true,
                          const size_t num_threads = 1);
  TriangularArray<std::vector<Eigen::MatrixXd>> hessian_full(
    Eigen::MatrixXd u,
    bool step_wise = true,
    const size_t num_threads = 1);
  Eigen::MatrixXd scores_cov(Eigen::MatrixXd u,
                             bool step_wise = true,
                             const size_t num_threads = 1);

  // Scores, gradient, and Hessian with per-observation parameters. `parameters`
  // is an n x npars matrix, one full-vine parameter vector per observation,
  // with columns in the (tree, edge, parameter) order of scores(). Continuous,
  // all-parametric models only.
  ScoresResult scores_full(Eigen::MatrixXd u,
                           const Eigen::MatrixXd& parameters,
                           bool step_wise = true,
                           const size_t num_threads = 1,
                           const bool keep_all = true);
  Eigen::MatrixXd scores(Eigen::MatrixXd u,
                         const Eigen::MatrixXd& parameters,
                         bool step_wise = true,
                         const size_t num_threads = 1);
  Eigen::VectorXd gradient(Eigen::MatrixXd u,
                           const Eigen::MatrixXd& parameters,
                           bool step_wise = true,
                           const size_t num_threads = 1);
  Eigen::MatrixXd hessian(Eigen::MatrixXd u,
                          const Eigen::MatrixXd& parameters,
                          bool step_wise = true,
                          const size_t num_threads = 1);
  TriangularArray<std::vector<Eigen::MatrixXd>> hessian_full(
    Eigen::MatrixXd u,
    const Eigen::MatrixXd& parameters,
    bool step_wise = true,
    const size_t num_threads = 1);
  Eigen::MatrixXd scores_cov(Eigen::MatrixXd u,
                             const Eigen::MatrixXd& parameters,
                             bool step_wise = true,
                             const size_t num_threads = 1);

private:
  struct ReorientationMap
  {
    RVineStructure structure;
    TriangularArray<RVineTrees::PairCopulaLocation> pair_copulas;
    bool identity{ false };
  };

  class VinecopView
  {
  public:
    explicit VinecopView(const Vinecop& vinecop);
    VinecopView(const Vinecop& vinecop, const ReorientationMap& reorientation);

    const RVineStructure& get_structure() const;
    BicopView get_pair_copula(size_t tree, size_t edge) const;

  private:
    const Vinecop* vinecop_;
    const ReorientationMap* reorientation_;
  };

  ReorientationMap make_reorientation_map(
    const std::vector<size_t>& conditioning_set) const;
  Eigen::MatrixXd simulate_conditional_impl(
    const Eigen::MatrixXd& u_cond,
    const std::vector<size_t>& conditioning_set,
    const VinecopView& view,
    bool qrng,
    size_t num_threads,
    const std::vector<int>& seeds) const;
  Eigen::MatrixXd rosenblatt_impl(Eigen::MatrixXd u,
                                  const VinecopView& view,
                                  size_t num_threads,
                                  bool randomize_discrete,
                                  std::vector<int> seeds) const;
  Eigen::MatrixXd inverse_rosenblatt_impl(const Eigen::MatrixXd& u,
                                          const VinecopView& view,
                                          size_t num_threads) const;

  // Per-edge derivative caches shared by the analytic score/gradient/Hessian
  // cascades. One forward walk over the vine (build_deriv_cache) fills them;
  // the cascades then only read them.
  //
  // Each edge's pair copula c(u1, u2) consumes two arguments produced by the
  // previous tree and produces (up to) two h-function outputs consumed by
  // deeper trees: hfunc1 = P(u2 <= . | u1) and hfunc2 = P(u1 <= . | u2).
  // A DerivLeaf holds the derivatives of ONE such h-function output; it is
  // what the chain rule propagates through when a parameter perturbation
  // travels from its own edge to the deeper edges that (indirectly) depend
  // on it. For an h-function output, `du1`/`du2` are ∂h/∂u1, ∂h/∂u2 (one of
  // them equals the copula density `c` by the identity ∂h2/∂u1 = ∂h1/∂u2 =
  // c).
  struct DerivLeaf
  {
    Eigen::VectorXd du1, du2;          // ∂h/∂u1, ∂h/∂u2
    std::vector<Eigen::VectorXd> dpar; // ∂h/∂θ_p (cascade seed)
    // second-order (only when requested):
    Eigen::VectorXd du1u1, du1u2, du2u2;           // ∂²h/∂{u1²,u1u2,u2²}
    std::vector<Eigen::VectorXd> dpar_u1, dpar_u2; // ∂²h/∂θ_p∂{u1,u2}
    std::vector<std::vector<Eigen::VectorXd>> dpar_par; // ∂²h/∂θ_p∂θ_q
    // whether any deeper tree consumes this h-function output (from the
    // structure's needed_hfunc1/2 masks); inactive leaves are left empty
    // and the cascades skip them
    bool active{ false };
  };
  // All derivative data of one edge: the log-density derivatives (du*/dpar*,
  // used for the score/Hessian contributions of the edge itself) plus the
  // two output leaves `h1` (its hfunc1 output) and `h2` (its hfunc2 output)
  // through which perturbations propagate to deeper trees. `arg2_col` is this
  // edge's min_array entry: the storage column of the previous tree that
  // provides the edge's second argument; `arg2_is_h2` says whether that
  // argument is that column's hfunc2 (true) or hfunc1 (false) value —
  // mirroring how the pdf/rosenblatt passes assemble their arguments.
  // (The du*/dpar* here are derivatives of `log c`; the identically named
  // fields of `DerivLeaf` are derivatives of an h-function `h`.)
  struct DerivCache
  {
    size_t np{ 0 }, arg2_col{ 0 };
    bool arg2_is_h2{ false };
    Eigen::VectorXd c, du1, du2;       // c, ∂logc/∂u1, ∂logc/∂u2
    std::vector<Eigen::VectorXd> dpar; // ∂logc/∂θ_p (step-wise score)
    // second-order log-density (only when requested):
    Eigen::VectorXd du1u1, du1u2, du2u2;
    std::vector<Eigen::VectorXd> dpar_u1, dpar_u2;
    std::vector<std::vector<Eigen::VectorXd>> dpar_par;
    DerivLeaf h1, h2; // the edge's hfunc1 and hfunc2 output leaves
  };
  // one forward walk over the rows [begin, begin + size) of `u`, filling the
  // per-edge derivative caches (second-order fields only when `second_order`).
  // When `per_obs_params` is non-empty (an n x npars matrix), each edge reads
  // its own per-observation parameters from the matching column block instead
  // of the pair copula's stored parameters; an empty matrix (the default)
  // leaves the fixed-parameter fast path unchanged.
  TriangularArray<DerivCache> build_deriv_cache(
    const Eigen::MatrixXd& u,
    size_t begin,
    size_t size,
    bool second_order,
    const Eigen::MatrixXd& per_obs_params = Eigen::MatrixXd()) const;

  // throws if any pair copula is nonparametric (differentiating w.r.t. an
  // interpolation grid is meaningless); `fn` names the calling method.
  void check_parametric(const char* fn) const;

  // validates a per-observation parameter matrix for the score/Hessian
  // overloads: rejects discrete variables and checks the n x npars shape.
  void check_per_obs_params(const Eigen::MatrixXd& u,
                            const Eigen::MatrixXd& per_obs_params) const;

protected:
  size_t d_{ 1 };
  RVineStructure rvine_structure_;
  mutable std::vector<std::vector<Bicop>> pair_copulas_;
  double threshold_{ 0.0 };
  double loglik_{ NAN };
  size_t nobs_{ 0 };
  std::vector<std::string> var_types_;
  int n_discrete_{ 0 };

  void check_data_dim(const Eigen::MatrixXd& data) const;
  void check_data(const Eigen::MatrixXd& data) const;
  void check_pair_copulas_rvine_structure(
    const std::vector<std::vector<Bicop>>& pair_copulas) const;
  double calculate_mbicv_penalty(const size_t nobs, const double psi0) const;
  void finalize_fit(const tools_select::VinecopSelector& selector);
  void check_conditioning_set(const std::vector<size_t>& conditioning_set,
                              const FitControlsVinecop& controls) const;
  static void check_tree_criterion_function(const FitControlsVinecop& controls);
  void check_weights_size(const Eigen::VectorXd& weights,
                          const Eigen::MatrixXd& data) const;
  void check_enough_data(const Eigen::MatrixXd& data) const;
  void check_fitted() const;
  void check_indices(const size_t tree, const size_t edge) const;
  void check_var_types(const std::vector<std::string>& var_types) const;
  void set_continuous_var_types();
  void set_var_types_internal(const std::vector<std::string>& var_types);
  int get_n_discrete() const;
  bool is_discrete() const;
  Eigen::MatrixXd collapse_data(const Eigen::MatrixXd& u) const;
  void collapse_data_inplace(Eigen::MatrixXd& u) const;
};
}

#include <vinecopulib/vinecop/implementation/class.ipp>
