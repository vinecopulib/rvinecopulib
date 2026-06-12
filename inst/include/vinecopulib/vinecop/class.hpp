// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
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
  Vinecop() {}

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

  explicit Vinecop(
    const Eigen::MatrixXd& data,
    const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& matrix =
      Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>(),
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

  DEPRECATED void select_all(
    const Eigen::MatrixXd& data,
    const FitControlsVinecop& controls = FitControlsVinecop());

  DEPRECATED void select_families(
    const Eigen::MatrixXd& data,
    const FitControlsVinecop& controls = FitControlsVinecop());

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
                               const size_t num_threads,
                               const bool keep_all = true) const;

  Eigen::VectorXd cdf(const Eigen::MatrixXd& u,
                      const size_t N = 1e4,
                      const size_t num_threads = 1,
                      std::vector<int> seeds = std::vector<int>()) const;

  Eigen::MatrixXd simulate(
    const size_t n,
    const bool qrng = false,
    const size_t num_threads = 1,
    const std::vector<int>& seeds = std::vector<int>()) const;

  Eigen::MatrixXd rosenblatt(Eigen::MatrixXd u,
                             const size_t num_threads = 1,
                             bool randomize_discrete = true,
                             std::vector<int> seeds = {}) const;
  Eigen::MatrixXd inverse_rosenblatt(const Eigen::MatrixXd& u,
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
  Eigen::MatrixXd scores(Eigen::MatrixXd u,
                         bool step_wise = true,
                         const size_t num_threads = 1);
  TriangularArray<std::vector<Eigen::MatrixXd>> hessian(
    Eigen::MatrixXd u,
    bool step_wise = true,
    const size_t num_threads = 1);
  Eigen::MatrixXd hessian_avg(Eigen::MatrixXd u,
                              bool step_wise = true,
                              const size_t num_threads = 1);
  Eigen::MatrixXd scores_cov(Eigen::MatrixXd u,
                             bool step_wise = true,
                             const size_t num_threads = 1);

protected:
  size_t d_{ 1 };
  RVineStructure rvine_structure_;
  mutable std::vector<std::vector<Bicop>> pair_copulas_;
  double threshold_{ 0.0 };
  double loglik_{ NAN };
  size_t nobs_{ 0 };
  mutable std::vector<std::string> var_types_;

  void check_data_dim(const Eigen::MatrixXd& data) const;
  void check_data(const Eigen::MatrixXd& data) const;
  void check_pair_copulas_rvine_structure(
    const std::vector<std::vector<Bicop>>& pair_copulas) const;
  double calculate_mbicv_penalty(const size_t nobs, const double psi0) const;
  void finalize_fit(const tools_select::VinecopSelector& selector);
  void check_weights_size(const Eigen::VectorXd& weights,
                          const Eigen::MatrixXd& data) const;
  void check_enough_data(const Eigen::MatrixXd& data) const;
  void check_fitted() const;
  void check_indices(const size_t tree, const size_t edge) const;
  void check_var_types(const std::vector<std::string>& var_types) const;
  void set_continuous_var_types() const;
  void set_var_types_internal(const std::vector<std::string>& var_types) const;
  int get_n_discrete() const;
  bool is_discrete() const;
  Eigen::MatrixXd collapse_data(const Eigen::MatrixXd& u) const;
};
}

#include <vinecopulib/vinecop/implementation/class.ipp>
