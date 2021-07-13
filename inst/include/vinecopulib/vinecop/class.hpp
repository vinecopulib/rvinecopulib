// Copyright © 2016-2021 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <Eigen/Dense>
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
//! A vine copula model is characterized by its structure (see
//! `RVineStructure` objects) and the pair-copulas.
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

  DEPRECATED void select_all(
    const Eigen::MatrixXd& data,
    const FitControlsVinecop& controls = FitControlsVinecop());

  DEPRECATED void select_families(
    const Eigen::MatrixXd& data,
    const FitControlsVinecop& controls = FitControlsVinecop());

  // Getters for a single pair copula
  Bicop get_pair_copula(size_t tree, size_t edge) const;

  BicopFamily get_family(size_t tree, size_t edge) const;

  int get_rotation(size_t tree, size_t edge) const;

  Eigen::MatrixXd get_parameters(size_t tree, size_t edge) const;

  double get_tau(size_t tree, size_t edge) const;

  size_t get_trunc_lvl() const;

  // Getters for all pair copulas
  std::vector<std::vector<Bicop>> get_all_pair_copulas() const;

  std::vector<std::vector<BicopFamily>> get_all_families() const;

  std::vector<std::vector<int>> get_all_rotations() const;

  std::vector<std::vector<Eigen::MatrixXd>> get_all_parameters() const;

  std::vector<std::vector<double>> get_all_taus() const;

  // Getters for the structure
  size_t get_dim() const;

  std::vector<size_t> get_order() const;

  RVineStructure get_rvine_structure() const;

  Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> get_matrix() const;

  TriangularArray<size_t> get_struct_array(bool natural_order = false) const;

  // getters for fit statistics
  double get_threshold() const;
  double get_loglik() const;
  size_t get_nobs() const;
  double get_aic() const;
  double get_bic() const;
  double get_mbicv(const double psi0 = 0.9) const;

  // Stats methods
  Eigen::VectorXd pdf(Eigen::MatrixXd u, const size_t num_threads = 1) const;

  Eigen::VectorXd cdf(const Eigen::MatrixXd& u,
                      const size_t N = 1e4,
                      const size_t num_threads = 1,
                      std::vector<int> seeds = std::vector<int>()) const;

  Eigen::MatrixXd simulate(
    const size_t n,
    const bool qrng = false,
    const size_t num_threads = 1,
    const std::vector<int>& seeds = std::vector<int>()) const;

  Eigen::MatrixXd rosenblatt(const Eigen::MatrixXd& u,
                             const size_t num_threads = 1) const;
  Eigen::MatrixXd inverse_rosenblatt(const Eigen::MatrixXd& u,
                                     const size_t num_threads = 1) const;

  void set_all_pair_copulas(
    const std::vector<std::vector<Bicop>>& pair_copulas);
  void set_var_types(const std::vector<std::string>& var_types);

  std::vector<std::string> get_var_types() const;

  // Fit statistics
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

  std::string str() const;

protected:
  size_t d_;
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
  Eigen::MatrixXd collapse_data(const Eigen::MatrixXd& u) const;
};
}

#include <vinecopulib/vinecop/implementation/class.ipp>
