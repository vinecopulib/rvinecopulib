// Copyright © 2016-2021 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <vinecopulib/bicop/fit_controls.hpp>
#include <vinecopulib/misc/nlohmann_json.hpp>

namespace vinecopulib {

// forward declaration of Abstract class
class AbstractBicop;
using BicopPtr = std::shared_ptr<AbstractBicop>;

//! @brief A class for bivariate copula models.
//!
//! The copula model is fully characterized by the family, rotation,
//! and parameters.
class Bicop
{

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

  explicit Bicop(const std::string& filename);

  explicit Bicop(const nlohmann::json& input);

  Bicop& operator=(Bicop other);

  // Serialize
  nlohmann::json to_json() const;

  void to_file(const std::string& filename) const;

  // Getters and setters
  BicopFamily get_family() const;

  std::string get_family_name() const;

  int get_rotation() const;

  Eigen::MatrixXd get_parameters() const;

  double get_tau() const;

  double get_npars() const;

  double get_loglik() const;
  size_t get_nobs() const;
  double get_aic() const;
  double get_bic() const;
  double get_mbic(const double psi0 = 0.9) const;

  void set_rotation(const int rotation);

  void set_parameters(const Eigen::MatrixXd& parameters);

  void set_var_types(const std::vector<std::string>& var_types = { "c", "c" });

  std::vector<std::string> get_var_types() const;

  // Stats methods
  Eigen::VectorXd pdf(const Eigen::MatrixXd& u) const;

  Eigen::VectorXd cdf(const Eigen::MatrixXd& u) const;

  Eigen::VectorXd hfunc1(const Eigen::MatrixXd& u) const;

  Eigen::VectorXd hfunc2(const Eigen::MatrixXd& u) const;

  Eigen::VectorXd hinv1(const Eigen::MatrixXd& u) const;

  Eigen::VectorXd hinv2(const Eigen::MatrixXd& u) const;

  Eigen::MatrixXd simulate(
    const size_t& n,
    const bool qrng = false,
    const std::vector<int>& seeds = std::vector<int>()) const;

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

  Eigen::MatrixXd tau_to_parameters(const double& tau) const;

  void flip();

  Eigen::MatrixXd get_parameters_lower_bounds() const;

  Eigen::MatrixXd get_parameters_upper_bounds() const;

  Bicop as_continuous() const;

private:
  Eigen::MatrixXd format_data(const Eigen::MatrixXd& u) const;

  void rotate_data(Eigen::MatrixXd& u) const;

  Eigen::MatrixXd prep_for_abstract(const Eigen::MatrixXd& u) const;

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
}

#include <vinecopulib/bicop/implementation/class.ipp>
