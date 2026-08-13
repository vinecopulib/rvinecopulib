// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <Eigen/Dense>
#include <vector>
#include <vinecopulib/bicop/family.hpp>
#include <vinecopulib/misc/fit_controls.hpp>

namespace vinecopulib {
//! @brief A class for controlling fits of bivariate copula models.
//!
class FitControlsBicop
{
public:
  // Constructor
  FitControlsBicop(std::vector<BicopFamily> family_set = bicop_families::all,
                   std::string parametric_method = "mle",
                   std::string nonparametric_method = "constant",
                   double nonparametric_mult = 1.0,
                   size_t nonparametric_grid_size = 30,
                   std::string selection_criterion = "aic",
                   const Eigen::VectorXd& weights = Eigen::VectorXd(),
                   double psi0 = 0.9,
                   bool preselect_families = true,
                   bool allow_rotations = true,
                   size_t num_threads = 1);

  explicit FitControlsBicop(std::string parametric_method);

  explicit FitControlsBicop(std::string nonparametric_method,
                            double nonparametric_mult = 1.0,
                            size_t nonparametric_grid_size = 30);

  explicit FitControlsBicop(const FitControlsConfig& config);

  // Getters

  //! @return the set of pair-copula families considered during selection.
  std::vector<BicopFamily> get_family_set() const;

  //! @return the parametric estimator used for parametric families
  //! (one of `"mle"`, `"itau"`).
  std::string get_parametric_method() const;

  //! @return the nonparametric estimator used for the `tll` family
  //! (one of `"constant"`, `"linear"`, `"quadratic"`).
  std::string get_nonparametric_method() const;

  //! @return the bandwidth multiplier applied during nonparametric fits.
  double get_nonparametric_mult() const;

  //! @return the grid size used during nonparametric fits (number of grid
  //! points per axis).
  size_t get_nonparametric_grid_size() const;

  //! @return the model-selection criterion (one of `"loglik"`, `"aic"`,
  //! `"bic"`, `"mbic"`).
  std::string get_selection_criterion() const;

  //! @return the per-observation weights used during fitting (empty when
  //! observations are weighted equally).
  Eigen::VectorXd get_weights() const;

  //! @return whether families are pre-screened by a quick fit before the
  //! full criterion comparison.
  bool get_preselect_families() const;

  //! @return the prior probability of a non-independence copula used by
  //! the `mbic` criterion.
  double get_psi0() const;

  //! @return the number of threads used during fitting.
  size_t get_num_threads() const;

  //! @return whether non-zero copula rotations are considered during
  //! family selection.
  bool get_allow_rotations() const;

  // Setters
  void set_family_set(std::vector<BicopFamily> family_set);

  void set_parametric_method(std::string parametric_method);

  void set_nonparametric_method(std::string nonparametric_method);

  void set_nonparametric_mult(double nonparametric_mult);

  void set_nonparametric_grid_size(size_t nonparametric_grid_size);

  void set_selection_criterion(std::string selection_criterion);

  void set_weights(const Eigen::VectorXd& weights);

  void set_preselect_families(bool preselect_families);

  void set_psi0(double psi0);

  void set_num_threads(size_t num_threads);

  void set_allow_rotations(bool allow_rotations);

  // Misc
  std::string str() const;

protected:
  std::string str_internal(bool print_threads = true) const;

private:
  std::vector<BicopFamily> family_set_;
  std::string parametric_method_;
  std::string nonparametric_method_;
  double nonparametric_mult_;
  size_t nonparametric_grid_size_;
  std::string selection_criterion_;
  Eigen::VectorXd weights_;
  bool preselect_families_;
  double psi0_;
  size_t num_threads_;
  bool allow_rotations_;

  void check_parametric_method(std::string parametric_method);

  void check_nonparametric_method(std::string nonparametric_method);

  void check_nonparametric_mult(double nonparametric_mult);

  void check_nonparametric_grid_size(size_t nonparametric_grid_size);

  void check_selection_criterion(std::string selection_criterion);

  void check_psi0(double psi0);

  size_t process_num_threads(size_t num_threads);
};
}

#include <vinecopulib/bicop/implementation/fit_controls.ipp>
