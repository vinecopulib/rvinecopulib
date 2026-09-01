// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <boost/random/mersenne_twister.hpp>
#include <limits>
#include <vinecopulib/bicop/fit_controls.hpp>

namespace vinecopulib {
//! @brief A class for controlling fits of vine copula models.
//!
class FitControlsVinecop : public FitControlsBicop
{
public:
  // Constructor
  FitControlsVinecop();

  explicit FitControlsVinecop(
    std::vector<BicopFamily> family_set,
    std::string parametric_method = "mle",
    std::string nonparametric_method = "constant",
    double nonparametric_mult = 1.0,
    size_t nonparametric_grid_size = 30,
    size_t trunc_lvl = std::numeric_limits<size_t>::max(),
    std::string tree_criterion = "tau",
    double threshold = 0.0,
    std::string selection_criterion = "aic",
    const Eigen::VectorXd& weights = Eigen::VectorXd(),
    double psi0 = 0.9,
    bool preselect_families = true,
    bool select_trunc_lvl = false,
    bool select_threshold = false,
    bool select_families = true,
    bool show_trace = false,
    size_t num_threads = 1,
    std::string tree_algorithm = "mst_prim",
    bool allow_rotations = true,
    std::vector<int> seeds = std::vector<int>());

  explicit FitControlsVinecop(
    const FitControlsBicop& controls,
    size_t trunc_lvl = std::numeric_limits<size_t>::max(),
    std::string tree_criterion = "tau",
    double threshold = 0.0,
    bool select_trunc_lvl = false,
    bool select_threshold = false,
    bool select_families = true,
    bool show_trace = false,
    std::string tree_algorithm = "mst_prim",
    std::vector<int> seeds = std::vector<int>());

  explicit FitControlsVinecop(const FitControlsConfig& config);

  // Getters
  //! @return the truncation level (the number of trees that will be fit;
  //! pair copulas above this level are forced to independence).
  size_t get_trunc_lvl() const;

  //! @return the edge-weighting criterion used to grow the structure
  //! (one of `"tau"`, `"rho"`, `"hoeffd"`, `"mcor"`, `"cxi"`, `"joe"`,
  //! `"custom"`).
  std::string get_tree_criterion() const;

  //! @return the custom edge-weight function used when `tree_criterion` is
  //! `"custom"` (empty otherwise). It is always called on the thread that
  //! starts the fit, so it need not be thread safe.
  TreeCriterionFunction get_tree_criterion_function() const;

  //! @return the absolute-dependence threshold below which pair copulas
  //! are set to independence during structure selection (`0` disables).
  double get_threshold() const;

  //! @return whether progress information is printed during fitting.
  bool get_show_trace() const;

  //! @return whether the truncation level is selected automatically via
  //! the mBICv criterion during fitting.
  bool get_select_trunc_lvl() const;

  //! @return whether the threshold is selected automatically during
  //! fitting.
  bool get_select_threshold() const;

  //! @return whether pair-copula families are selected during fitting
  //! (when `false`, pre-specified families are used).
  bool get_select_families() const;

  //! @return whether sparse selection (truncation or thresholding) is
  //! enabled in the current configuration.
  bool needs_sparse_select() const;

  //! @return the bicop fit controls used for each pair-copula fit.
  FitControlsBicop get_fit_controls_bicop() const;

  //! @return the structure-selection algorithm (one of `"mst_prim"` for
  //! Dissmann's greedy heuristic or `"random_weighted"` for Wilson-weighted
  //! random spanning trees).
  std::string get_tree_algorithm() const;

  //! @return the random seeds used by the structure-selection RNG (empty
  //! to use a non-reproducible seed).
  std::vector<int> get_seeds() const;

  //! @return the conditioning set (1-based variable indices) placed at the
  //! tail of the vine order during structure selection; empty for
  //! unconditional selection.
  std::vector<size_t> get_conditioning_set() const;

  boost::random::mt19937 get_rng() const;

  // Setters
  void set_trunc_lvl(size_t trunc_lvl);

  void set_tree_criterion(std::string tree_criterion);

  //! Sets the custom edge-weight function used when `tree_criterion` is
  //! `"custom"`. The two may be set in either order, but a fit is rejected
  //! unless both are set. The function is always called on the thread that
  //! starts the fit, so it need not be thread safe; the pair-copula fits still
  //! use `num_threads` threads.
  //! @param tree_criterion_function a callable mapping a two-column matrix of
  //! pair-copula data and a vector of weights to a scalar dependence value.
  void set_tree_criterion_function(
    TreeCriterionFunction tree_criterion_function);

  void set_threshold(double threshold);

  void set_show_trace(bool show_trace);

  void set_select_trunc_lvl(bool select_trunc_lvl);

  void set_select_threshold(bool select_threshold);

  void set_select_families(bool select_families);

  void set_fit_controls_bicop(const FitControlsBicop& controls);

  void set_tree_algorithm(std::string tree_algorithm);

  void set_seeds(std::vector<int> seeds);

  //! Sets the conditioning set for conditioning-aware structure selection.
  //! @param conditioning_set 1-based variable indices to place at the tail of
  //! the vine order (so they can be conditioned on with
  //! `Vinecop::simulate_conditional()`); empty for unconditional selection.
  void set_conditioning_set(std::vector<size_t> conditioning_set);

  // Misc
  std::string str() const;

private:
  // In-class defaults so that construction from a partial `FitControlsConfig`
  // (which only assigns the fields present in the config) leaves every member
  // in a well-defined state.
  size_t trunc_lvl_ = std::numeric_limits<size_t>::max();
  std::string tree_criterion_ = "tau";
  TreeCriterionFunction tree_criterion_function_ = {};
  double threshold_ = 0.0;
  bool show_trace_ = false;
  bool select_trunc_lvl_ = false;
  bool select_threshold_ = false;
  bool select_families_ = true;
  std::string tree_algorithm_ = "mst_prim";
  std::vector<int> seeds_;
  std::vector<size_t> conditioning_set_ = {};
  boost::random::mt19937 rng_;

  void check_tree_criterion(std::string tree_criterion);

  void check_threshold(double threshold);

  void check_conditioning_set(const std::vector<size_t>& conditioning_set);
};
}

#include <vinecopulib/vinecop/implementation/fit_controls.ipp>
