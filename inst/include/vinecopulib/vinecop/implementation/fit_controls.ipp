// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <boost/random/seed_seq.hpp>
#include <random>
#include <stdexcept>
#include <utility>
#include <vinecopulib/misc/tools_stl.hpp>

//! @file vinecop/implementation/fit_controls.ipp
//! @brief Fit controls for Vinecop class (Implementation).

//! Tools for bivariate and vine copula modeling
namespace vinecopulib {
//! Instantiates default controls for fitting vine copula models.
inline FitControlsVinecop::FitControlsVinecop()
  : FitControlsBicop()
{
  trunc_lvl_ = std::numeric_limits<size_t>::max();
  threshold_ = 0.0;
  tree_criterion_ = "tau";
  select_trunc_lvl_ = false;
  select_threshold_ = false;
  select_families_ = true;
  show_trace_ = false;
  tree_algorithm_ = "mst_prim";
  set_seeds(std::vector<int>());
}

//! @brief Instantiates custom controls for fitting vine copula models.
//! @param family_set The set of copula families to consider (if empty, then
//!     all families are included).
//! @param parametric_method The fit method for parametric families;
//!     possible choices: `"mle"`, `"itau"`.
//! @param nonparametric_method The fit method for the local-likelihood
//!     nonparametric family (TLLs); possible choices: `"constant"`,
//!     `"linear"`, `"quadratic"`.
//! @param nonparametric_mult A factor with which the smoothing parameters
//!     are multiplied.
//! @param nonparametric_grid_size The size of the interpolation grid backing
//!     the local-likelihood nonparametric family (TLLs).
//! @param trunc_lvl Truncation level for truncated vines.
//! @param tree_criterion The criterion for selecting the spanning
//!     tree (`"tau"`, `"rho"`, `"hoeffd"`, `"mcor"`, `"cxi"`, `"joe"`, or
//!     `"custom"`) during the tree-wise structure selection. `"custom"` uses
//!     the callable set through `set_tree_criterion_function()`.
//! @param threshold For thresholded vines (0 = no threshold).
//! @param selection_criterion The selection criterion (`"loglik"`, `"aic"`,
//!     `"bic"`, `"mbic"`, or `"mbicv"`) for the pair copula families.
//! @param weights A vector of weights for the observations.
//! @param psi0 Only for `selection_criterion = "mbic"`, prior probability of
//!     non-independence.
//! @param preselect_families Whether to exclude families before fitting
//!     based on symmetry properties of the data.
//! @param select_trunc_lvl Whether the truncation shall be selected
//!     automatically.
//! @param select_threshold Whether the threshold parameter shall be
//!     selected automatically.
//! @param select_families Whether the families shall be selected
//! automatically, or should the method simply update the parameters for
//! the pair copulas already present in the model.
//! @param show_trace Whether to show a trace of the building progress.
//! @param num_threads Number of concurrent threads to use while fitting
//!     pair copulas within a tree; never uses more than the number
//!     of concurrent threads supported by the implementation. A custom
//!     tree criterion is exempt and always runs on the calling thread.
//! @param tree_algorithm The algorithm for building the spanning
//!     tree (`"mst_prim"`, `"mst_kruskal"`, `"random_weighted"`, or
//!     `"random_unweighted"`) during the tree-wise structure selection.
//!     `"mst_prim"` and `"mst_kruskal"` use Prim's and Kruskal's algorithms
//!     respectively to select the maximum spanning tree, maximizing
//!     the sum of the edge weights (i.e., `tree_criterion`).
//!     `"random_weighted"` and `"random_unweighted"` use Wilson's
//!     algorithm to generate a random spanning tree, either with probability
//!     proportional to the product of the edge weights (weighted) or
//!     uniformly (unweighted).
//! @param allow_rotations Allow rotations for the families when doing
//!     model selection (default: true).
//! @param seeds A vector of random seeds for the random number generator
//!     for parts of the algorithm that are randomized (e.g.,
//!     random tree selection).
inline FitControlsVinecop::FitControlsVinecop(
  std::vector<BicopFamily> family_set,
  std::string parametric_method,
  std::string nonparametric_method,
  double nonparametric_mult,
  size_t nonparametric_grid_size,
  size_t trunc_lvl,
  std::string tree_criterion,
  double threshold,
  std::string selection_criterion,
  const Eigen::VectorXd& weights,
  double psi0,
  bool preselect_families,
  bool select_trunc_lvl,
  bool select_threshold,
  bool select_families,
  bool show_trace,
  size_t num_threads,
  std::string tree_algorithm,
  bool allow_rotations,
  std::vector<int> seeds)
  : FitControlsBicop(std::move(family_set),
                     std::move(parametric_method),
                     std::move(nonparametric_method),
                     nonparametric_mult,
                     nonparametric_grid_size,
                     std::move(selection_criterion),
                     weights,
                     psi0,
                     preselect_families,
                     allow_rotations,
                     num_threads)
{
  set_trunc_lvl(trunc_lvl);
  set_tree_criterion(std::move(tree_criterion));
  set_threshold(threshold);
  set_select_trunc_lvl(select_trunc_lvl);
  set_select_threshold(select_threshold);
  set_select_families(select_families);
  set_show_trace(show_trace);
  set_tree_algorithm(std::move(tree_algorithm));
  set_seeds(std::move(seeds));
}

//! @brief Instantiates custom controls for fitting vine copula models.
//!
//! @param controls See `FitControlsBicop()`.
//! @param trunc_lvl Truncation level for truncated vines.
//! @param tree_criterion The criterion for selecting the spanning
//!     tree (`"tau"`, `"rho"`, `"hoeffd"`, `"mcor"`, `"cxi"`, `"joe"`, or
//!     `"custom"`) during the tree-wise structure selection. `"custom"` uses
//!     the callable set through `set_tree_criterion_function()`.
//! @param threshold For thresholded vines (`0` = no threshold).
//! @param show_trace Whether to show a trace of the building progress.
//! @param select_trunc_lvl Whether the truncation shall be selected
//!     automatically.
//! @param select_threshold Whether the threshold parameter shall be
//!     selected automatically.
//! @param select_families Whether the families shall be selected
//! automatically, or should the method simply update the parameters for
//! the pair copulas already present in the model.
//! @param tree_algorithm The algorithm for building the spanning
//!     tree (`"mst_prim"`, `"mst_kruskal"`, `"random_weighted"`, or
//!     `"random_unweighted"`) during the tree-wise structure selection.
//!     `"mst_prim"` and `"mst_kruskal"` use Prim's and Kruskal's algorithms
//!     respectively to select the maximum spanning tree, maximizing
//!     the sum of the edge weights (i.e., `tree_criterion`).
//!     `"random_weighted"` and `"random_unweighted"` use Wilson's
//!     algorithm to generate a random spanning tree, either with probability
//!     proportional to the product of the edge weights (weighted) or
//!     uniformly (unweighted).
//! @param seeds A vector of random seeds for the random number generator
//!     for parts of the algorithm that are randomized (e.g.,
//!     random tree selection).
inline FitControlsVinecop::FitControlsVinecop(const FitControlsBicop& controls,
                                              size_t trunc_lvl,
                                              std::string tree_criterion,
                                              double threshold,
                                              bool select_trunc_lvl,
                                              bool select_threshold,
                                              bool select_families,
                                              bool show_trace,
                                              std::string tree_algorithm,
                                              std::vector<int> seeds)
  : FitControlsBicop(controls)
{
  set_trunc_lvl(trunc_lvl);
  set_tree_criterion(std::move(tree_criterion));
  set_threshold(threshold);
  set_select_trunc_lvl(select_trunc_lvl);
  set_select_threshold(select_threshold);
  set_select_families(select_families);
  set_show_trace(show_trace);
  set_tree_algorithm(std::move(tree_algorithm));
  set_seeds(std::move(seeds));
}

//! @brief Instantiates the controls from a configuration object.
//! @param config The configuration object.
inline FitControlsVinecop::FitControlsVinecop(const FitControlsConfig& config)
  : FitControlsBicop(config)
{
  if (config.trunc_lvl.has_value()) {
    set_trunc_lvl(*config.trunc_lvl);
  }
  if (config.tree_criterion.has_value()) {
    set_tree_criterion(*config.tree_criterion);
  }
  if (config.tree_criterion_function.has_value()) {
    set_tree_criterion_function(*config.tree_criterion_function);
  }
  if (config.threshold.has_value()) {
    set_threshold(*config.threshold);
  }
  if (config.select_trunc_lvl.has_value()) {
    set_select_trunc_lvl(*config.select_trunc_lvl);
  }
  if (config.select_threshold.has_value()) {
    set_select_threshold(*config.select_threshold);
  }
  if (config.select_families.has_value()) {
    set_select_families(*config.select_families);
  }
  if (config.show_trace.has_value()) {
    set_show_trace(*config.show_trace);
  }
  if (config.tree_algorithm.has_value()) {
    set_tree_algorithm(*config.tree_algorithm);
  }
  if (config.seeds.has_value()) {
    set_seeds(*config.seeds);
  }
  if (config.conditioning_set.has_value()) {
    set_conditioning_set(*config.conditioning_set);
  }
}

//! @name Sanity checks
//! @{
inline void
FitControlsVinecop::check_tree_criterion(std::string tree_criterion)
{
  if (!tools_stl::is_member(
        std::move(tree_criterion),
        { "tau", "rho", "joe", "hoeffd", "mcor", "cxi", "custom" })) {
    throw std::runtime_error("tree_criterion must be one of "
                             "'tau', 'rho', 'hoeffd', 'mcor', 'cxi', 'joe', "
                             "or 'custom'");
  }
}

inline void
FitControlsVinecop::check_threshold(double threshold)
{
  if (threshold < 0 || threshold > 1) {
    throw std::runtime_error("threshold should be in [0,1]");
  }
}

inline void
FitControlsVinecop::check_conditioning_set(
  const std::vector<size_t>& conditioning_set)
{
  // Dimension-free checks only (the vine dimension d is not known here; the
  // upper bounds max(set) <= d and |set| <= d - 1 are checked in
  // Vinecop::select()).
  for (auto v : conditioning_set) {
    if (v < 1) {
      throw std::runtime_error(
        "conditioning_set entries must be >= 1 (1-based variable indices)");
    }
  }
  auto sorted = conditioning_set;
  std::sort(sorted.begin(), sorted.end());
  if (std::adjacent_find(sorted.begin(), sorted.end()) != sorted.end()) {
    throw std::runtime_error("conditioning_set must not contain duplicates");
  }
}
//! @}

//! @name Getters and setters.
//! @{

//! @brief Gets the truncation level.
inline size_t
FitControlsVinecop::get_trunc_lvl() const
{
  return trunc_lvl_;
}

//! @brief Sets the truncation level.
inline void
FitControlsVinecop::set_trunc_lvl(size_t trunc_lvl)
{
  trunc_lvl_ = trunc_lvl;
}

//! @brief Gets whether to select the truncation level automatically.
inline bool
FitControlsVinecop::get_select_trunc_lvl() const
{
  return select_trunc_lvl_;
}

//! @brief Sets whether to select the truncation level automatically.
inline void
FitControlsVinecop::set_select_trunc_lvl(bool select_trunc_lvl)
{
  select_trunc_lvl_ = select_trunc_lvl;
}

//! @brief Gets whether to select the families automatically.
inline bool
FitControlsVinecop::get_select_families() const
{
  return select_families_;
}

//! @brief Sets whether to select the families automatically.
inline void
FitControlsVinecop::set_select_families(bool select_families)
{
  select_families_ = select_families;
}

//! @brief Gets the criterion for tree selection.
inline std::string
FitControlsVinecop::get_tree_criterion() const
{
  return tree_criterion_;
}

//! @brief Sets the criterion for tree selection.
inline void
FitControlsVinecop::set_tree_criterion(std::string tree_criterion)
{
  check_tree_criterion(tree_criterion);
  tree_criterion_ = std::move(tree_criterion);
}

//! @brief Gets the custom criterion function for tree selection.
inline TreeCriterionFunction
FitControlsVinecop::get_tree_criterion_function() const
{
  return tree_criterion_function_;
}

//! @brief Sets the custom criterion function for tree selection.
inline void
FitControlsVinecop::set_tree_criterion_function(
  TreeCriterionFunction tree_criterion_function)
{
  tree_criterion_function_ = std::move(tree_criterion_function);
}

//! @brief Gets the threshold parameter.
inline double
FitControlsVinecop::get_threshold() const
{
  return threshold_;
}

//! @brief Sets the threshold parameter.
inline void
FitControlsVinecop::set_threshold(double threshold)
{
  check_threshold(threshold);
  threshold_ = threshold;
}

//! @brief Gets whether to show a trace is during fitting.
inline bool
FitControlsVinecop::get_show_trace() const
{
  return show_trace_;
}

//! @brief Gets whether to show a trace is during fitting.
inline void
FitControlsVinecop::set_show_trace(bool show_trace)
{
  show_trace_ = show_trace;
}

//! @brief Gets whether to select the threshold automatically.
inline bool
FitControlsVinecop::get_select_threshold() const
{
  return select_threshold_;
}

//! @brief Gets the maximum spanning tree algorithm.
inline std::string
FitControlsVinecop::get_tree_algorithm() const
{
  return tree_algorithm_;
}

//! @brief Gets the random seeds for the random number generator.
inline std::vector<int>
FitControlsVinecop::get_seeds() const
{
  return seeds_;
}

//! @brief Gets the conditioning set for conditioning-aware selection.
inline std::vector<size_t>
FitControlsVinecop::get_conditioning_set() const
{
  return conditioning_set_;
}

//! @brief Sets the conditioning set for conditioning-aware selection.
inline void
FitControlsVinecop::set_conditioning_set(std::vector<size_t> conditioning_set)
{
  check_conditioning_set(conditioning_set);
  conditioning_set_ = std::move(conditioning_set);
}

//! @brief Gets the random number generator.
inline boost::random::mt19937
FitControlsVinecop::get_rng() const
{
  return rng_;
}

//! @brief Sets whether to select the threshold automatically.
inline void
FitControlsVinecop::set_select_threshold(bool select_threshold)
{
  select_threshold_ = select_threshold;
}

inline bool
FitControlsVinecop::needs_sparse_select() const
{
  return (select_trunc_lvl_ | select_threshold_);
}

//! @brief Gets the fit controls for bivariate fitting.
inline FitControlsBicop
FitControlsVinecop::get_fit_controls_bicop() const
{
  FitControlsBicop controls_bicop(get_family_set(),
                                  get_parametric_method(),
                                  get_nonparametric_method(),
                                  get_nonparametric_mult(),
                                  get_nonparametric_grid_size(),
                                  get_selection_criterion(),
                                  get_weights(),
                                  get_psi0(),
                                  get_preselect_families());
  return controls_bicop;
}

//! @brief Sets the fit controls for bivariate fitting.
inline void
FitControlsVinecop::set_fit_controls_bicop(const FitControlsBicop& controls)
{
  set_family_set(controls.get_family_set());
  set_parametric_method(controls.get_parametric_method());
  set_nonparametric_method(controls.get_nonparametric_method());
  set_nonparametric_mult(controls.get_nonparametric_mult());
  set_nonparametric_grid_size(controls.get_nonparametric_grid_size());
  set_selection_criterion(get_selection_criterion());
  set_preselect_families(controls.get_preselect_families());
}

//! @brief Sets the maximum spanning tree algorithm.
inline void
FitControlsVinecop::set_tree_algorithm(std::string tree_algorithm)
{
  if (!tools_stl::is_member(tree_algorithm,
                            { "mst_prim",
                              "mst_kruskal",
                              "random_weighted",
                              "random_unweighted" })) {
    throw std::runtime_error(
      "tree_algorithm must be one of 'mst_prim', 'mst_kruskal', "
      "'random_weighted', or 'random_unweighted'");
  }
  tree_algorithm_ = std::move(tree_algorithm);
}

//! @brief Sets the random seeds for the random number generator.
inline void
FitControlsVinecop::set_seeds(std::vector<int> seeds)
{
  if (seeds.empty()) {
    // no seeds provided, seed randomly
    std::random_device rd{};
    seeds = std::vector<int>(20);
    std::generate(
      seeds.begin(), seeds.end(), [&]() { return static_cast<int>(rd()); });
  }
  seeds_ = seeds;
  boost::random::seed_seq seq(seeds.begin(), seeds.end());
  rng_.seed(seq);
}

//! @}

//! @brief Summarizes the controls into a string (can be used for printing).
inline std::string
FitControlsVinecop::str() const
{
  std::stringstream controls_str;

  controls_str << str_internal(false);
  controls_str << "Truncation level: "
               << (get_trunc_lvl() == std::numeric_limits<size_t>::max()
                     ? "none (default)"
                     : std::to_string(get_trunc_lvl()))
               << std::endl;
  controls_str << "Tree criterion: " << get_tree_criterion() << std::endl;
  controls_str << "Threshold: " << get_threshold() << std::endl;
  controls_str << "Select truncation level: "
               << static_cast<std::string>(get_select_trunc_lvl() ? "yes"
                                                                  : "no")
               << std::endl;
  controls_str << "Select threshold: "
               << static_cast<std::string>(get_select_threshold() ? "yes"
                                                                  : "no")
               << std::endl;
  controls_str << "Select families: "
               << static_cast<std::string>(get_select_families() ? "yes" : "no")
               << std::endl;
  controls_str << "Show trace: "
               << static_cast<std::string>(get_show_trace() ? "yes" : "no")
               << std::endl;
  controls_str << "Number of threads: "
               << (get_num_threads() == 0 ? 1 : get_num_threads()) << std::endl;
  controls_str << "MST algorithm: " << get_tree_algorithm() << std::endl;
  controls_str << "Conditioning set: ";
  if (get_conditioning_set().empty()) {
    controls_str << "none (default)";
  } else {
    for (auto v : get_conditioning_set())
      controls_str << v << " ";
  }
  controls_str << std::endl;
  return controls_str.str();
}

}
