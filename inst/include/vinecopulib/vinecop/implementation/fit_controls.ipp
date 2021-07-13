// Copyright © 2016-2021 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <stdexcept>
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
  show_trace_ = false;
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
//! @param trunc_lvl Truncation level for truncated vines.
//! @param tree_criterion The criterion for selecting the maximum spanning
//!     tree ("tau", "hoeffd", "rho", and "mcor" implemented so far).
//! @param threshold For thresholded vines (0 = no threshold).
//! @param selection_criterion The selection criterion (`"loglik"`, `"aic"`
//!     or `"bic"`).
//! @param weights A vector of weights for the observations.
//! @param psi0 Only for `selection_criterion = "mbic", prior probability of
//!     non-independence.
//! @param preselect_families Whether to exclude families before fitting
//!     based on symmetry properties of the data.
//! @param select_trunc_lvl Whether the truncation shall be selected
//!     automatically.
//! @param select_threshold Whether the threshold parameter shall be
//!     selected automatically.
//! @param show_trace Whether to show a trace of the building progress.
//! @param num_threads Number of concurrent threads to use while fitting
//!     pair copulas within a tree; never uses more than the number
//!     of concurrent threads supported by the implementation.
inline FitControlsVinecop::FitControlsVinecop(
  std::vector<BicopFamily> family_set,
  std::string parametric_method,
  std::string nonparametric_method,
  double nonparametric_mult,
  size_t trunc_lvl,
  std::string tree_criterion,
  double threshold,
  std::string selection_criterion,
  const Eigen::VectorXd& weights,
  double psi0,
  bool preselect_families,
  bool select_trunc_lvl,
  bool select_threshold,
  bool show_trace,
  size_t num_threads)
  : FitControlsBicop(family_set,
                     parametric_method,
                     nonparametric_method,
                     nonparametric_mult,
                     selection_criterion,
                     weights,
                     psi0,
                     preselect_families)
{
  set_trunc_lvl(trunc_lvl);
  set_tree_criterion(tree_criterion);
  set_threshold(threshold);
  set_select_trunc_lvl(select_trunc_lvl);
  set_select_threshold(select_threshold);
  set_show_trace(show_trace);
  set_num_threads(num_threads);
}

//! @brief Instantiates custom controls for fitting vine copula models.
//! @param trunc_lvl Truncation level for truncated vines.
//! @param tree_criterion The criterion for selecting the maximum spanning
//!     tree ("tau", "hoeffd" and "rho" implemented so far).
//! @param threshold For thresholded vines (0 = no threshold).
//! @param show_trace Whether to show a trace of the building progress.
//! @param select_trunc_lvl Whether the truncation shall be selected
//!     automatically.
//! @param select_threshold Whether the threshold parameter shall be
//!     selected automatically.
//! @param controls See FitControlsBicop.
//! @param num_threads Number of concurrent threads to use while fitting
//!     pair copulas within a tree; never uses more than the number returned
//!     by `std::thread::hardware_concurrency()``.
inline FitControlsVinecop::FitControlsVinecop(const FitControlsBicop& controls,
                                              size_t trunc_lvl,
                                              std::string tree_criterion,
                                              double threshold,
                                              bool select_trunc_lvl,
                                              bool select_threshold,
                                              bool show_trace,
                                              size_t num_threads)
  : FitControlsBicop(controls)
{
  set_trunc_lvl(trunc_lvl);
  set_tree_criterion(tree_criterion);
  set_threshold(threshold);
  set_select_trunc_lvl(select_trunc_lvl);
  set_select_threshold(select_threshold);
  set_show_trace(show_trace);
  set_num_threads(num_threads);
}

//! @name Sanity checks
//! @{
inline void
FitControlsVinecop::check_tree_criterion(std::string tree_criterion)
{
  if (!tools_stl::is_member(tree_criterion,
                            { "tau", "rho", "joe", "hoeffd", "mcor" })) {
    throw std::runtime_error("tree_criterion must be one of "
                             "'tau', 'rho', 'hoeffd', 'mcor', or 'joe'");
  }
}

inline void
FitControlsVinecop::check_threshold(double threshold)
{
  if (threshold < 0 || threshold > 1) {
    throw std::runtime_error("threshold should be in [0,1]");
  }
}
//! @}

//! @name Getters and setters.
//! @{

//! returns the truncation level.
inline size_t
FitControlsVinecop::get_trunc_lvl() const
{
  return trunc_lvl_;
}

//! Sets the truncation level.
inline void
FitControlsVinecop::set_trunc_lvl(size_t trunc_lvl)
{
  trunc_lvl_ = trunc_lvl;
}

//! returns whether to select the truncation level automatically.
inline bool
FitControlsVinecop::get_select_trunc_lvl() const
{
  return select_trunc_lvl_;
}

//! Sets whether to select the truncation level automatically.
inline void
FitControlsVinecop::set_select_trunc_lvl(bool select_trunc_lvl)
{
  select_trunc_lvl_ = select_trunc_lvl;
}

//! returns the criterion for tree selection.
inline std::string
FitControlsVinecop::get_tree_criterion() const
{
  return tree_criterion_;
}

//! Sets the criterion for tree selection.
inline void
FitControlsVinecop::set_tree_criterion(std::string tree_criterion)
{
  check_tree_criterion(tree_criterion);
  tree_criterion_ = tree_criterion;
}

//! returns the threshold parameter.
inline double
FitControlsVinecop::get_threshold() const
{
  return threshold_;
}

//! Sets the threshold parameter.
inline void
FitControlsVinecop::set_threshold(double threshold)
{
  check_threshold(threshold);
  threshold_ = threshold;
}

//! returns whether to show a trace is during fitting.
inline bool
FitControlsVinecop::get_show_trace() const
{
  return show_trace_;
}

//! returns whether to show a trace is during fitting.
inline void
FitControlsVinecop::set_show_trace(bool show_trace)
{
  show_trace_ = show_trace;
}

//! returns whether to select the threshold automatically.
inline bool
FitControlsVinecop::get_select_threshold() const
{
  return select_threshold_;
}

//! Sets whether to select the threshold automatically.
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

//! Returns the fit controls for bivariate fitting.
inline FitControlsBicop
FitControlsVinecop::get_fit_controls_bicop() const
{
  FitControlsBicop controls_bicop(get_family_set(),
                                  get_parametric_method(),
                                  get_nonparametric_method(),
                                  get_nonparametric_mult(),
                                  get_selection_criterion(),
                                  get_weights(),
                                  get_psi0(),
                                  get_preselect_families());
  return controls_bicop;
}

//! Sets the fit controls for bivariate fitting.
inline void
FitControlsVinecop::set_fit_controls_bicop(FitControlsBicop controls)
{
  set_family_set(controls.get_family_set());
  set_parametric_method(controls.get_parametric_method());
  set_selection_criterion(get_selection_criterion());
  set_preselect_families(controls.get_preselect_families());
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
               << static_cast<std::string>(get_select_trunc_lvl() ? "yes"
                                                                  : "no")
               << std::endl;
  controls_str << "Show trace: "
               << static_cast<std::string>(get_show_trace() ? "yes" : "no")
               << std::endl;
  controls_str << "Number of threads: "
               << (get_num_threads() == 0 ? 1 : get_num_threads()) << std::endl;
  return controls_str.str().c_str();
}

}
