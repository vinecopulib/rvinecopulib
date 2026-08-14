// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <utility>
#include <vinecopulib/bicop/class.hpp>
#include <vinecopulib/misc/tools_interface.hpp>
#include <vinecopulib/misc/tools_serialization.hpp>
#include <vinecopulib/misc/tools_stats.hpp>
#include <vinecopulib/misc/tools_stl.hpp>
#include <vinecopulib/vinecop/tools_select.hpp>

#include <stdexcept>

namespace vinecopulib {

//! @brief Instantiates a D-vine with all pair-copulas set to independence.
//! @param d The dimension (= number of variables) of the model.
inline Vinecop::Vinecop(const size_t d)
  : Vinecop(RVineStructure(d, static_cast<size_t>(0)))
{
}

//! @brief Instantiates an arbitrary vine copula model.
//! @param structure An `RVineStructure` object specifying the vine structure.
//! @param pair_copulas `Bicop` objects specifying the pair-copulas, namely
//!     a nested list such that `pc_store[t][e]` contains a `Bicop`
//!     object for the pair copula corresponding to tree `t` and edge `e`.
//! @param var_types Strings specifying the types of the variables,
//!   e.g., `("c", "d")` means first variable continuous, second discrete.
//!   If empty, then all variables are set as continuous.
inline Vinecop::Vinecop(const RVineStructure& structure,
                        const std::vector<std::vector<Bicop>>& pair_copulas,
                        const std::vector<std::string>& var_types)
  : d_(structure.get_dim())
  , rvine_structure_(structure)
{

  if (!pair_copulas.empty()) {
    set_all_pair_copulas(pair_copulas);
  }

  if (!var_types.empty()) {
    set_var_types(var_types);
  } else {
    set_continuous_var_types();
  }
}

//! @brief Instantiates an arbitrary vine copula model.
//! @param matrix An R-vine matrix specifying the vine structure.
//! @param pair_copulas `Bicop` objects specifying the pair-copulas, namely
//!     a nested list such that `pc_store[t][e]` contains a `Bicop`
//!     object for the pair copula corresponding to tree `t` and edge `e`.
//! @param var_types Strings specifying the types of the variables,
//!   e.g., `("c", "d")` means first variable continuous, second discrete.
//!   If empty, then all variables are set as continuous.
inline Vinecop::Vinecop(
  const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& matrix,
  const std::vector<std::vector<Bicop>>& pair_copulas,
  const std::vector<std::string>& var_types)
  : Vinecop(RVineStructure(matrix), pair_copulas, var_types)
{
}

//! @brief Instantiates from data.
//!
//! @details Equivalent to creating a default `Vinecop()` and then
//! selecting the model using `select()`.
//!
//! @param data An \f$ n \times d \f$ matrix of observations.
//! @param structure An RVineStructure object specifying the vine structure.
//!    If empty, then it is selected as part of the fit.
//! @param var_types Strings specifying the types of the variables,
//!   e.g., `("c", "d")` means first variable continuous, second discrete.
//!   If empty, then all variables are set as continuous.
//! @param controls See `FitControlsVinecop()`.
inline Vinecop::Vinecop(const Eigen::MatrixXd& data,
                        const RVineStructure& structure,
                        const std::vector<std::string>& var_types,
                        const FitControlsVinecop& controls)
{
  check_enough_data(data);
  if (structure.get_dim() > 1) {
    d_ = structure.get_dim();
    rvine_structure_ = structure;
  } else {
    if (!var_types.empty()) {
      d_ = var_types.size();
    } else {
      d_ = data.cols();
    }
    rvine_structure_ = RVineStructure(d_, static_cast<size_t>(0));
  }
  if (var_types.empty()) {
    set_continuous_var_types();
  } else {
    set_var_types(var_types);
  }
  check_weights_size(controls.get_weights(), data);
  select(data, controls);
}

//! @brief Instantiates from data.
//!
//! @details Equivalent to creating a default `Vinecop()` and
//! then selecting the model using `select()`.
//!
//! @param data An \f$ n \times d \f$ matrix of observations.
//! @param matrix Either an R-vine structure matrix, see `select()`, or an
//!     empty matrix, in which case the structure is selected as part of the
//!     fit. To select the structure, prefer `Vinecop(data)`.
//! @param var_types Strings specifying the types of the variables,
//!   e.g., `("c", "d")` means first variable continuous, second discrete.
//!   If empty, then all variables are set as continuous.
//! @param controls See `FitControlsVinecop()`.
inline Vinecop::Vinecop(
  const Eigen::MatrixXd& data,
  const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& matrix,
  const std::vector<std::string>& var_types,
  const FitControlsVinecop& controls)
  : Vinecop(data, RVineStructure(matrix), var_types, controls)
{
}

//! @brief Instantiates from a nlohmann::json object.
//! @param input The nlohmann::json object to convert from
//! (see `to_json()` for the structure of the input).
//! @param check Whether to check if the `"structure"` node represents
//!      a valid R-vine structure.
inline Vinecop::Vinecop(const nlohmann::json& input, const bool check)
{
  rvine_structure_ = RVineStructure(input["structure"], check);
  d_ = static_cast<size_t>(rvine_structure_.get_dim());
  size_t trunc_lvl = rvine_structure_.get_trunc_lvl();

  nlohmann::json pcs_json = input["pair copulas"];
  for (size_t tree = 0; tree < std::min(d_, trunc_lvl); ++tree) {
    nlohmann::json tree_json;
    try {
      tree_json = pcs_json["tree" + std::to_string(tree)];
    } catch (...) {
      break; // vine was truncated, no more trees to parse
    }
    // reserve space for pair copulas of this tree
    pair_copulas_.resize(tree + 1);
    pair_copulas_[tree].resize(d_ - tree - 1);

    for (size_t edge = 0; edge < d_ - tree - 1; ++edge) {
      nlohmann::json pc_json = tree_json["pc" + std::to_string(edge)];
      pair_copulas_[tree][edge] = Bicop(pc_json);
    }
  }

  // try block for backwards compatibility
  try {
    var_types_ =
      tools_serialization::json_to_vector<std::string>(input["var_types"]);
    n_discrete_ = 0;
    for (const auto& t : var_types_) {
      n_discrete_ += (t == "d");
    }
    nobs_ = static_cast<size_t>(input["nobs_"]);
    threshold_ = static_cast<double>(input["threshold"]);
    loglik_ = static_cast<double>(input["loglik"]);
  } catch (...) {
  }
}

//! @brief Instantiates from a JSON or CBOR file.
//!
//! @details Files ending in `.cbor` are read as CBOR. All other filenames are
//! read as JSON for backwards compatibility. The input contains 2 attributes:
//! `"structure"` for the vine structure, which itself contains attributes
//! `"array"` for the structure triangular array and `"order"` for the order
//! vector, and `"pair copulas"`.
//! `"pair copulas"` contains a list of attributes for the trees
//! (`"tree1"`, `"tree2"`, etc), each containing
//! a list of attributes for the edges (`"pc1"`, `"pc2"`, etc).
//! See the corresponding method of `Bicop` objects for the encoding of
//! pair-copulas.
//!
//! @param filename The name of the file to read.
//! @param check Whether to check if the `"structure"` node of the input
//! represents a valid R-vine structure.
inline Vinecop::Vinecop(const std::string& filename, const bool check)
  : Vinecop(tools_serialization::file_to_json(filename), check)
{
}

//! @brief Converts the copula into a nlohmann::json object.
//!
//! @details The `nlohmann::json` object contains two nodes : `"structure"`
//! for the vine
//! structure, which itself contains nodes `"array"` for the structure
//! triangular array and `"order"` for the order vector, and `"pair copulas"`.
//! The former two encode the R-Vine structure and the latter is a list of
//! child nodes for the trees (`"tree1"`, `"tree2"`, etc), each containing
//! a list of child nodes for the edges (`"pc1"`, `"pc2"`, etc).
//! See Bicop::to_json() for the encoding of pair-copulas.
//!
//! @return the nlohmann::json object containing the copula.
inline nlohmann::json
Vinecop::to_json() const
{
  nlohmann::json pair_copulas;
  for (size_t tree = 0; tree < pair_copulas_.size(); ++tree) {
    nlohmann::json tree_json;
    for (size_t edge = 0; edge < d_ - tree - 1; ++edge) {
      tree_json["pc" + std::to_string(edge)] =
        pair_copulas_[tree][edge].to_json();
    }
    pair_copulas["tree" + std::to_string(tree)] = tree_json;
  }

  nlohmann::json output;
  output["pair copulas"] = pair_copulas;
  auto structure_json = rvine_structure_.to_json();
  output["structure"] = structure_json;
  output["var_types"] = tools_serialization::vector_to_json(var_types_);
  output["nobs_"] = nobs_;
  output["threshold"] = threshold_;
  output["loglik"] = loglik_;

  return output;
}

//! @brief Writes the copula object into a JSON or CBOR file.
//!
//! @details Filenames ending in `.cbor` are written as CBOR. All other
//! filenames are written as JSON for backwards compatibility. The output
//! contains 2 attributes: `"structure"` for the vine structure, which itself
//! contains attributes `"array"` for the structure triangular array and
//! `"order"` for the order vector, and `"pair copulas"`.
//! `"pair copulas"` contains a list of attributes for the trees
//! (`"tree1"`, `"tree2"`, etc), each containing
//! a list of attributes for the edges (`"pc1"`, `"pc2"`, etc).
//! See `Bicop::to_file()` objects for the encoding of
//! pair-copulas.
//!
//! @param filename The name of the JSON file to write.
inline void
Vinecop::to_file(const std::string& filename) const
{
  tools_serialization::json_to_file(filename, this->to_json());
}

//! @brief Initializes object for storing pair copulas.
//!
//! @param d Dimension of the vine copula.
//! @param trunc_lvl A truncation level (optional).
//! @return A nested list such that `pc_store[t][e]` contains a Bicop.
//!     object for the pair copula corresponding to tree `t` and edge `e`.
inline std::vector<std::vector<Bicop>>
Vinecop::make_pair_copula_store(const size_t d, const size_t trunc_lvl)
{
  return tools_select::VinecopSelector::make_pair_copula_store(d, trunc_lvl);
}

//! @brief Automatically fits and selects a vine copula model.
//!
//! @details This method can be used to select either the pair-copulas only,
//! or the pair-copulas and the structure. The latter is done by specifying
//! a truncation level in the controls. The method then selects the structure
//! and fits the pair-copulas for all trees up to the truncation level.

//! In other words, `select()` behaves differently depending on its current
//! truncation level and the truncation level specified in the controls,
//! respectively called `trunc_lvl` and `controls.trunc_lvl` in what follows.
//! Essentially, `controls.trunc_lvl` defines the object's truncation level
//! after calling `select()`:
//!
//!   - If `controls.trunc_lvl <= trunc_lvl`, the families and parameters for
//!     all pairs in trees smaller or equal to `controls.trunc_lvl`
//!     are selected, using the current structure.
//!   - If `controls.trunc_lvl > trunc_lvl`, `select()` behaves as above for
//!     all trees that are smaller or equal to `trunc_lvl`, and then it selects
//!     the structure for higher trees along with the families and parameters.
//!     This includes the case where `trunc_lvl = 0`, namely where the
//!     structure is fully unspecified.
//!
//! Selection of the structure is performed using the algorithm of
//! Dissmann, J. F., E. C. Brechmann, C. Czado, and D. Kurowicka (2013).
//! *Selecting and estimating regular vine copulae and application to
//! financial returns.* Computational Statistics & Data Analysis, 59 (1),
//! 52-69.
//! The dependence measure used to select trees (default: Kendall's tau) is
//! corrected for ties (see the [wdm](https://github.com/tnagler/wdm) library).
//! The dependence measure can be changed using `controls.tree_criterion`,
//! which can be set to `"tau"`, `"rho"`, `"hoeffd"`, `"mcor"`, `"joe"`, or
//! `"custom"`. The last one uses the callable supplied through
//! `controls.tree_criterion_function`, which is always called on the thread
//! that starts the fit, so it need not be thread safe; the pair-copula fits
//! still use `controls.num_threads` threads.
//! Both Prim's (default: `"mst_prim"`) and Kruskal's (`"mst_kruskal"`)
//! algorithms are available through `controls.tree_algorithm` for the
//! maximum spanning tree selection.
//! An alternative to the maximum spanning tree selection is to use random
//! spanning trees, which can be selected using `controls.tree_algorithm` and
//! come in two flavors, both using Wilson's algorithm loop erased random walks:
//!
//!   - `"random_weighted"` generates a random spanning tree with probability
//!     proportional to the product of the weights (i.e., the dependence) of
//!     the edges in the tree.
//!   - `"random_unweighted"` generates a random spanning tree uniformly over
//!     all spanning trees satisfying the proximity condition.
//!
//! If the `controls` object has been instantiated with
//! `select_families = false`, then the method simply updates the parameters of
//! the pair-copulas without selecting the families or the structure.
//! In this case, this is equivalent to calling `fit()` for each pair-copula,
//! albeit potentially in parallel if `num_threads > 1`.
//!
//! When at least one variable is discrete, two types of
//! "observations" are required: the first \f$ n \times d \f$ block contains
//! realizations of \f$ F_Y(Y), F_X(X) \f$; the second \f$ n \times d \f$ block
//! contains realizations of \f$ F_Y(Y^-), F_X(X^-), ... \f$. The minus
//! indicates a left-sided limit of the cdf. For continuous variables the left
//! limit and the cdf itself coincide. For, e.g., an integer-valued variable, it
//! holds \f$ F_Y(Y^-) = F_Y(Y - 1) \f$. Continuous variables in the second
//! block can be omitted.
//!
//! If there are missing data (i.e., NaN entries), incomplete observations are
//! discarded before fitting a pair-copula. This is done on a pair-by-pair basis
//! so that the maximal available information is used.
//!
//!
//! @param data \f$ n \times 2d \f$ or \f$ n \times (d + k) \f$ matrix of
//!   observations, where \f$ k \f$ is the number of discrete variables.
//! @param controls The controls to the algorithm (see `FitControlsVinecop()`).
inline void
Vinecop::select(const Eigen::MatrixXd& data, const FitControlsVinecop& controls)
{
  if (controls.get_select_families()) {
    check_tree_criterion_function(controls);
    check_data(data);
    if (d_ == 1) {
      loglik_ = 0;
      nobs_ = data.rows();
      return;
    }

    auto conditioning_set = controls.get_conditioning_set();
    if (!conditioning_set.empty()) {
      check_conditioning_set(conditioning_set, controls);
    }

    Eigen::MatrixXd u = collapse_data(data);

    tools_select::VinecopSelector selector(
      u, rvine_structure_, controls, var_types_);
    if (controls.needs_sparse_select()) {
      selector.sparse_select_all_trees(u);
    } else {
      selector.select_all_trees(u);
    }
    finalize_fit(selector);

    // Relabel the fitted vine (a value-preserving re-orientation, no re-fit) so
    // the conditioning set becomes the tail of the order and can be conditioned
    // on with simulate_conditional(). The constrained first tree makes this
    // achievable.
    if (!conditioning_set.empty()) {
      reorient(conditioning_set);
    }
  } else {
    fit(data, controls.get_fit_controls_bicop(), controls.get_num_threads());
  }
}

//! @brief Validates the conditioning set against the vine dimension and
//! selection controls (called from `select()` where `d_` is known).
inline void
Vinecop::check_conditioning_set(const std::vector<size_t>& conditioning_set,
                                const FitControlsVinecop& controls) const
{
  if (conditioning_set.size() >= d_) {
    throw std::runtime_error(
      "conditioning_set must contain at most d - 1 variables.");
  }
  for (auto v : conditioning_set) {
    if ((v < 1) || (v > d_)) {
      throw std::runtime_error(
        "conditioning_set entries must be in 1, ..., d.");
    }
  }
  auto algo = controls.get_tree_algorithm();
  if ((algo != "mst_prim") && (algo != "mst_kruskal")) {
    throw std::runtime_error(
      "conditioning-aware selection requires an MST tree_algorithm "
      "('mst_prim' or 'mst_kruskal').");
  }
  // v1: re-orientation operates on the full R-vine matrix, so structural
  // truncation is not supported (thresholding, which keeps the structure full,
  // is fine).
  if (controls.get_select_trunc_lvl() || (controls.get_trunc_lvl() < d_ - 1)) {
    throw std::runtime_error(
      "conditioning-aware selection does not support truncation "
      "(trunc_lvl / select_trunc_lvl) in this version.");
  }
}

//! @brief Validates the custom edge-weight function against the criterion it
//! belongs to (called from `select()`, the single entry point to structure
//! selection; the two fields can be set in either order, so the pairing can
//! only be checked once the fit starts).
inline void
Vinecop::check_tree_criterion_function(const FitControlsVinecop& controls)
{
  bool is_custom = (controls.get_tree_criterion() == "custom");
  bool has_function = static_cast<bool>(controls.get_tree_criterion_function());
  if (is_custom && !has_function) {
    throw std::runtime_error("tree_criterion = \"custom\" requires a "
                             "tree_criterion_function callable");
  }
  if (has_function && !is_custom) {
    throw std::runtime_error(
      "a tree_criterion_function was provided, but tree_criterion is \"" +
      controls.get_tree_criterion() +
      "\"; set tree_criterion = \"custom\" to use it");
  }
}

inline Vinecop::VinecopView::VinecopView(const Vinecop& vinecop)
  : vinecop_(&vinecop)
  , reorientation_(nullptr)
{
}

inline Vinecop::VinecopView::VinecopView(const Vinecop& vinecop,
                                         const ReorientationMap& reorientation)
  : vinecop_(&vinecop)
  , reorientation_(reorientation.identity ? nullptr : &reorientation)
{
}

inline const RVineStructure&
Vinecop::VinecopView::get_structure() const
{
  return reorientation_ ? reorientation_->structure
                        : vinecop_->rvine_structure_;
}

inline BicopView
Vinecop::VinecopView::get_pair_copula(size_t tree, size_t edge) const
{
  if (!reorientation_)
    return BicopView(vinecop_->pair_copulas_[tree][edge]);

  auto location = reorientation_->pair_copulas(tree, edge);
  return BicopView(vinecop_->pair_copulas_[location.tree][location.edge],
                   location.flipped);
}

inline Vinecop::ReorientationMap
Vinecop::make_reorientation_map(
  const std::vector<size_t>& conditioning_set) const
{
  size_t k = conditioning_set.size();

  // ---- validation ----
  if ((k == 0) || (k >= d_)) {
    throw std::runtime_error(
      "conditioning_set must contain between 1 and d - 1 variables.");
  }
  auto sorted = conditioning_set;
  std::sort(sorted.begin(), sorted.end());
  if (std::adjacent_find(sorted.begin(), sorted.end()) != sorted.end()) {
    throw std::runtime_error("conditioning_set must not contain duplicates.");
  }
  if ((sorted.front() < 1) || (sorted.back() > d_)) {
    throw std::runtime_error("conditioning_set entries must be in 1, ..., d.");
  }
  if (rvine_structure_.get_trunc_lvl() != d_ - 1) {
    throw std::runtime_error(
      "reorient() requires a non-truncated vine (trunc_lvl == d - 1).");
  }

  std::vector<char> in_b(d_ + 1, 0);
  for (auto v : conditioning_set)
    in_b[v] = 1;

  auto tail_is_b = [&](const std::vector<size_t>& order) {
    for (size_t i = d_ - k; i < d_; ++i)
      if (!in_b[order[i]])
        return false;
    return true;
  };

  // Preserve the exact representation if the requested variables already
  // form the sampling-order tail.
  if (tail_is_b(rvine_structure_.get_order())) {
    ReorientationMap identity;
    identity.identity = true;
    return identity;
  }

  // Peel the structure while steering the diagonal so that the conditioning
  // set lands at the tail. No fitted pair copulas are copied here.
  auto to_tail =
    [&](size_t col,
        const std::vector<std::vector<size_t>>& leaf_edges) noexcept -> size_t {
    bool want_cond = (col >= d_ - k);
    for (const auto& options : leaf_edges)
      for (auto v : options)
        if (want_cond ? (in_b[v] != 0) : (in_b[v] == 0))
          return v;
    return leaf_edges.front().front();
  };

  const std::string inadmissible =
    "conditioning set is not admissible as a sampling-order tail of this "
    "vine; fit with FitControlsVinecop::set_conditioning_set() or condition "
    "on an admissible set.";
  RVineTrees::Decomposition decomposition;
  try {
    decomposition = rvine_structure_.get_trees().to_struct_array_map(to_tail);
  } catch (const std::exception&) {
    throw std::runtime_error(inadmissible);
  }
  if (!tail_is_b(decomposition.order)) {
    throw std::runtime_error(inadmissible);
  }

  RVineStructure structure(
    decomposition.order, decomposition.struct_array, false, true);
  return { std::move(structure),
           std::move(decomposition.pair_copula_locations),
           false };
}

//! @brief Relabels the vine to an equivalent one whose order tail equals
//! `conditioning_set`.
//!
//! @details A value-preserving re-orientation: the fitted model is unchanged
//! (both the density `pdf` and the log-likelihood `loglik` are invariant), only
//! its sampling-order representation changes, so that the variables in
//! `conditioning_set` are drawn first and can be conditioned on with
//! `simulate_conditional`. It chooses, per tree, which conditioned variable
//! sits on the matrix diagonal, and flips each pair copula whose stored
//! orientation no longer matches its new position. Throws if the current
//! structure admits no sampling order ending in `conditioning_set` (fit with a
//! conditioning set via `FitControlsVinecop` to guarantee one exists).
//!
//! @param conditioning_set 1-based variable indices to place at the tail of the
//! variable order.
inline void
Vinecop::reorient(const std::vector<size_t>& conditioning_set)
{
  auto reorientation = make_reorientation_map(conditioning_set);
  if (reorientation.identity)
    return;

  std::vector<std::vector<Bicop>> pair_copulas(d_ - 1);
  for (size_t tree = 0; tree < d_ - 1; ++tree) {
    pair_copulas[tree].reserve(d_ - 1 - tree);
    for (size_t edge = 0; edge < d_ - 1 - tree; ++edge) {
      auto location = reorientation.pair_copulas(tree, edge);
      pair_copulas[tree].emplace_back(
        pair_copulas_[location.tree][location.edge]);
      if (location.flipped)
        pair_copulas[tree].back().flip();
    }
  }

  rvine_structure_ = std::move(reorientation.structure);
  pair_copulas_ = std::move(pair_copulas);
}

//! @brief Fits the parameters of a pre-specified vine copula model.
//!
//! @details This method fits the pair-copulas of a vine copula model. It is
//! assumed that the structure  and pair-copula families are already set.
//! The method is equivalent to calling `fit()` for each pair-copula in the
//! model. The same can be achieved by calling `select()` with the same data
//! and a `FitControlsVinecop` object instantiated
//! with `select_families = false`.
//!
//! @param data \f$ n \times 2d \f$ or \f$ n \times (d + k) \f$ matrix of
//!   observations, where \f$ k \f$ is the number of discrete variables.
//! @param controls The controls for each bivariate fit (see
//! `FitControlsBicop()`).
//! @param num_threads The number of threads to use for parallel computation.
inline void
Vinecop::fit(const Eigen::MatrixXd& data,
             const FitControlsBicop& controls,
             const size_t num_threads)
{
  check_data(data);
  auto u = collapse_data(data);

  // info about the vine structure (reverse rows (!) for more natural indexing)
  size_t trunc_lvl = rvine_structure_.get_trunc_lvl();
  if (trunc_lvl == 0)
    return;

  auto order = rvine_structure_.get_order();
  auto disc_cols = tools_select::get_disc_cols(var_types_);
  size_t n = u.rows();

  // temporary storage objects (all data must be in (0, 1))
  Eigen::MatrixXd hfunc1, hfunc2, hfunc1_sub, hfunc2_sub;
  hfunc1 = Eigen::MatrixXd::Zero(n, d_);
  hfunc2 = Eigen::MatrixXd::Zero(n, d_);
  if (get_n_discrete() > 0) {
    hfunc1_sub = hfunc1;
    hfunc2_sub = hfunc2;
  }

  // set up thread pool
  tools_thread::ThreadPool pool((num_threads == 1) ? 0 : num_threads);

  // fill first row of hfunc2 matrix with observed data;
  // points have to be reordered to correspond to natural order
  for (size_t j = 0; j < d_; ++j) {
    hfunc2.col(j) = u.col(order[j] - 1);
    if (var_types_[order[j] - 1] == "d") {
      hfunc2_sub.col(j) = u.col(d_ + disc_cols[order[j] - 1]);
    }
  }

  for (size_t tree = 0; tree < trunc_lvl; ++tree) {
    tools_interface::check_user_interrupt();
    // scale down the per-fit thread budget: the edges of this tree already
    // run concurrently on the pool, so nested threading would oversubscribe
    FitControlsBicop tree_controls = controls;
    if (num_threads > 1) {
      tree_controls.set_num_threads(
        std::max<size_t>(1, num_threads / std::max<size_t>(1, d_ - tree - 1)));
    }
    auto fit_edge = [&](size_t edge) {
      tools_interface::check_user_interrupt(edge % 5 == 0);
      // extract evaluation point from hfunction matrices (have been
      // computed in previous tree level)
      Bicop* edge_copula = &pair_copulas_[tree][edge];
      auto var_types = edge_copula->get_var_types();
      size_t m = rvine_structure_.min_array(tree, edge);

      Eigen::MatrixXd u_e(n, 2), u_e_sub;
      u_e.col(0) = hfunc2.col(edge);
      if (m == rvine_structure_.struct_array(tree, edge, true)) {
        u_e.col(1) = hfunc2.col(m - 1);
      } else {
        u_e.col(1) = hfunc1.col(m - 1);
      }

      if ((var_types[0] == "d") || (var_types[1] == "d")) {
        u_e.conservativeResize(n, 4);
        u_e.col(2) = hfunc2_sub.col(edge);
        if (m == rvine_structure_.struct_array(tree, edge, true)) {
          u_e.col(3) = hfunc2_sub.col(m - 1);
        } else {
          u_e.col(3) = hfunc1_sub.col(m - 1);
        }
      }

      edge_copula->fit(u_e, tree_controls);

      // h-functions are only evaluated if needed in next tree
      if (rvine_structure_.needed_hfunc1(tree, edge)) {
        hfunc1.col(edge) = edge_copula->hfunc1(u_e);
        if (var_types[1] == "d") {
          u_e_sub = u_e;
          u_e_sub.col(1) = u_e.col(3);
          hfunc1_sub.col(edge) = edge_copula->hfunc1(u_e_sub);
        }
      }
      if (rvine_structure_.needed_hfunc2(tree, edge)) {
        hfunc2.col(edge) = edge_copula->hfunc2(u_e);
        if (var_types[0] == "d") {
          u_e_sub = u_e;
          u_e_sub.col(0) = u_e.col(2);
          hfunc2_sub.col(edge) = edge_copula->hfunc2(u_e_sub);
        }
      }
    };

    pool.map(fit_edge, tools_stl::seq_int(0, d_ - tree - 1));
    pool.wait();
  }
  pool.join();

  loglik_ = 0, nobs_ = n;
  for (size_t tree = 0; tree < trunc_lvl; ++tree) {
    for (size_t edge = 0; edge < d_ - tree - 1; ++edge) {
      loglik_ += pair_copulas_[tree][edge].get_loglik();
    }
  }
}

//! @name Getters and setters
//! @{

//! @brief Gets a pair copula.
//!
//! @param tree Tree index (starting with 0).
//! @param edge Edge index (starting with 0).
inline Bicop
Vinecop::get_pair_copula(const size_t tree, const size_t edge) const
{
  this->check_indices(tree, edge);
  if (tree >= pair_copulas_.size()) {
    return Bicop(); // vine is truncated
  }
  return pair_copulas_[tree][edge];
}

//! @brief Gets all pair copulas.
//!
//! @return a nested std::vector with entry `[t][e]` corresponding to
//! edge `e` in tree `t`.
inline std::vector<std::vector<Bicop>>
Vinecop::get_all_pair_copulas() const
{
  return pair_copulas_;
}

//! @brief Gets the family of a pair copula.
//!
//! @param tree Tree index (starting with 0).
//! @param edge Edge index (starting with 0).
inline BicopFamily
Vinecop::get_family(const size_t tree, const size_t edge) const
{
  this->check_indices(tree, edge);
  if (tree >= pair_copulas_.size()) {
    return BicopFamily::indep; // vine is truncated
  }
  return pair_copulas_[tree][edge].get_family();
}

//! @brief Gets the families of all pair copulas.
//!
//! @return a nested std::vector with entry `[t][e]` corresponding to
//! edge `e` in tree `t`.
inline std::vector<std::vector<BicopFamily>>
Vinecop::get_all_families() const
{
  std::vector<std::vector<BicopFamily>> families(pair_copulas_.size());
  for (size_t tree = 0; tree < pair_copulas_.size(); ++tree) {
    families[tree].resize(d_ - 1 - tree);
    for (size_t edge = 0; edge < d_ - 1 - tree; ++edge) {
      families[tree][edge] = pair_copulas_[tree][edge].get_family();
    }
  }

  return families;
}

//! @brief Gets the rotation of a pair copula.
//!
//! @param tree Tree index (starting with 0).
//! @param edge Edge index (starting with 0).
inline int
Vinecop::get_rotation(const size_t tree, const size_t edge) const
{
  this->check_indices(tree, edge);
  if (tree >= pair_copulas_.size()) {
    return 0; // vine is truncated
  }
  return pair_copulas_[tree][edge].get_rotation();
}

//! @brief Gets the rotations of all pair copulas.
//!
//! @return a nested std::vector with entry `[t][e]` corresponding to
//! edge `e` in tree `t`.
inline std::vector<std::vector<int>>
Vinecop::get_all_rotations() const
{
  std::vector<std::vector<int>> rotations(pair_copulas_.size());
  for (size_t tree = 0; tree < pair_copulas_.size(); ++tree) {
    rotations[tree].resize(d_ - 1 - tree);
    for (size_t edge = 0; edge < d_ - 1 - tree; ++edge) {
      rotations[tree][edge] = pair_copulas_[tree][edge].get_rotation();
    }
  }

  return rotations;
}

//! @brief Gets the parameters of a pair copula.
//!
//! @param tree Tree index (starting with 0).
//! @param edge Edge index (starting with 0).
inline Eigen::MatrixXd
Vinecop::get_parameters(const size_t tree, const size_t edge) const
{
  this->check_indices(tree, edge);
  if (tree >= pair_copulas_.size()) {
    return Eigen::MatrixXd(); // vine is truncated
  }
  return pair_copulas_[tree][edge].get_parameters();
}

//! @brief Gets the Kendall's \f$ tau \f$ of a pair copula.
//!
//! @param tree Tree index (starting with 0).
//! @param edge Edge index (starting with 0).
inline double
Vinecop::get_tau(const size_t tree, const size_t edge) const
{
  this->check_indices(tree, edge);
  if (tree >= pair_copulas_.size()) {
    return 0; // vine is truncated
  }
  return pair_copulas_[tree][edge].get_tau();
}

inline size_t
Vinecop::get_trunc_lvl() const
{
  return rvine_structure_.get_trunc_lvl();
}

//! @brief Gets the parameters of all pair copulas.
//!
//! @return a nested std::vector with entry `[t][e]` corresponding to
//! edge `e` in tree `t`.
inline std::vector<std::vector<Eigen::MatrixXd>>
Vinecop::get_all_parameters() const
{
  std::vector<std::vector<Eigen::MatrixXd>> parameters(pair_copulas_.size());
  for (size_t tree = 0; tree < parameters.size(); ++tree) {
    parameters[tree].resize(d_ - 1 - tree);
    for (size_t edge = 0; edge < d_ - 1 - tree; ++edge) {
      parameters[tree][edge] = pair_copulas_[tree][edge].get_parameters();
    }
  }

  return parameters;
}

//! @brief Gets the Kendall's \f$ tau \f$s of all pair copulas.
//!
//! @return a nested std::vector with entry `[t][e]` corresponding to
//! edge `e` in tree `t`.
inline std::vector<std::vector<double>>
Vinecop::get_all_taus() const
{
  std::vector<std::vector<double>> taus(pair_copulas_.size());
  for (size_t tree = 0; tree < taus.size(); ++tree) {
    taus[tree].resize(d_ - 1 - tree);
    for (size_t edge = 0; edge < d_ - 1 - tree; ++edge) {
      taus[tree][edge] = pair_copulas_[tree][edge].get_tau();
    }
  }

  return taus;
}

//! @brief Gets the dimension of the vine copula model.
inline size_t
Vinecop::get_dim() const
{
  return d_;
}

//! @brief Gets the order vector of the vine copula model.
inline std::vector<size_t>
Vinecop::get_order() const
{
  return rvine_structure_.get_order();
}

//! @brief Gets the structure matrix of the vine copula model.
inline RVineStructure
Vinecop::get_rvine_structure() const
{
  return rvine_structure_;
}

//! @brief Gets the structure matrix of the vine copula model.
inline Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>
Vinecop::get_matrix() const
{
  return rvine_structure_.get_matrix();
}

//! @brief Gets the above diagonal coefficients of the vine copula model.
inline TriangularArray<size_t>
Vinecop::get_struct_array(bool natural_order) const
{
  return rvine_structure_.get_struct_array(natural_order);
}

//! @brief Gets the vine as a list-of-trees decomposition, each edge carrying
//! its fitted pair-copula (stored with its first argument on the diagonal).
inline RVineTrees
Vinecop::get_trees() const
{
  // decompose in original labels: the diagonal `get_order()`, the original-
  // label structure array, and the pair-copulas share the same labelling.
  return RVineTrees(rvine_structure_.get_order(),
                    rvine_structure_.get_struct_array(false),
                    pair_copulas_);
}

//! @brief Gets the log-likelihood (throws an error if model has not been.
//! fitted to data).
inline double
Vinecop::get_loglik() const
{
  check_fitted();
  return loglik_;
}

//! @brief Gets the number of observations used for the fit.
//!
//! The function throws an error if model has not been fitted to data.
inline size_t
Vinecop::get_nobs() const
{
  check_fitted();
  return nobs_;
}

//! @brief Gets the AIC.
//!
//! The function throws an error if model has not been fitted to data.
inline double
Vinecop::get_aic() const
{
  check_fitted();
  return -2 * loglik_ + 2 * get_npars();
}

//! @brief Gets the BIC.
//!
//! The function throws an error if model has not been fitted to data.
inline double
Vinecop::get_bic() const
{
  check_fitted();
  return -2 * loglik_ + get_npars() * std::log(nobs_);
}

//! @brief Gets the log-likelihood.
//!
//! The function throws an error if model has not been fitted to data.
inline double
Vinecop::get_mbicv(const double psi0) const
{
  check_fitted();
  return -2 * loglik_ + this->calculate_mbicv_penalty(nobs_, psi0);
}

//! @brief Computes the penalty term for mBICV.
inline double
Vinecop::calculate_mbicv_penalty(const size_t nobs, const double psi0) const
{
  if ((psi0 <= 0.0) || (psi0 >= 1.0)) {
    throw std::runtime_error("psi0 must be in the interval (0, 1)");
  }
  auto all_fams = get_all_families();
  Eigen::Matrix<size_t, Eigen::Dynamic, 1> non_indeps(d_ - 1);
  non_indeps.setZero();
  for (size_t t = 0; t < d_ - 1; t++) {
    if (t == all_fams.size()) {
      break;
    }
    for (size_t e = 0; e < d_ - 1 - t; e++) {
      if (all_fams[t][e] != BicopFamily::indep) {
        non_indeps(t)++;
      }
    }
  }
  auto sq0 = tools_stl::seq_int(1, d_ - 1);
  Eigen::Matrix<size_t, Eigen::Dynamic, 1> sq(d_ - 1);
  auto psis = Eigen::VectorXd(d_ - 1);
  for (size_t i = 0; i < d_ - 1; i++) {
    sq(i) = sq0[i];
    psis(i) = std::pow(psi0, sq0[i]);
  }
  double npars = this->get_npars();
  double log_prior = (non_indeps.cast<double>().array() * psis.array().log() +
                      (d_ - non_indeps.array() - sq.array()).cast<double>() *
                        (1 - psis.array()).log())
                       .sum();

  return std::log(nobs) * npars - 2.0 * log_prior;
}

//! @brief Gets the threshold.
//!
//! Usually zero except `select_threshold == TRUE` in `FitControlsVinecop()`).
inline double
Vinecop::get_threshold() const
{
  return threshold_;
}

//! @brief Sets variable types.
inline void
Vinecop::set_var_types(const std::vector<std::string>& var_types)
{
  check_var_types(var_types);
  set_var_types_internal(var_types);
}

//! @brief Sets all pair-copulas.
inline void
Vinecop::set_all_pair_copulas(
  const std::vector<std::vector<Bicop>>& pair_copulas)
{
  check_pair_copulas_rvine_structure(pair_copulas);
  pair_copulas_ = pair_copulas;
  rvine_structure_.truncate(pair_copulas.size());
}

inline void
Vinecop::check_var_types(const std::vector<std::string>& var_types) const
{
  std::stringstream msg;
  if (var_types.size() > d_) {
    msg << "more var_types (" << var_types.size() << ") "
        << "than variables (" << d_ << ")" << std::endl;
    throw std::runtime_error(msg.str());
  }
  for (const auto& t : var_types) {
    if (!tools_stl::is_member(t, { "c", "d" })) {
      msg << "variable type must be 'c' or 'd' (not '" << t << "')."
          << std::endl;
      throw std::runtime_error(msg.str());
    }
  }
}

//! @brief Sets variable types.
//! @param var_types A vector specifying the types of the variables,
//!   e.g., `{"c", "d"}` means first varible continuous, second discrete.
inline void
Vinecop::set_var_types_internal(const std::vector<std::string>& var_types)
{
  var_types_ = var_types;
  n_discrete_ = 0;
  for (const auto& t : var_types_) {
    n_discrete_ += (t == "d");
  }
  if (pair_copulas_.empty()) {
    return;
  }

  // set new var_types for all pair-copulas
  const auto order = rvine_structure_.get_order();
  std::vector<std::string> natural_types(d_), pair_types(2);
  for (size_t j = 0; j < d_; ++j) {
    natural_types[j] = var_types[order[j] - 1];
  }
  // we set the first tree explicitly and deduce later trees
  for (size_t e = 0; e < d_ - 1; ++e) {
    pair_types[0] = natural_types[e];
    pair_types[1] =
      natural_types[rvine_structure_.struct_array(0, e, true) - 1];
    pair_copulas_[0][e].set_var_types(pair_types);
  }

  for (size_t t = 1; t < pair_copulas_.size(); ++t) {
    for (size_t e = 0; e < d_ - t - 1; ++e) {
      size_t m = rvine_structure_.min_array(t, e);
      pair_types[0] = pair_copulas_[t - 1][e].get_var_types()[0];
      if (m == rvine_structure_.struct_array(t, e, true)) {
        pair_types[1] = pair_copulas_[t - 1][m - 1].get_var_types()[0];
      } else {
        pair_types[1] = pair_copulas_[t - 1][m - 1].get_var_types()[1];
      }
      pair_copulas_[t][e].set_var_types(pair_types);
    }
  }
}

//! @brief Gets the variable types.
inline std::vector<std::string>
Vinecop::get_var_types() const
{
  return var_types_;
}

//! @}

//! @brief Evaluates the per-pair copula density and h-functions.
//!
//! @details The copula density is defined as joint density divided by marginal
//! densities, irrespective of variable types.
//!
//! When at least one variable is discrete, two types of
//! "observations" are required in `u`: the first \f$ n \; x \; d \f$ block
//! contains realizations of \f$ F_{X_j}(X_j) \f$.
//! The second \f$ n \; x \; d \f$
//! block contains realizations of \f$ F_{X_j}(X_j^-) \f$. The minus indicates a
//! left-sided limit of the cdf. For, e.g., an integer-valued variable, it holds
//! \f$ F_{X_j}(X_j^-) = F_{X_j}(X_j - 1) \f$. For continuous variables the left
//! limit and the cdf itself coincide. Respective columns can be omitted in the
//! second block.
//!
//! @param u An \f$ n \times 2d \f$ or \f$ n \times (d + k) \f$ matrix of
//!   evaluation points, where \f$ k \f$ is the number of discrete variables
//!   (see `Vinecop::select()`).
//! @param num_threads The number of threads to use for computations; if greater
//!   than 1, the function will be applied concurrently to `num_threads` batches
//!   of `u`.
//! @param keep_all Whether to keep and return per-edge pdfs and h-functions.
//! @return A struct containing:
//!   - `pdf`: the copula density evaluated at `u`.
//! If `keep_all = true`, the struct also contains the following fields:
//!   - `pdf_edges`: a triangular array of vectors containing
//!     the per-edge copula densities evaluated at `u`.
//!   - `hfunc1`:  a triangular array of vectors containing the first h-function
//!   of each edge evaluated at `u`.
//!   - `hfunc2`:  a triangular array of vectors containing the second
//!   h-function of each edge evaluated at `u`.
//!   - `hfunc1_sub`: a triangular array of vectors containing the first
//!   h-function of each edge evaluated at the second block of `u` (i.e., the
//!   left-sided limits), if at least one variable is discrete.
//!   - `hfunc2_sub`: a triangular array of vectors containing the second
//!   h-function of each edge evaluated at the second block of `u` (i.e., the
//!   left-sided limits), if at least one variable is discrete.
inline Vinecop::PdfWithHfuncsResult
Vinecop::pdf_full(Eigen::MatrixXd u,
                  const size_t num_threads,
                  const bool keep_all) const
{
  return pdf_full(std::move(u), Eigen::MatrixXd(), num_threads, keep_all);
}

//! @brief Evaluates the copula density (and per-edge quantities) with
//! per-observation parameters.
//!
//! Same as `pdf_full()`, but each observation uses its own full-vine parameter
//! vector, supplied as an \f$ n \times \mathrm{npars} \f$ matrix `parameters`
//! whose columns follow the (tree, edge, parameter) order of `scores()`.
//! Continuous, all-parametric models only (discrete variables and nonparametric
//! pair copulas are rejected).
inline Vinecop::PdfWithHfuncsResult
Vinecop::pdf_full(Eigen::MatrixXd u,
                  const Eigen::MatrixXd& parameters,
                  const size_t num_threads,
                  const bool keep_all) const
{
  check_data(u);
  collapse_data_inplace(u);

  const bool per_obs = parameters.size() > 0;
  if (per_obs) {
    check_parametric("pdf()/loglik() with per-observation parameters");
    check_per_obs_params(u, parameters);
  }

  // info about the vine structure (reverse rows (!) for more natural indexing)
  size_t trunc_lvl = rvine_structure_.get_trunc_lvl();
  auto order = rvine_structure_.get_order();
  auto disc_cols = tools_select::get_disc_cols(var_types_);
  const bool discrete = is_discrete();

  PdfWithHfuncsResult result;

  if (keep_all) {
    result.pdf_edges = TriangularArray<Eigen::VectorXd>(d_, trunc_lvl);
    result.hfunc1 = TriangularArray<Eigen::VectorXd>(d_, trunc_lvl);
    result.hfunc2 = TriangularArray<Eigen::VectorXd>(d_, trunc_lvl);
    if (discrete) {
      result.hfunc1_sub = TriangularArray<Eigen::VectorXd>(d_, trunc_lvl);
      result.hfunc2_sub = TriangularArray<Eigen::VectorXd>(d_, trunc_lvl);
    }
    for (size_t tree = 0; tree < trunc_lvl; ++tree) {
      for (size_t edge = 0; edge < d_ - tree - 1; ++edge) {
        result.pdf_edges(tree, edge) = Eigen::VectorXd::Zero(u.rows());
        if (rvine_structure_.needed_hfunc1(tree, edge)) {
          result.hfunc1(tree, edge) = Eigen::VectorXd::Zero(u.rows());
          if (discrete) {
            result.hfunc1_sub(tree, edge) = Eigen::VectorXd::Zero(u.rows());
          }
        }
        if (rvine_structure_.needed_hfunc2(tree, edge)) {
          result.hfunc2(tree, edge) = Eigen::VectorXd::Zero(u.rows());
          if (discrete) {
            result.hfunc2_sub(tree, edge) = Eigen::VectorXd::Zero(u.rows());
          }
        }
      }
    }
  }

  // initial value must be 1.0 for multiplication
  result.pdf = Eigen::VectorXd::Constant(u.rows(), 1.0);

  auto do_batch = [&](const tools_batch::Batch& b) {
    // temporary storage objects (all data must be in (0, 1))
    Eigen::MatrixXd hfunc1, hfunc2, u_e, hfunc1_sub, hfunc2_sub, u_e_sub;
    hfunc1 = Eigen::MatrixXd::Zero(b.size, d_);
    hfunc2 = Eigen::MatrixXd::Zero(b.size, d_);
    if (discrete) {
      hfunc1_sub = hfunc1;
      hfunc2_sub = hfunc2;
    }

    // fill first row of hfunc2 matrix with evaluation points;
    // points have to be reordered to correspond to natural order
    for (size_t j = 0; j < d_; ++j) {
      hfunc2.col(j) = u.block(b.begin, order[j] - 1, b.size, 1);
      if (var_types_[order[j] - 1] == "d") {
        hfunc2_sub.col(j) =
          u.block(b.begin, d_ + disc_cols[order[j] - 1], b.size, 1);
      }
    }

    // running column offset of the current edge's parameters within the flat
    // (tree, edge, parameter) order of `parameters` (per-observation path only)
    size_t par_offset = 0;

    for (size_t tree = 0; tree < trunc_lvl; ++tree) {
      tools_interface::check_user_interrupt(
        static_cast<double>(u.rows()) * static_cast<double>(d_) > 1e5);
      for (size_t edge = 0; edge < d_ - tree - 1; ++edge) {
        tools_interface::check_user_interrupt(edge % 100 == 0);
        // extract evaluation point from hfunction matrices (have been
        // computed in previous tree level)
        Bicop* edge_copula = &pair_copulas_[tree][edge];
        auto var_types = edge_copula->get_var_types();
        size_t m = rvine_structure_.min_array(tree, edge);

        u_e.resize(b.size, 2);
        u_e.col(0) = hfunc2.col(edge);
        if (m == rvine_structure_.struct_array(tree, edge, true)) {
          u_e.col(1) = hfunc2.col(m - 1);
        } else {
          u_e.col(1) = hfunc1.col(m - 1);
        }

        if ((var_types[0] == "d") || (var_types[1] == "d")) {
          u_e.conservativeResize(b.size, 4);
          u_e.col(2) = hfunc2_sub.col(edge);
          if (m == rvine_structure_.struct_array(tree, edge, true)) {
            u_e.col(3) = hfunc2_sub.col(m - 1);
          } else {
            u_e.col(3) = hfunc1_sub.col(m - 1);
          }
        }

        // per-edge dispatch: on the per-observation path each helper slices
        // this edge's n x np parameter block and routes to the per-row Bicop
        // overload; on the fixed path it calls the stored-parameter method
        // unchanged. Discrete never arises here (per-observation requires
        // continuous), so the _sub h-functions stay on the scalar path.
        Eigen::MatrixXd pars_e;
        if (per_obs) {
          size_t np = static_cast<size_t>(edge_copula->get_parameters().size());
          pars_e = parameters.block(b.begin, par_offset, b.size, np);
          par_offset += np;
        }
        auto ec_pdf = [&]() {
          return per_obs ? edge_copula->pdf(u_e, pars_e)
                         : edge_copula->pdf(u_e);
        };
        auto ec_hfunc1 = [&]() {
          return per_obs ? edge_copula->hfunc1(u_e, pars_e)
                         : edge_copula->hfunc1(u_e);
        };
        auto ec_hfunc2 = [&]() {
          return per_obs ? edge_copula->hfunc2(u_e, pars_e)
                         : edge_copula->hfunc2(u_e);
        };

        Eigen::VectorXd edge_pdf = ec_pdf();
        result.pdf.segment(b.begin, b.size) =
          result.pdf.segment(b.begin, b.size).cwiseProduct(edge_pdf);

        // h-functions are only evaluated if needed in next step
        if (rvine_structure_.needed_hfunc1(tree, edge)) {
          hfunc1.col(edge) = ec_hfunc1();
          if (var_types[1] == "d") {
            u_e_sub = u_e;
            u_e_sub.col(1) = u_e.col(3);
            hfunc1_sub.col(edge) = edge_copula->hfunc1(u_e_sub);
          }
        }
        if (rvine_structure_.needed_hfunc2(tree, edge)) {
          hfunc2.col(edge) = ec_hfunc2();
          if (var_types[0] == "d") {
            u_e_sub = u_e;
            u_e_sub.col(0) = u_e.col(2);
            hfunc2_sub.col(edge) = edge_copula->hfunc2(u_e_sub);
          }
        }

        if (keep_all) {
          result.pdf_edges(tree, edge).segment(b.begin, b.size) = edge_pdf;
          if (rvine_structure_.needed_hfunc1(tree, edge)) {
            result.hfunc1(tree, edge).segment(b.begin, b.size) =
              hfunc1.col(edge);
            if (discrete) {
              result.hfunc1_sub(tree, edge).segment(b.begin, b.size) =
                hfunc1_sub.col(edge);
            }
          }
          if (rvine_structure_.needed_hfunc2(tree, edge)) {
            result.hfunc2(tree, edge).segment(b.begin, b.size) =
              hfunc2.col(edge);
            if (discrete) {
              result.hfunc2_sub(tree, edge).segment(b.begin, b.size) =
                hfunc2_sub.col(edge);
            }
          }
        }
      }
    }
  };

  if (trunc_lvl > 0) {
    tools_thread::ThreadPool pool((num_threads == 1) ? 0 : num_threads);
    pool.map(do_batch, tools_batch::create_batches(u.rows(), num_threads));
    pool.join();
  }

  return result;
}

//! @brief Evaluates the copula density.
//!
//! @details The copula density is defined as joint density divided by marginal
//! densities, irrespective of variable types.
//!
//! When at least one variable is discrete, two types of
//! "observations" are required in `u`: the first \f$ n \; x \; d \f$ block
//! contains realizations of \f$ F_{X_j}(X_j) \f$.
//! The second \f$ n \; x \; d \f$
//! block contains realizations of \f$ F_{X_j}(X_j^-) \f$. The minus indicates a
//! left-sided limit of the cdf. For, e.g., an integer-valued variable, it holds
//! \f$ F_{X_j}(X_j^-) = F_{X_j}(X_j - 1) \f$. For continuous variables the left
//! limit and the cdf itself coincide. Respective columns can be omitted in the
//! second block.
//!
//! @param u An \f$ n \times 2d \f$ or \f$ n \times (d + k) \f$ matrix of
//!   evaluation points, where \f$ k \f$ is the number of discrete variables
//!   (see `Vinecop::select()`).
//! @param num_threads The number of threads to use for computations; if greater
//!   than 1, the function will be applied concurrently to `num_threads` batches
//!   of `u`.
//! @return A vector of length `n` containing the copula density values.
inline Eigen::VectorXd
Vinecop::pdf(Eigen::MatrixXd u, const size_t num_threads) const
{
  return pdf_full(std::move(u), num_threads, false).pdf;
}

//! @brief Evaluates the copula density with per-observation parameters.
//!
//! Per-observation counterpart of `pdf()`; see the per-observation
//! `pdf_full()` overload for the `parameters` layout and restrictions.
inline Eigen::VectorXd
Vinecop::pdf(Eigen::MatrixXd u,
             const Eigen::MatrixXd& parameters,
             const size_t num_threads) const
{
  return pdf_full(std::move(u), parameters, num_threads, false).pdf;
}

//! throws if the model has a nonparametric pair copula (see scores()).
inline void
Vinecop::check_parametric(const char* fn) const
{
  size_t trunc_lvl = rvine_structure_.get_trunc_lvl();
  for (size_t t = 0; t < trunc_lvl; t++) {
    for (size_t e = 0; e < d_ - 1 - t; e++) {
      if (!tools_stl::is_member(pair_copulas_[t][e].get_family(),
                                bicop_families::parametric)) {
        throw std::runtime_error(
          std::string(fn) +
          " is only available for models with parametric pair copulas");
      }
    }
  }
}

//! validates a per-observation parameter matrix: the analytic cascade requires
//! continuous variables (discrete would need finite differences that mutate
//! shared parameters), and `parameters` must be n x npars, one full-vine
//! parameter vector per observation in the (tree, edge, param) order of
//! scores().
inline void
Vinecop::check_per_obs_params(const Eigen::MatrixXd& u,
                              const Eigen::MatrixXd& per_obs_params) const
{
  if (get_n_discrete() > 0) {
    throw std::runtime_error(
      "per-observation parameters are only supported for continuous "
      "variables");
  }
  if (per_obs_params.rows() != u.rows()) {
    throw std::runtime_error("parameters must have one row per row of u "
                             "(parameters.rows() must equal u.rows()).");
  }
  if (per_obs_params.cols() != static_cast<Eigen::Index>(this->get_npars())) {
    throw std::runtime_error(
      "parameters must have get_npars() columns, in the (tree, edge, "
      "parameter) order of scores().");
  }
}

//! first-order chain-rule push-forward of an argument perturbation through
//! an h-function output: `∂h/∂u1 · v1 + ∂h/∂u2 · v2` (shared by the gradient
//! cascade and the Hessian's hat/til propagations).
inline Eigen::VectorXd
propagate_first_order(const Eigen::VectorXd& du1,
                      const Eigen::VectorXd& du2,
                      const Eigen::VectorXd& v1,
                      const Eigen::VectorXd& v2)
{
  return du1.cwiseProduct(v1) + du2.cwiseProduct(v2);
}

//! @brief Builds the per-edge derivative caches by one forward walk over the
//! rows `[begin, begin + size)` of `u`.
//!
//! This is the single forward pass shared by the analytic gradient, joint
//! Hessian, and step-wise Hessian: it assembles each edge's arguments `u_e`
//! from the previous tree's h-functions and evaluates the pair copulas'
//! derivative leaves. Continuous, all-parametric models only (the callers
//! reject nonparametric families and route discrete/step-wise=true to finite
//! differences). Second-order leaves are computed only when `second_order`.
inline TriangularArray<Vinecop::DerivCache>
Vinecop::build_deriv_cache(const Eigen::MatrixXd& u,
                           size_t begin,
                           size_t size,
                           bool second_order,
                           const Eigen::MatrixXd& per_obs_params) const
{
  size_t trunc_lvl = rvine_structure_.get_trunc_lvl();
  auto order = rvine_structure_.get_order();
  const size_t m = size;
  TriangularArray<DerivCache> cache(d_, trunc_lvl);

  Eigen::MatrixXd hfunc1 = Eigen::MatrixXd::Zero(m, d_);
  Eigen::MatrixXd hfunc2 = Eigen::MatrixXd::Zero(m, d_);
  for (size_t j = 0; j < d_; ++j) {
    hfunc2.col(j) = u.block(begin, order[j] - 1, m, 1);
  }
  auto sel = [](size_t p) { return "par" + std::to_string(p + 1); };

  // per-observation parameters: when supplied, each edge reads its own n x np
  // column block (in the same (t, e) order the flat scores() index follows);
  // `par_offset` tracks that column start. Empty => fixed-parameter fast path.
  const bool per_obs = per_obs_params.size() > 0;
  size_t par_offset = 0;

  Eigen::MatrixXd u_e;
  for (size_t t = 0; t < trunc_lvl; ++t) {
    tools_interface::check_user_interrupt(
      static_cast<double>(u.rows()) * static_cast<double>(d_) > 1e3);
    for (size_t e = 0; e < d_ - 1 - t; ++e) {
      tools_interface::check_user_interrupt(e % 100 == 0);
      const Bicop& ec = pair_copulas_[t][e]; // const ref: no per-edge copy
      DerivCache& ce = cache(t, e);
      size_t np = static_cast<size_t>(ec.get_parameters().size());
      ce.np = np;
      ce.arg2_col = rvine_structure_.min_array(t, e);
      ce.arg2_is_h2 =
        (ce.arg2_col == rvine_structure_.struct_array(t, e, true));

      u_e = Eigen::MatrixXd(m, 2);
      u_e.col(0) = hfunc2.col(e);
      u_e.col(1) = ce.arg2_is_h2 ? hfunc2.col(ce.arg2_col - 1)
                                 : hfunc1.col(ce.arg2_col - 1);

      // per-edge dispatch: on the fixed path (per_obs == false) each helper
      // calls the pair copula's stored-parameter method exactly as before; on
      // the per-observation path it forwards this edge's n x np parameter slice
      // to the matching per-row overload. Every evaluation (density,
      // h-functions, and their derivatives) must use the slice so the deeper
      // trees see pseudo-observations consistent with the per-obs parameters.
      Eigen::MatrixXd pars_e;
      if (per_obs) {
        pars_e = per_obs_params.block(begin, par_offset, m, np);
      }
      auto ec_pdf = [&]() {
        return per_obs ? ec.pdf(u_e, pars_e) : ec.pdf(u_e);
      };
      auto ec_logpdf_deriv = [&](const std::string& s) {
        return per_obs ? ec.logpdf_deriv(u_e, s, pars_e)
                       : ec.logpdf_deriv(u_e, s);
      };
      auto ec_pdf_deriv = [&](const std::string& s) {
        return per_obs ? ec.pdf_deriv(u_e, s, pars_e) : ec.pdf_deriv(u_e, s);
      };
      auto ec_logpdf_deriv2 = [&](const std::string& s) {
        return per_obs ? ec.logpdf_deriv2(u_e, s, pars_e)
                       : ec.logpdf_deriv2(u_e, s);
      };
      auto ec_hfunc1 = [&]() {
        return per_obs ? ec.hfunc1(u_e, pars_e) : ec.hfunc1(u_e);
      };
      auto ec_hfunc2 = [&]() {
        return per_obs ? ec.hfunc2(u_e, pars_e) : ec.hfunc2(u_e);
      };
      auto ec_hfunc1_deriv = [&](const std::string& s) {
        return per_obs ? ec.hfunc1_deriv(u_e, s, pars_e)
                       : ec.hfunc1_deriv(u_e, s);
      };
      auto ec_hfunc2_deriv = [&](const std::string& s) {
        return per_obs ? ec.hfunc2_deriv(u_e, s, pars_e)
                       : ec.hfunc2_deriv(u_e, s);
      };
      auto ec_hfunc1_deriv2 = [&](const std::string& s) {
        return per_obs ? ec.hfunc1_deriv2(u_e, s, pars_e)
                       : ec.hfunc1_deriv2(u_e, s);
      };
      auto ec_hfunc2_deriv2 = [&](const std::string& s) {
        return per_obs ? ec.hfunc2_deriv2(u_e, s, pars_e)
                       : ec.hfunc2_deriv2(u_e, s);
      };

      ce.c = ec_pdf();
      ce.du1 = ec_logpdf_deriv("u1");
      ce.du2 = ec_logpdf_deriv("u2");
      ce.dpar.resize(np);
      for (size_t p = 0; p < np; ++p) {
        ce.dpar[p] = ec_logpdf_deriv(sel(p));
      }

      // ∂c/∂u1, ∂c/∂u2, ∂c/∂θ are shared by the two h-outputs' 2nd-order
      // fields (via ∂h2/∂u1 = ∂h1/∂u2 = c); only needed for second order
      Eigen::VectorXd c_u1, c_u2;
      std::vector<Eigen::VectorXd> c_par;
      if (second_order) {
        c_u1 = ec_pdf_deriv("u1");
        c_u2 = ec_pdf_deriv("u2");
        c_par.resize(np);
        for (size_t p = 0; p < np; ++p) {
          c_par[p] = ec_pdf_deriv(sel(p));
        }
        ce.du1u1 = ec_logpdf_deriv2("u1u1");
        ce.du1u2 = ec_logpdf_deriv2("u1u2");
        ce.du2u2 = ec_logpdf_deriv2("u2u2");
        ce.dpar_u1.resize(np);
        ce.dpar_u2.resize(np);
        ce.dpar_par.assign(np, std::vector<Eigen::VectorXd>(np));
        for (size_t p = 0; p < np; ++p) {
          ce.dpar_u1[p] = ec_logpdf_deriv2(sel(p) + "u1");
          ce.dpar_u2[p] = ec_logpdf_deriv2(sel(p) + "u2");
          for (size_t q = p; q < np; ++q) {
            ce.dpar_par[p][q] = ec_logpdf_deriv2(sel(p) + sel(q));
            ce.dpar_par[q][p] = ce.dpar_par[p][q];
          }
        }
      }

      bool has_deeper_tree = (t + 1 < trunc_lvl);
      // hfunc2 output: ∂h2/∂u1 = c (identity), ∂h2/∂u2 conditioning
      if (has_deeper_tree && rvine_structure_.needed_hfunc2(t, e)) {
        DerivLeaf& leaf = ce.h2;
        leaf.active = true;
        leaf.du1 = ce.c;
        leaf.du2 = ec_hfunc2_deriv("u2");
        leaf.dpar.resize(np);
        for (size_t p = 0; p < np; ++p) {
          leaf.dpar[p] = ec_hfunc2_deriv(sel(p));
        }
        if (second_order) {
          leaf.du1u1 = c_u1;
          leaf.du1u2 = c_u2;
          leaf.du2u2 = ec_hfunc2_deriv2("u2u2");
          leaf.dpar_u1.resize(np);
          leaf.dpar_u2.resize(np);
          leaf.dpar_par.assign(np, std::vector<Eigen::VectorXd>(np));
          for (size_t p = 0; p < np; ++p) {
            leaf.dpar_u1[p] = c_par[p]; // ∂²h2/∂θ∂u1 = ∂c/∂θ
            leaf.dpar_u2[p] = ec_hfunc2_deriv2(sel(p) + "u2");
            for (size_t q = p; q < np; ++q) {
              leaf.dpar_par[p][q] = ec_hfunc2_deriv2(sel(p) + sel(q));
              leaf.dpar_par[q][p] = leaf.dpar_par[p][q];
            }
          }
        }
      }
      // hfunc1 output: ∂h1/∂u2 = c (identity), ∂h1/∂u1 conditioning
      if (has_deeper_tree && rvine_structure_.needed_hfunc1(t, e)) {
        DerivLeaf& leaf = ce.h1;
        leaf.active = true;
        leaf.du1 = ec_hfunc1_deriv("u1");
        leaf.du2 = ce.c;
        leaf.dpar.resize(np);
        for (size_t p = 0; p < np; ++p) {
          leaf.dpar[p] = ec_hfunc1_deriv(sel(p));
        }
        if (second_order) {
          leaf.du1u1 = ec_hfunc1_deriv2("u1u1");
          leaf.du1u2 = c_u1; // ∂²h1/∂u1∂u2 = ∂c/∂u1
          leaf.du2u2 = c_u2; // ∂²h1/∂u2²   = ∂c/∂u2
          leaf.dpar_u1.resize(np);
          leaf.dpar_u2.resize(np);
          leaf.dpar_par.assign(np, std::vector<Eigen::VectorXd>(np));
          for (size_t p = 0; p < np; ++p) {
            leaf.dpar_u1[p] = ec_hfunc1_deriv2(sel(p) + "u1");
            leaf.dpar_u2[p] = c_par[p]; // ∂²h1/∂θ∂u2 = ∂c/∂θ
            for (size_t q = p; q < np; ++q) {
              leaf.dpar_par[p][q] = ec_hfunc1_deriv2(sel(p) + sel(q));
              leaf.dpar_par[q][p] = leaf.dpar_par[p][q];
            }
          }
        }
      }

      if (rvine_structure_.needed_hfunc1(t, e)) {
        hfunc1.col(e) = ec_hfunc1();
      }
      if (rvine_structure_.needed_hfunc2(t, e)) {
        hfunc2.col(e) = ec_hfunc2();
      }

      par_offset += np;
    }
  }
  return cache;
}

//! @brief Evaluates the score function together with the per-edge
//! derivative caches.
//!
//! The score function is defined as the gradient of the log-likelihood
//! with respect to the parameters.
//!
//! The scores are computed analytically from the pair copulas' derivatives
//! (`Bicop::logpdf_deriv()` and, for the full gradient, the h-function
//! derivatives propagated through the vine by the chain rule); models with
//! discrete variables use central finite differences instead (in which case
//! the caches below stay empty). Models with nonparametric pair copulas are
//! rejected (differentiating w.r.t. the interpolation grid is meaningless).
//!
//! @param u An \f$ n \times 2d \f$ or \f$ n \times (d + k) \f$ matrix of
//!   evaluation points, where \f$ k \f$ is the number of discrete variables
//!   (see `select()`).
//! @param step_wise if `false`, full gradient of the log-likelihood; if `true`,
//!   score function of the step-wise MLE (gradients computed per pair-copula,
//!   treating the pseudo-observations as fixed).
//! @param num_threads The number of threads to use for computations; if greater
//!   than 1, the function will be applied concurrently to `num_threads` batches
//!   of `u`.
//! @param keep_all Whether to keep and return the per-edge derivative caches.
//! @return A struct containing:
//!   - `scores`: the \f$ n \times \mathrm{npars} \f$ matrix of
//!     per-observation scores.
//! If `keep_all = true` and the analytic path applies, the struct also
//! contains triangular arrays of per-edge quantities (indexed by tree and
//! edge, evaluated at the edge's pseudo-observations):
//!   - `pdf_edges`: the edge densities (full gradient only);
//!   - `logpdf_deriv_pars`: the log-density derivative w.r.t. each of the
//!     edge's parameters (= the step-wise scores);
//!   - `hfunc1_deriv_pars`, `hfunc2_deriv_pars`: the h-function derivatives
//!     w.r.t. each parameter, where a deeper tree consumes them (full
//!     gradient only);
//!   - `logpdf_deriv_u1`, `logpdf_deriv_u2`: the log-density derivatives
//!     w.r.t. the edge's arguments (full gradient only, trees > 0);
//!   - `hfunc1_deriv_u1`, `hfunc2_deriv_u2`: the h-function derivatives
//!     w.r.t. their conditioning argument, where a deeper tree consumes
//!     them (full gradient only).
//!
//! @literature
//! Stoeber, J. and Schepsmeier, U. (2013). Estimating standard errors in
//! regular vine copula models. Computational Statistics, 28 (6), 2679-2707.
inline Vinecop::ScoresResult
Vinecop::scores_full(Eigen::MatrixXd u,
                     bool step_wise,
                     const size_t num_threads,
                     const bool keep_all)
{
  return scores_full(
    std::move(u), Eigen::MatrixXd(), step_wise, num_threads, keep_all);
}

//! @brief Evaluates the score function with per-observation parameters.
//!
//! Same as `scores_full()`, but each observation uses its own full-vine
//! parameter vector, supplied as an \f$ n \times \mathrm{npars} \f$ matrix
//! `parameters` whose columns follow the (tree, edge, parameter) order of
//! `scores()`. Continuous, all-parametric models only; discrete variables are
//! rejected (they would require finite differences that mutate shared
//! parameters).
inline Vinecop::ScoresResult
Vinecop::scores_full(Eigen::MatrixXd u,
                     const Eigen::MatrixXd& per_obs_params,
                     bool step_wise,
                     const size_t num_threads,
                     const bool keep_all)
{
  check_data(u);
  u = collapse_data(u);

  // info about the vine structure (reverse rows (!) for more natural indexing)
  size_t trunc_lvl = rvine_structure_.get_trunc_lvl();
  auto order = rvine_structure_.get_order();
  auto disc_cols = tools_select::get_disc_cols(var_types_);
  const size_t n = static_cast<size_t>(u.rows());

  check_parametric("scores()");

  const bool per_obs = per_obs_params.size() > 0;
  if (per_obs) {
    check_per_obs_params(u, per_obs_params);
  }

  ScoresResult result;
  result.scores =
    Eigen::MatrixXd::Zero(n, static_cast<size_t>(this->get_npars()));
  if (trunc_lvl == 0) {
    return result;
  }

  if (!step_wise) {
    // the analytic cascade requires continuous variables (nonparametric
    // families are already rejected above); discrete variables keep the
    // finite-difference path
    bool analytic_ok = (get_n_discrete() == 0);

    if (!analytic_ok) {
      // fall back to finite differences of the whole-vine log-density
      size_t ipar = 0;
      for (size_t t = 0; t < trunc_lvl; t++) {
        for (size_t e = 0; e < d_ - 1 - t; e++) {
          auto pars = pair_copulas_[t][e].get_parameters();
          auto ub = pair_copulas_[t][e].get_parameters_upper_bounds();
          auto lb = pair_copulas_[t][e].get_parameters_lower_bounds();
          for (size_t p = 0; p < static_cast<size_t>(pars.size()); p++) {
            auto pars_tmp = pars;
            double eps = 0;

            pars_tmp(p) = std::min(pars(p) + 1e-3, ub(p));
            eps += pars_tmp(p) - pars(p);
            pair_copulas_[t][e].set_parameters(pars_tmp);
            Eigen::VectorXd f1 =
              this->pdf(u, num_threads).array().max(1e-20).log();

            pars_tmp(p) = std::max(pars(p) - 1e-3, lb(p));
            eps -= pars_tmp(p) - pars(p);
            pair_copulas_[t][e].set_parameters(pars_tmp);
            Eigen::VectorXd f2 =
              this->pdf(u, num_threads).array().max(1e-20).log();

            result.scores.col(ipar++) = (f1 - f2) / eps;
            pair_copulas_[t][e].set_parameters(pars);
          }
        }
      }

      return result;
    }

    // Analytic full gradient: the log-density of the vine is a sum over
    // edges whose arguments are h-functions of earlier trees, so a
    // parameter of edge (t0, e0) also enters all deeper edges fed by it.
    // A parameter's cascade seeds the perturbation of its own edge's
    // h-functions (dh/dtheta) and propagates it to deeper edges by the chain
    // rule, accumulating
    //   score      += (dc/du1)/c * darg1 + (dc/du2)/c * darg2,
    //   v_h2' = c * darg1 + (dh2/du2) * darg2,
    //   v_h1' = (dh1/du1) * darg1 + c * darg2,
    // where dh1(u2|u1)/du2 = dh2(u1|u2)/du1 = c was used. All per-edge
    // derivatives come from build_deriv_cache().
    //
    // Reference: Stoeber, J. and U. Schepsmeier (2013), Estimating standard
    // errors in regular vine copula models, Computational Statistics, 28
    // (6), 2679-2707 (the C implementation is VineCopula's rvinederiv.c).
    if (keep_all) {
      // per-edge caches returned to the caller (filled per batch below);
      // argument log-derivatives exist only for trees > 0 (where a
      // perturbation can arrive) and h-function derivatives only where a
      // deeper non-truncated tree consumes the output
      result.pdf_edges = TriangularArray<Eigen::VectorXd>(d_, trunc_lvl);
      result.logpdf_deriv_pars =
        TriangularArray<std::vector<Eigen::VectorXd>>(d_, trunc_lvl);
      result.hfunc1_deriv_pars =
        TriangularArray<std::vector<Eigen::VectorXd>>(d_, trunc_lvl);
      result.hfunc2_deriv_pars =
        TriangularArray<std::vector<Eigen::VectorXd>>(d_, trunc_lvl);
      result.logpdf_deriv_u1 = TriangularArray<Eigen::VectorXd>(d_, trunc_lvl);
      result.logpdf_deriv_u2 = TriangularArray<Eigen::VectorXd>(d_, trunc_lvl);
      result.hfunc1_deriv_u1 = TriangularArray<Eigen::VectorXd>(d_, trunc_lvl);
      result.hfunc2_deriv_u2 = TriangularArray<Eigen::VectorXd>(d_, trunc_lvl);
      for (size_t t = 0; t < trunc_lvl; ++t) {
        for (size_t e = 0; e < d_ - 1 - t; ++e) {
          size_t np =
            static_cast<size_t>(pair_copulas_[t][e].get_parameters().size());
          result.pdf_edges(t, e) = Eigen::VectorXd(n);
          result.logpdf_deriv_pars(t, e).resize(np);
          result.hfunc1_deriv_pars(t, e).resize(np);
          result.hfunc2_deriv_pars(t, e).resize(np);
          for (size_t p = 0; p < np; ++p) {
            result.logpdf_deriv_pars(t, e)[p] = Eigen::VectorXd(n);
          }
          if (t > 0) {
            result.logpdf_deriv_u1(t, e) = Eigen::VectorXd(n);
            result.logpdf_deriv_u2(t, e) = Eigen::VectorXd(n);
          }
          bool has_deeper_tree = (t + 1 < trunc_lvl);
          if (has_deeper_tree && rvine_structure_.needed_hfunc1(t, e)) {
            result.hfunc1_deriv_u1(t, e) = Eigen::VectorXd(n);
            for (size_t p = 0; p < np; ++p) {
              result.hfunc1_deriv_pars(t, e)[p] = Eigen::VectorXd(n);
            }
          }
          if (has_deeper_tree && rvine_structure_.needed_hfunc2(t, e)) {
            result.hfunc2_deriv_u2(t, e) = Eigen::VectorXd(n);
            for (size_t p = 0; p < np; ++p) {
              result.hfunc2_deriv_pars(t, e)[p] = Eigen::VectorXd(n);
            }
          }
        }
      }
    }

    auto do_batch = [&](const tools_batch::Batch& b) {
      auto cache = build_deriv_cache(u, b.begin, b.size, false, per_obs_params);

      if (keep_all) {
        for (size_t t = 0; t < trunc_lvl; ++t) {
          for (size_t e = 0; e < d_ - 1 - t; ++e) {
            const DerivCache& ce = cache(t, e);
            result.pdf_edges(t, e).segment(b.begin, b.size) = ce.c;
            if (t > 0) {
              result.logpdf_deriv_u1(t, e).segment(b.begin, b.size) = ce.du1;
              result.logpdf_deriv_u2(t, e).segment(b.begin, b.size) = ce.du2;
            }
            for (size_t p = 0; p < ce.np; ++p) {
              result.logpdf_deriv_pars(t, e)[p].segment(b.begin, b.size) =
                ce.dpar[p];
            }
            if (ce.h1.active) {
              result.hfunc1_deriv_u1(t, e).segment(b.begin, b.size) = ce.h1.du1;
              for (size_t p = 0; p < ce.np; ++p) {
                result.hfunc1_deriv_pars(t, e)[p].segment(b.begin, b.size) =
                  ce.h1.dpar[p];
              }
            }
            if (ce.h2.active) {
              result.hfunc2_deriv_u2(t, e).segment(b.begin, b.size) = ce.h2.du2;
              for (size_t p = 0; p < ce.np; ++p) {
                result.hfunc2_deriv_pars(t, e)[p].segment(b.begin, b.size) =
                  ce.h2.dpar[p];
              }
            }
          }
        }
      }

      // Per-parameter cascade over the cache: for each parameter θ (of edge
      // (t0, e0)), walk the deeper trees and accumulate the full score
      // ∂/∂θ Σ_edges log c via the chain rule.
      //
      // State per tree level, indexed by storage column e:
      //  - dh1.col(e) / dh2.col(e): the current first-order
      //    perturbations ∂hfunc1/∂θ and ∂hfunc2/∂θ of column e's h-function
      //    values (hfunc1 and hfunc2, as in the pdf
      //    pass). They are seeded with ∂h/∂θ at the parameter's own edge
      //    and updated tree by tree via ce.h1/ce.h2 (the edge's output
      //    derivative leaves).
      //  - h1_affected[e] / h2_affected[e]: flags marking whether column e's
      //    hfunc1/hfunc2 values are affected by θ at the current level
      //    (i.e. whether the corresponding tilde column is valid). Edges
      //    whose two input columns are both unaffected are skipped and
      //    contribute nothing.
      Eigen::MatrixXd dh1(b.size, d_), dh2(b.size, d_);
      std::vector<char> h1_affected(d_), h2_affected(d_);
      size_t ipar = 0;
      for (size_t t0 = 0; t0 < trunc_lvl; ++t0) {
        for (size_t e0 = 0; e0 < d_ - t0 - 1; ++e0) {
          const DerivCache& seed = cache(t0, e0);
          for (size_t p = 0; p < seed.np; ++p) {
            Eigen::VectorXd score = seed.dpar[p];
            std::fill(h1_affected.begin(), h1_affected.end(), 0);
            std::fill(h2_affected.begin(), h2_affected.end(), 0);
            if (seed.h1.active) {
              dh1.col(e0) = seed.h1.dpar[p];
              h1_affected[e0] = 1;
            }
            if (seed.h2.active) {
              dh2.col(e0) = seed.h2.dpar[p];
              h2_affected[e0] = 1;
            }

            for (size_t t = t0 + 1; t < trunc_lvl; ++t) {
              for (size_t e = 0; e < d_ - t - 1; ++e) {
                const DerivCache& ce = cache(t, e);
                size_t arg2_col = ce.arg2_col;
                // is either input argument of this edge perturbed by θ?
                // (first argument = column e's hfunc2; second argument =
                // column arg2_col - 1's hfunc2 or hfunc1, depending on
                // `arg2_is_h2`). The proximity condition guarantees arg2_col -
                // 1 > e, so columns read here are not yet overwritten at this
                // tree level.
                bool arg1_affected = (h2_affected[e] != 0);
                bool arg2_affected = ce.arg2_is_h2
                                       ? (h2_affected[arg2_col - 1] != 0)
                                       : (h1_affected[arg2_col - 1] != 0);
                if (!arg1_affected && !arg2_affected) {
                  // θ does not reach this edge; its outputs are unperturbed
                  h1_affected[e] = 0;
                  h2_affected[e] = 0;
                  continue;
                }
                // darg1/darg2 = ∂(arg1)/∂θ and ∂(arg2)/∂θ for this edge
                Eigen::VectorXd darg1 = Eigen::VectorXd::Zero(b.size);
                Eigen::VectorXd darg2 = Eigen::VectorXd::Zero(b.size);
                if (arg1_affected) {
                  darg1 = dh2.col(e);
                  // chain rule: this edge's log-density term responds to
                  // the perturbation of its first argument
                  score += ce.du1.cwiseProduct(darg1);
                }
                if (arg2_affected) {
                  darg2 = ce.arg2_is_h2 ? dh2.col(arg2_col - 1)
                                        : dh1.col(arg2_col - 1);
                  score += ce.du2.cwiseProduct(darg2);
                }
                // has_deeper_tree the argument perturbations through this
                // edge's h-function outputs for the next tree level: ∂h/∂θ =
                // ∂h/∂u1 · darg1 + ∂h/∂u2 · darg2
                if (ce.h2.active) {
                  dh2.col(e) =
                    propagate_first_order(ce.h2.du1, ce.h2.du2, darg1, darg2);
                }
                if (ce.h1.active) {
                  dh1.col(e) =
                    propagate_first_order(ce.h1.du1, ce.h1.du2, darg1, darg2);
                }
                h1_affected[e] = ce.h1.active ? 1 : 0;
                h2_affected[e] = ce.h2.active ? 1 : 0;
              }
            }
            result.scores.col(ipar++).segment(b.begin, b.size) = score;
          }
        }
      }
    };

    tools_thread::ThreadPool pool((num_threads == 1) ? 0 : num_threads);
    pool.map(do_batch, tools_batch::create_batches(n, num_threads));
    pool.join();

    return result;
  }

  // step-wise scores: only the per-edge log-density derivatives are needed
  if (keep_all) {
    result.logpdf_deriv_pars =
      TriangularArray<std::vector<Eigen::VectorXd>>(d_, trunc_lvl);
    for (size_t t = 0; t < trunc_lvl; ++t) {
      for (size_t e = 0; e < d_ - 1 - t; ++e) {
        size_t npars =
          static_cast<size_t>(pair_copulas_[t][e].get_parameters().size());
        result.logpdf_deriv_pars(t, e).resize(npars);
        for (size_t p = 0; p < npars; ++p) {
          result.logpdf_deriv_pars(t, e)[p] = Eigen::VectorXd(n);
        }
      }
    }
  }

  auto do_batch = [&](const tools_batch::Batch& b) {
    // temporary storage objects (all data must be in (0, 1))
    Eigen::MatrixXd hfunc1, hfunc2, u_e, hfunc1_sub, hfunc2_sub, u_e_sub;
    hfunc1 = Eigen::MatrixXd::Zero(b.size, d_);
    hfunc2 = Eigen::MatrixXd::Zero(b.size, d_);
    if (get_n_discrete() > 0) {
      hfunc1_sub = hfunc1;
      hfunc2_sub = hfunc2;
    }

    // fill first row of hfunc2 matrix with evaluation points;
    // points have to be reordered to correspond to natural order
    for (size_t j = 0; j < d_; ++j) {
      hfunc2.col(j) = u.block(b.begin, order[j] - 1, b.size, 1);
      if (var_types_[order[j] - 1] == "d") {
        hfunc2_sub.col(j) =
          u.block(b.begin, d_ + disc_cols[order[j] - 1], b.size, 1);
      }
    }

    size_t ipar = 0;
    for (size_t tree = 0; tree < trunc_lvl; ++tree) {
      tools_interface::check_user_interrupt(
        static_cast<double>(u.rows()) * static_cast<double>(d_) > 1e3);
      for (size_t edge = 0; edge < d_ - tree - 1; ++edge) {
        tools_interface::check_user_interrupt(edge % 100 == 0);
        // extract evaluation point from hfunction matrices (have been
        // computed in previous tree level)
        Bicop edge_copula = get_pair_copula(tree, edge);
        auto var_types = edge_copula.get_var_types();
        size_t m = rvine_structure_.min_array(tree, edge);

        u_e = Eigen::MatrixXd(b.size, 2);
        u_e.col(0) = hfunc2.col(edge);
        if (m == rvine_structure_.struct_array(tree, edge, true)) {
          u_e.col(1) = hfunc2.col(m - 1);
        } else {
          u_e.col(1) = hfunc1.col(m - 1);
        }

        if ((var_types[0] == "d") || (var_types[1] == "d")) {
          u_e.conservativeResize(b.size, 4);
          u_e.col(2) = hfunc2_sub.col(edge);
          if (m == rvine_structure_.struct_array(tree, edge, true)) {
            u_e.col(3) = hfunc2_sub.col(m - 1);
          } else {
            u_e.col(3) = hfunc1_sub.col(m - 1);
          }
        }

        auto pars = edge_copula.get_parameters();
        // this edge's per-observation parameter slice (its columns start at the
        // current flat index `ipar`, which equals edge_par0(tree, edge) here)
        Eigen::MatrixXd pars_e;
        if (per_obs) {
          pars_e = per_obs_params.block(b.begin, ipar, b.size, pars.size());
        }
        if (var_types == std::vector<std::string>{ "c", "c" }) {
          // analytic per-edge gradient of the log-density (closed forms for
          // bicop_families::analytic_derivs, internal finite differences of
          // the density otherwise); nonparametric edges were rejected above
          for (size_t p = 0; p < static_cast<size_t>(pars.size()); p++) {
            Eigen::VectorXd col =
              per_obs
                ? edge_copula.logpdf_deriv(
                    u_e, "par" + std::to_string(p + 1), pars_e)
                : edge_copula.logpdf_deriv(u_e, "par" + std::to_string(p + 1));
            result.scores.col(ipar).segment(b.begin, b.size) = col;
            if (keep_all) {
              result.logpdf_deriv_pars(tree, edge)[p].segment(b.begin, b.size) =
                col;
            }
            ipar++;
          }
        } else {
          // discrete edges: finite differences of the edge log-density
          auto ub = edge_copula.get_parameters_upper_bounds();
          auto lb = edge_copula.get_parameters_lower_bounds();
          for (size_t p = 0; p < static_cast<size_t>(pars.size()); p++) {
            auto pars_tmp = pars;
            double eps = 0;

            pars_tmp(p) = std::min(pars(p) + 1e-3, ub(p));
            eps += pars_tmp(p) - pars(p);
            edge_copula.set_parameters(pars_tmp);
            Eigen::VectorXd f1 = edge_copula.pdf(u_e).array().max(1e-20).log();

            pars_tmp(p) = std::max(pars(p) - 1e-3, lb(p));
            eps -= pars_tmp(p) - pars(p);
            edge_copula.set_parameters(pars_tmp);
            Eigen::VectorXd f2 = edge_copula.pdf(u_e).array().max(1e-20).log();

            Eigen::VectorXd col = (f1 - f2) / eps;
            result.scores.col(ipar).segment(b.begin, b.size) = col;
            if (keep_all) {
              result.logpdf_deriv_pars(tree, edge)[p].segment(b.begin, b.size) =
                col;
            }
            ipar++;
            edge_copula.set_parameters(pars);
          }
        }

        // h-functions are only evaluated if needed in next step (per_obs uses
        // the edge slice so deeper trees stay consistent; discrete sub-columns
        // never arise on the per_obs path, which requires continuous variables)
        if (rvine_structure_.needed_hfunc1(tree, edge)) {
          hfunc1.col(edge) =
            per_obs ? edge_copula.hfunc1(u_e, pars_e) : edge_copula.hfunc1(u_e);
          if (var_types[1] == "d") {
            u_e_sub = u_e;
            u_e_sub.col(1) = u_e.col(3);
            hfunc1_sub.col(edge) = edge_copula.hfunc1(u_e_sub);
          }
        }
        if (rvine_structure_.needed_hfunc2(tree, edge)) {
          hfunc2.col(edge) =
            per_obs ? edge_copula.hfunc2(u_e, pars_e) : edge_copula.hfunc2(u_e);
          if (var_types[0] == "d") {
            u_e_sub = u_e;
            u_e_sub.col(0) = u_e.col(2);
            hfunc2_sub.col(edge) = edge_copula.hfunc2(u_e_sub);
          }
        }
      }
    }
  };

  tools_thread::ThreadPool pool((num_threads == 1) ? 0 : num_threads);
  pool.map(do_batch, tools_batch::create_batches(n, num_threads));
  pool.join();

  return result;
}

//! @brief Evaluates the score function.
//!
//! The score function is defined as the gradient of the log-likelihood
//! with respect to the parameters. This is a thin wrapper around
//! `scores_full()`; see there for the computational details.
//!
//! @param u An \f$ n \times 2d \f$ or \f$ n \times (d + k) \f$ matrix of
//!   evaluation points, where \f$ k \f$ is the number of discrete variables
//!   (see `select()`).
//! @param step_wise if `false`, full gradient of the log-likelihood; if `true`,
//!   score function of the step-wise MLE (gradients computed per pair-copula,
//!   treating the pseudo-observations as fixed).
//! @param num_threads The number of threads to use for computations; if greater
//!   than 1, the function will be applied concurrently to `num_threads` batches
//!   of `u`.
inline Eigen::MatrixXd
Vinecop::scores(Eigen::MatrixXd u, bool step_wise, const size_t num_threads)
{
  return scores_full(std::move(u), step_wise, num_threads, false).scores;
}

//! @brief Evaluates the score function with per-observation parameters.
//!
//! Per-observation counterpart of `scores()`; see the per-observation
//! `scores_full()` overload for the `parameters` layout and restrictions.
inline Eigen::MatrixXd
Vinecop::scores(Eigen::MatrixXd u,
                const Eigen::MatrixXd& parameters,
                bool step_wise,
                const size_t num_threads)
{
  return scores_full(std::move(u), parameters, step_wise, num_threads, false)
    .scores;
}

//! @brief Evaluates the gradient of the average log-likelihood.
//!
//! Returns the observation-average of `scores()` as a vector of length
//! `npars`, mirroring how `hessian()` averages `hessian_full()`.
//!
//! @param u An \f$ n \times 2d \f$ or \f$ n \times (d + k) \f$ matrix of
//!   evaluation points, where \f$ k \f$ is the number of discrete variables
//!   (see `select()`).
//! @param step_wise if `false`, full gradient of the log-likelihood; if `true`,
//!   score function of the step-wise MLE (gradients computed per pair-copula,
//!   treating the pseudo-observations as fixed).
//! @param num_threads The number of threads to use for computations; if greater
//!   than 1, the function will be applied concurrently to `num_threads` batches
//!   of `u`.
inline Eigen::VectorXd
Vinecop::gradient(Eigen::MatrixXd u, bool step_wise, const size_t num_threads)
{
  return this->scores(std::move(u), step_wise, num_threads)
    .colwise()
    .mean()
    .transpose();
}

//! @brief Evaluates the gradient of the average log-likelihood with
//! per-observation parameters.
//!
//! Per-observation counterpart of `gradient()`; the return is the
//! observation-average of the per-observation `scores()`.
inline Eigen::VectorXd
Vinecop::gradient(Eigen::MatrixXd u,
                  const Eigen::MatrixXd& parameters,
                  bool step_wise,
                  const size_t num_threads)
{
  return this->scores(std::move(u), parameters, step_wise, num_threads)
    .colwise()
    .mean()
    .transpose();
}

//! @brief Evaluates the hessian per observation.
//!
//! Hessian is meant loosely as "gradients of each component of the score
//! function", i.e. `hess(t, e)[p](i, a) = ∂² log-likelihood_i / (∂θ_{t,e,p}
//! ∂θ_a)`. For a continuous model it is computed analytically from the pair
//! copulas' first and second derivatives: the joint (non-step-wise) Hessian
//! by a second-order cascade through the vine (the second derivative of the
//! RVineGrad-style gradient cascade in `scores()`), and the step-wise Hessian
//! by a first-order cascade of the step-wise score's argument derivatives.
//! Models with discrete variables use central finite differences of
//! `scores()` instead; models with nonparametric pair copulas are rejected.
//!
//! @param u An \f$ n \times 2d \f$ or \f$ n \times (d + k) \f$ matrix of
//!   evaluation points, where \f$ k \f$ is the number of discrete variables
//!   (see `select()`).
//! @param step_wise if `false`, full gradient of the log-likelihood; if `true`,
//!   score function of the step-wise MLE (gradients computed per pair-copula).
//! @param num_threads The number of threads to use for computations; if greater
//!   than 1, the function will be applied concurrently to `num_threads` batches
//!   of `u`.
//!
//! @literature
//! Stoeber, J. and Schepsmeier, U. (2013). Estimating standard errors in
//! regular vine copula models. Computational Statistics, 28 (6), 2679-2707.
inline TriangularArray<std::vector<Eigen::MatrixXd>>
Vinecop::hessian_full(Eigen::MatrixXd u,
                      bool step_wise,
                      const size_t num_threads)
{
  return hessian_full(std::move(u), Eigen::MatrixXd(), step_wise, num_threads);
}

//! @brief Evaluates the per-observation Hessians with per-observation
//! parameters.
//!
//! Per-observation counterpart of `hessian_full()`; see the per-observation
//! `scores_full()` overload for the `parameters` layout and restrictions.
inline TriangularArray<std::vector<Eigen::MatrixXd>>
Vinecop::hessian_full(Eigen::MatrixXd u,
                      const Eigen::MatrixXd& per_obs_params,
                      bool step_wise,
                      const size_t num_threads)
{
  check_data(u);
  u = collapse_data(u);

  size_t trunc_lvl = rvine_structure_.get_trunc_lvl();

  check_parametric("hessian()");

  const bool per_obs = per_obs_params.size() > 0;
  if (per_obs) {
    check_per_obs_params(u, per_obs_params);
  }

  TriangularArray<std::vector<Eigen::MatrixXd>> hess(d_);

  if (get_n_discrete() > 0) {
    // central finite differences of the (analytic-where-possible) scores;
    // used for models with discrete variables (both step-wise and joint)
    for (size_t t = 0; t < trunc_lvl; t++) {
      for (size_t e = 0; e < d_ - 1 - t; e++) {
        auto pars = pair_copulas_[t][e].get_parameters();
        auto ub = pair_copulas_[t][e].get_parameters_upper_bounds();
        auto lb = pair_copulas_[t][e].get_parameters_lower_bounds();
        hess(t, e).resize(pars.size());
        for (size_t p = 0; p < static_cast<size_t>(pars.size()); p++) {
          auto pars_tmp = pars;
          double eps = 0;

          pars_tmp(p) = std::min(pars(p) + 1e-3, ub(p));
          eps += pars_tmp(p) - pars(p);
          pair_copulas_[t][e].set_parameters(pars_tmp);
          Eigen::MatrixXd f1 = this->scores(u, step_wise, num_threads);

          pars_tmp(p) = std::max(pars(p) - 1e-3, lb(p));
          eps -= pars_tmp(p) - pars(p);
          pair_copulas_[t][e].set_parameters(pars_tmp);
          Eigen::MatrixXd f2 = this->scores(u, step_wise, num_threads);

          hess(t, e)[p] = (f1 - f2) / eps;
          pair_copulas_[t][e].set_parameters(pars);
        }
      }
    }
    return hess;
  }

  // Analytic Hessian (continuous, all-parametric). Both variants use the
  // per-edge derivative caches from build_deriv_cache() and the shared
  // first-order propagation; see hessian_full()'s @literature reference.
  const size_t npars = static_cast<size_t>(this->get_npars());

  // flat parameter index <-> (tree, edge, param), the first flat index of
  // each edge, and the per-edge per-observation output blocks
  std::vector<std::array<size_t, 3>> par_of(npars);
  TriangularArray<size_t> edge_par0(d_, trunc_lvl);
  {
    size_t idx = 0;
    for (size_t t = 0; t < trunc_lvl; t++) {
      for (size_t e = 0; e < d_ - 1 - t; e++) {
        edge_par0(t, e) = idx;
        size_t np =
          static_cast<size_t>(pair_copulas_[t][e].get_parameters().size());
        hess(t, e).resize(np);
        for (size_t p = 0; p < np; p++) {
          hess(t, e)[p] = Eigen::MatrixXd::Zero(u.rows(), npars);
          par_of[idx++] = { t, e, p };
        }
      }
    }
  }

  if (step_wise) {
    // Analytic step-wise Hessian: the step-wise score of edge e', parameter
    // p', is s = ∂logc_{e'}/∂θ_{p'} evaluated at fixed pseudo-observations,
    // so its derivative w.r.t. a parameter α of edge (tα, eα) is only
    //   ∂²logc_{eα}/∂θ_p∂θ_{p'}          (α's own edge, u_e fixed), or
    //   ∂²logc_{e'}/∂θ_{p'}∂u1 · P^α + ∂²logc/∂θ_{p'}∂u2 · Q^α   (deeper e'),
    // where P^α, Q^α are the first-order perturbations of e''s arguments,
    // propagated by the same cascade as the gradient. No second-order
    // ("bar") propagation is needed. Upper-triangular in the flat order.
    auto do_batch = [&](const tools_batch::Batch& b) {
      const size_t m = b.size;
      auto cache = build_deriv_cache(u, b.begin, m, true, per_obs_params);

      // dh1/dh2 and h1_affected/h2_affected play the same role as in the
      // gradient cascade of scores_full: per storage column, the current
      // first-order perturbations of the column's hfunc1/hfunc2 values w.r.t.
      // the parameter, and flags marking which columns are affected (see the
      // DerivLeaf/DerivCache documentation in class.hpp)
      Eigen::MatrixXd dh1(m, d_), dh2(m, d_);
      std::vector<char> h1_affected(d_), h2_affected(d_);
      for (size_t a = 0; a < npars; ++a) {
        size_t ta = par_of[a][0], ea = par_of[a][1], pa = par_of[a][2];
        const DerivCache& seed = cache(ta, ea);

        // own edge: ∂ s_{(ta,ea,cp)} / ∂θ_{(ta,ea,pa)} = ∂²logc/∂θ_pa∂θ_cp
        for (size_t cp = 0; cp < seed.np; ++cp) {
          hess(ta, ea)[pa].block(b.begin, edge_par0(ta, ea) + cp, m, 1) =
            seed.dpar_par[pa][cp];
        }

        // propagate α's first-order perturbation to the deeper edges
        std::fill(h1_affected.begin(), h1_affected.end(), 0);
        std::fill(h2_affected.begin(), h2_affected.end(), 0);
        if (seed.h1.active) {
          dh1.col(ea) = seed.h1.dpar[pa];
          h1_affected[ea] = 1;
        }
        if (seed.h2.active) {
          dh2.col(ea) = seed.h2.dpar[pa];
          h2_affected[ea] = 1;
        }
        for (size_t t = ta + 1; t < trunc_lvl; ++t) {
          for (size_t e = 0; e < d_ - t - 1; ++e) {
            const DerivCache& ce = cache(t, e);
            size_t arg2_col = ce.arg2_col;
            bool arg1_affected = (h2_affected[e] != 0);
            bool arg2_affected = ce.arg2_is_h2
                                   ? (h2_affected[arg2_col - 1] != 0)
                                   : (h1_affected[arg2_col - 1] != 0);
            if (!arg1_affected && !arg2_affected) {
              h1_affected[e] = 0;
              h2_affected[e] = 0;
              continue;
            }
            Eigen::VectorXd darg1 = Eigen::VectorXd::Zero(m);
            Eigen::VectorXd darg2 = Eigen::VectorXd::Zero(m);
            if (arg1_affected) {
              darg1 = dh2.col(e);
            }
            if (arg2_affected) {
              darg2 =
                ce.arg2_is_h2 ? dh2.col(arg2_col - 1) : dh1.col(arg2_col - 1);
            }
            for (size_t cp = 0; cp < ce.np; ++cp) {
              hess(ta, ea)[pa].block(b.begin, edge_par0(t, e) + cp, m, 1) =
                ce.dpar_u1[cp].cwiseProduct(darg1) +
                ce.dpar_u2[cp].cwiseProduct(darg2);
            }
            if (ce.h2.active) {
              dh2.col(e) =
                propagate_first_order(ce.h2.du1, ce.h2.du2, darg1, darg2);
            }
            if (ce.h1.active) {
              dh1.col(e) =
                propagate_first_order(ce.h1.du1, ce.h1.du2, darg1, darg2);
            }
            h1_affected[e] = ce.h1.active ? 1 : 0;
            h2_affected[e] = ce.h2.active ? 1 : 0;
          }
        }
      }
    };

    if (trunc_lvl > 0) {
      tools_thread::ThreadPool pool((num_threads == 1) ? 0 : num_threads);
      pool.map(do_batch, tools_batch::create_batches(u.rows(), num_threads));
      pool.join();
    }
    return hess;
  }

  // Analytic joint Hessian: differentiate the gradient cascade a second
  // time. For a parameter pair (α, β) a forward walk propagates the first
  // derivatives of each edge's arguments w.r.t. α ("hat") and β ("tilde")
  // and their mixed second derivative ("bar"), accumulating
  //   H_αβ = Σ_e [ ∂²logc/∂p²·P^α P^β + ∂²logc/∂p∂q·(P^α Q^β + Q^α P^β)
  //               + ∂²logc/∂q²·Q^α Q^β + ∂logc/∂p·P^{αβ} + ∂logc/∂q·Q^{αβ}
  //               + θ-argument cross terms and the ∂²logc/∂θ² seed ],
  // where the bar quantities propagate with the h-function's own second
  // derivatives.
  auto do_batch = [&](const tools_batch::Batch& b) {
    const size_t m = b.size;
    // one forward walk fills every per-edge first- and second-order leaf
    auto cache = build_deriv_cache(u, b.begin, m, true, per_obs_params);

    // Second-order cascade for each parameter pair (only the upper triangle;
    // H is symmetric). Per storage column e (suffix 1 = the column's hfunc1
    // values, 2 = its hfunc2 values, as in the gradient cascade):
    //   hat*[e] = ∂(column value)/∂α,
    //   til*[e] = ∂(column value)/∂β,
    //   bar*[e] = ∂²(column value)/∂α∂β,
    // seeded at the two parameters' own edges and pushed through each
    // deeper edge's output leaves ce.h1/ce.h2 (see class.hpp). There are no
    // aff flags here; unaffected columns simply carry zero vectors.
    std::vector<Eigen::VectorXd> dh1_a(d_), dh2_a(d_), dh1_b(d_), dh2_b(d_),
      d2h1(d_), d2h2(d_);
    for (size_t a = 0; a < npars; ++a) {
      size_t ta = par_of[a][0], ea = par_of[a][1], pa = par_of[a][2];
      for (size_t bfl = a; bfl < npars; ++bfl) {
        size_t tb = par_of[bfl][0], eb = par_of[bfl][1], pb = par_of[bfl][2];

        for (size_t col = 0; col < d_; ++col) {
          dh1_a[col] = Eigen::VectorXd::Zero(m);
          dh2_a[col] = Eigen::VectorXd::Zero(m);
          dh1_b[col] = Eigen::VectorXd::Zero(m);
          dh2_b[col] = Eigen::VectorXd::Zero(m);
          d2h1[col] = Eigen::VectorXd::Zero(m);
          d2h2[col] = Eigen::VectorXd::Zero(m);
        }
        Eigen::VectorXd H_ab = Eigen::VectorXd::Zero(m);

        for (size_t t = 0; t < trunc_lvl; ++t) {
          for (size_t e = 0; e < d_ - 1 - t; ++e) {
            const DerivCache& ce = cache(t, e);
            size_t arg2_col = ce.arg2_col;
            bool arg2_is_h2 = ce.arg2_is_h2;
            bool is_edge_a = (t == ta) && (e == ea);
            bool is_edge_b = (t == tb) && (e == eb);

            // input-argument perturbations (p_e from column e's hfunc2,
            // q_e from column arg2_col-1's hfunc2/hfunc1). p_e's column is
            // overwritten by the propagation below, so copy it by value;
            // q_e's column (arg2_col-1 > e) is untouched this step.
            Eigen::VectorXd darg1_a = dh2_a[e];
            Eigen::VectorXd darg1_b = dh2_b[e];
            Eigen::VectorXd d2arg1 = d2h2[e];
            const Eigen::VectorXd& darg2_a =
              arg2_is_h2 ? dh2_a[arg2_col - 1] : dh1_a[arg2_col - 1];
            const Eigen::VectorXd& darg2_b =
              arg2_is_h2 ? dh2_b[arg2_col - 1] : dh1_b[arg2_col - 1];
            const Eigen::VectorXd& d2arg2 =
              arg2_is_h2 ? d2h2[arg2_col - 1] : d2h1[arg2_col - 1];

            // accumulate the Hessian contribution of this edge
            H_ab += ce.du1u1.cwiseProduct(darg1_a).cwiseProduct(darg1_b);
            H_ab += ce.du1u2.cwiseProduct(darg1_a.cwiseProduct(darg2_b) +
                                          darg2_a.cwiseProduct(darg1_b));
            H_ab += ce.du2u2.cwiseProduct(darg2_a).cwiseProduct(darg2_b);
            H_ab += ce.du1.cwiseProduct(d2arg1) + ce.du2.cwiseProduct(d2arg2);
            if (is_edge_a) {
              H_ab += ce.dpar_u1[pa].cwiseProduct(darg1_b) +
                      ce.dpar_u2[pa].cwiseProduct(darg2_b);
            }
            if (is_edge_b) {
              H_ab += ce.dpar_u1[pb].cwiseProduct(darg1_a) +
                      ce.dpar_u2[pb].cwiseProduct(darg2_a);
            }
            if (is_edge_a && is_edge_b) {
              H_ab += ce.dpar_par[pa][pb];
            }

            // propagate the perturbations to this edge's h-function outputs:
            // hat/til are first order (propagate_first_order), bar is second
            auto propagate = [&](const DerivLeaf& leaf,
                                 Eigen::VectorXd& out_a,
                                 Eigen::VectorXd& out_b,
                                 Eigen::VectorXd& out_ab) {
              out_a =
                propagate_first_order(leaf.du1, leaf.du2, darg1_a, darg2_a);
              out_b =
                propagate_first_order(leaf.du1, leaf.du2, darg1_b, darg2_b);
              out_ab = leaf.du1u1.cwiseProduct(darg1_a).cwiseProduct(darg1_b) +
                       leaf.du1u2.cwiseProduct(darg1_a.cwiseProduct(darg2_b) +
                                               darg2_a.cwiseProduct(darg1_b)) +
                       leaf.du2u2.cwiseProduct(darg2_a).cwiseProduct(darg2_b) +
                       leaf.du1.cwiseProduct(d2arg1) +
                       leaf.du2.cwiseProduct(d2arg2);
              if (is_edge_a) {
                out_a += leaf.dpar[pa];
                out_ab += leaf.dpar_u1[pa].cwiseProduct(darg1_b) +
                          leaf.dpar_u2[pa].cwiseProduct(darg2_b);
              }
              if (is_edge_b) {
                out_b += leaf.dpar[pb];
                out_ab += leaf.dpar_u1[pb].cwiseProduct(darg1_a) +
                          leaf.dpar_u2[pb].cwiseProduct(darg2_a);
              }
              if (is_edge_a && is_edge_b) {
                out_ab += leaf.dpar_par[pa][pb];
              }
            };
            if (ce.h2.active) {
              propagate(ce.h2, dh2_a[e], dh2_b[e], d2h2[e]);
            }
            if (ce.h1.active) {
              propagate(ce.h1, dh1_a[e], dh1_b[e], d2h1[e]);
            }
          }
        }

        hess(ta, ea)[pa].block(b.begin, bfl, m, 1) = H_ab;
        if (a != bfl) {
          hess(tb, eb)[pb].block(b.begin, a, m, 1) = H_ab;
        }
      }
    }
  };

  if (trunc_lvl > 0) {
    tools_thread::ThreadPool pool((num_threads == 1) ? 0 : num_threads);
    pool.map(do_batch, tools_batch::create_batches(u.rows(), num_threads));
    pool.join();
  }

  return hess;
}

//! @brief Evaluates the averaged Hessian of the log-likelihood.
//!
//! Hessian is meant loosely as "gradients of each component of the score
//! function". The Hessian is averaged over all samples in `u`, yielding an
//! \f$ \mathrm{npars} \times \mathrm{npars} \f$ matrix (use `hessian_full()`
//! for the per-observation decomposition).
//!
//! @param u An \f$ n \times 2d \f$ or \f$ n \times (d + k) \f$ matrix of
//!   evaluation points, where \f$ k \f$ is the number of discrete variables
//!   (see `select()`).
//! @param step_wise if `false`, full gradient of the log-likelihood; if `true`,
//!   score function of the step-wise MLE (gradients computed per pair-copula).
//! @param num_threads The number of threads to use for computations; if greater
//!   than 1, the function will be applied concurrently to `num_threads` batches
//!   of `u`.
inline Eigen::MatrixXd
Vinecop::hessian(Eigen::MatrixXd u, bool step_wise, const size_t num_threads)
{
  return hessian(std::move(u), Eigen::MatrixXd(), step_wise, num_threads);
}

//! @brief Evaluates the averaged Hessian with per-observation parameters.
//!
//! Per-observation counterpart of `hessian()`; see the per-observation
//! `scores_full()` overload for the `parameters` layout and restrictions.
inline Eigen::MatrixXd
Vinecop::hessian(Eigen::MatrixXd u,
                 const Eigen::MatrixXd& parameters,
                 bool step_wise,
                 const size_t num_threads)
{
  const bool per_obs = parameters.size() > 0;
  // validate up front so the per-chunk parameters.middleRows() slices below
  // are always in range (hessian_full re-validates each chunk, but only after
  // the slice is taken)
  if (per_obs) {
    check_per_obs_params(u, parameters);
  }
  const size_t npars = static_cast<size_t>(this->get_npars());
  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(npars, npars);
  const size_t n = static_cast<size_t>(u.rows());
  if (n == 0 || npars == 0) {
    return H;
  }
  size_t trunc_lvl = rvine_structure_.get_trunc_lvl();

  // ponytail: process observations in row-chunks so peak memory is
  // O(chunk * npars^2) rather than O(n * npars^2) -- hessian_full()
  // materialises the per-observation Hessian, which hessian() only needs in
  // aggregate. The chunk caps that at ~4M doubles; call hessian_full()
  // directly (or raise the cap) if you need the full per-observation tensor.
  // The overall analytic cascade is O(n * d^6) (npars^2 pairs x O(d^2)
  // edges) for one-parameter families.
  const size_t chunk =
    std::max<size_t>(1, 4000000 / std::max<size_t>(1, npars * npars));
  for (size_t begin = 0; begin < n; begin += chunk) {
    size_t size = std::min(chunk, n - begin);
    auto hess =
      per_obs
        ? this->hessian_full(u.middleRows(begin, size),
                             parameters.middleRows(begin, size),
                             step_wise,
                             num_threads)
        : this->hessian_full(u.middleRows(begin, size), step_wise, num_threads);
    size_t ipar = 0;
    for (size_t t = 0; t < trunc_lvl; t++) {
      for (size_t e = 0; e < d_ - 1 - t; e++) {
        for (const auto& block : hess(t, e)) {
          H.row(ipar++) += block.colwise().sum();
        }
      }
    }
  }
  return H / static_cast<double>(n);
}

//! @brief Computes the covariance matrix of scores.
//!
//! Returns the (mean-centered, divided by `n`) covariance of the
//! per-observation scores as an \f$ \mathrm{npars} \times \mathrm{npars} \f$
//! matrix. Together with `hessian()` this forms the sandwich estimator of the
//! asymptotic covariance.
//!
//! @param u An \f$ n \times 2d \f$ or \f$ n \times (d + k) \f$ matrix of
//!   evaluation points, where \f$ k \f$ is the number of discrete variables
//!   (see `select()`).
//! @param step_wise if `false`, full gradient of the log-likelihood; if `true`,
//!   score function of the step-wise MLE (gradients computed per pair-copula).
//! @param num_threads The number of threads to use for computations; if greater
//!   than 1, the function will be applied concurrently to `num_threads` batches
//!   of `u`.
inline Eigen::MatrixXd
Vinecop::scores_cov(Eigen::MatrixXd u, bool step_wise, const size_t num_threads)
{
  auto s = this->scores(std::move(u), step_wise, num_threads);
  // materialize the centered scores; a lazy expression would be evaluated
  // twice by the product below
  Eigen::MatrixXd sc = s.rowwise() - s.colwise().mean();
  return (sc.adjoint() * sc) / static_cast<double>(s.rows());
}

//! @brief Computes the covariance matrix of scores with per-observation
//! parameters.
//!
//! Per-observation counterpart of `scores_cov()`; see the per-observation
//! `scores_full()` overload for the `parameters` layout and restrictions.
inline Eigen::MatrixXd
Vinecop::scores_cov(Eigen::MatrixXd u,
                    const Eigen::MatrixXd& parameters,
                    bool step_wise,
                    const size_t num_threads)
{
  auto s = this->scores(std::move(u), parameters, step_wise, num_threads);
  Eigen::MatrixXd sc = s.rowwise() - s.colwise().mean();
  return (sc.adjoint() * sc) / static_cast<double>(s.rows());
}

//! @brief Evaluates the copula distribution.
//!
//! @details Because no closed-form expression is available, the distribution is
//! estimated numerically using Monte Carlo integration. The function uses
//! quasi-random numbers from the vine model to do so.
//!
//! When at least one variable is discrete, two types of
//! "observations" are required in `u`: the first \f$ n \; x \; d \f$ block
//! contains realizations of \f$ F_{X_j}(X_j) \f$.
//! The second \f$ n \; x \; d \f$
//! block contains realizations of \f$ F_{X_j}(X_j^-) \f$. The minus indicates a
//! left-sided limit of the cdf. For, e.g., an integer-valued variable, it holds
//! \f$ F_{X_j}(X_j^-) = F_{X_j}(X_j - 1) \f$. For continuous variables the left
//! limit and the cdf itself coincide. Respective columns can be omitted in the
//! second block.
//!
//! @param u An \f$ n \times 2d \f$ or \f$ n \times (d + k) \f$ matrix of
//!   evaluation points, where \f$ k \f$ is the number of discrete variables
//!   (see `Vinecop::select()`).
//! @param N Integer for the number of quasi-random numbers to draw
//! to evaluate the distribution (default: 1e4).
//! @param num_threads The number of threads to use for computations; if greater
//!   than 1, the function will generate `n` samples concurrently in
//!   `num_threads` batches.
//! @param seeds Seeds to scramble the quasi-random numbers; if empty (default),
//!   the random number quasi-generator is seeded randomly.
//! @return A vector of length `n` containing the copula distribution values.
inline Eigen::VectorXd
Vinecop::cdf(const Eigen::MatrixXd& u,
             const size_t N,
             const size_t num_threads,
             const std::vector<int>& seeds) const
{
  if (d_ > 21201) {
    std::stringstream message;
    message << "cumulative distribution available for models of "
            << "dimension 21201 or less. This model's dimension: " << d_
            << std::endl;
    throw std::runtime_error(message.str().c_str());
  }
  check_data(u);

  // Simulate N quasi-random numbers from the vine model
  auto u_sim = simulate(N, true, num_threads, seeds);

  size_t n = u.rows();
  Eigen::VectorXd vine_distribution(n);
  auto do_batch = [&](const tools_batch::Batch& b) {
    // lazy dominance count per query point (no N-sized temporary)
    Eigen::RowVectorXd temp(d_);
    for (size_t i = b.begin; i < b.begin + b.size; i++) {
      tools_interface::check_user_interrupt(i % 1000 == 0);
      temp = u.block(i, 0, 1, d_);
      vine_distribution(i) = static_cast<double>(
        ((u_sim.rowwise() - temp).rowwise().maxCoeff().array() <= 0.0).count());
    }
  };
  tools_thread::ThreadPool pool((num_threads == 1) ? 0 : num_threads);
  pool.map(do_batch, tools_batch::create_batches(n, num_threads));
  pool.join();
  return vine_distribution / static_cast<double>(N);
}

//! @brief Simulates from a vine copula model, see `inverse_rosenblatt()`.
//!
//! @details Simulated data is always a continous \f$ n \times d \f$ matrix.
//! Sampling from a vine copula model is done by first generating
//! \f$ n \times d \f$ uniform random numbers and then applying the inverse
//! Rosenblatt transformation.
//!
//! @param n Number of observations.
//! @param qrng Set to true for quasi-random numbers.
//! @param num_threads The number of threads to use for computations; if greater
//!   than 1, the function will generate `n` samples concurrently in
//!   `num_threads` batches.
//! @param seeds Seeds of the random number generator; if empty (default),
//!   the random number generator is seeded randomly.
//! @return An \f$ n \times d \f$ matrix of samples from the copula model.
inline Eigen::MatrixXd
Vinecop::simulate(const size_t n,
                  const bool qrng,
                  const size_t num_threads,
                  const std::vector<int>& seeds) const
{
  auto u = tools_stats::simulate_uniform(n, d_, qrng, seeds);
  return inverse_rosenblatt(u, num_threads);
  ;
}

//! @brief Simulates from the conditional distribution of a subset of variables
//! given fixed values of the remaining variables.
//!
//! @details The conditioning variables are the last `k` variables of the vine
//! order (they are drawn first by the vine and form a self-contained sub-vine).
//! Each row of `u_cond` is one
//! conditioning point; the corresponding output row is drawn from the
//! remaining variables' distribution conditional on that point.
//!
//! Discrete conditioning variables are supported. As for the other data-taking
//! methods, a discrete variable requires its left-limit CDF \f$ F(x^-) \f$ in
//! addition to \f$ F(x) \f$; `u_cond` therefore carries one extra column per
//! discrete conditioning variable (see the `u_cond` layout below). The
//! conditioned variables are then drawn using the randomized transform of
//! `rosenblatt()` (see its `randomize_discrete` note).
//!
//! @param u_cond An \f$ n \times (k + k_d) \f$ matrix of conditioning values,
//!   where `k` is the number of conditioning variables and `k_d` the number of
//!   discrete ones among them. Row `j` is the conditioning point for output
//!   sample `j` (`n` is the number of rows of `u_cond`). The first `k` columns
//!   hold the values \f$ F(x) \f$; column `i` corresponds to the
//!   `(d - k + i)`-th variable of the vine order, i.e. the columns correspond,
//!   left to right, to the last `k` entries of the order. The next `k_d`
//!   columns hold the left-limits \f$ F(x^-) \f$ of the discrete conditioning
//!   variables, in the order those variables appear among the first `k`
//!   columns. For all-continuous conditioning this is simply an
//!   \f$ n \times k \f$ matrix. `k` is inferred from the column count;
//!   supplying only `k` columns when a conditioning variable is discrete may be
//!   silently reinterpreted as a different `k`. To draw many samples at a
//!   single conditioning point, pass that point repeated over `n` rows.
//!   The overload taking `conditioning_set` also accepts the expanded
//!   \f$ n \times 2k \f$ layout, with one left-limit column per conditioning
//!   variable. Use that overload whenever the expanded layout is ambiguous
//!   with a compact layout for a different number of conditioning variables.
//! @param qrng Set to true for quasi-random numbers (over the conditioned
//!   variables).
//! @param num_threads The number of threads to use for computations.
//! @param seeds Seeds of the random number generator; if empty (default),
//!   the random number generator is seeded randomly.
//! @return An \f$ n \times d \f$ matrix; the remaining columns are draws from
//!   the conditional distribution. The conditioning columns reproduce `u_cond`
//!   exactly for continuous variables. A discrete conditioning variable is
//!   reproduced up to its atom (the column lands in \f$ [F(x^-), F(x)] \f$);
//!   with several discrete conditioning variables the later-drawn ones may land
//!   slightly outside their atom, but the draws of the remaining (conditioned)
//!   variables remain correct.
inline Eigen::MatrixXd
Vinecop::simulate_conditional(const Eigen::MatrixXd& u_cond,
                              const bool qrng,
                              const size_t num_threads,
                              const std::vector<int>& seeds) const
{
  size_t n_cols = static_cast<size_t>(u_cond.cols());
  auto order = rvine_structure_.get_order();

  size_t k = 0;
  for (size_t kk = 1; kk <= d_ - 1; ++kk) {
    size_t kd = 0;
    for (size_t i = 0; i < kk; ++i)
      if (var_types_[order[d_ - kk + i] - 1] == "d")
        ++kd;
    if (n_cols == kk + kd) {
      k = kk;
      break;
    }
  }
  if (k == 0) {
    throw std::runtime_error(
      "u_cond has an invalid number of columns; expected k + (number of "
      "discrete conditioning variables) columns for some number of "
      "conditioning variables k in 1, ..., d - 1.");
  }

  std::vector<size_t> conditioning_set(k);
  for (size_t i = 0; i < k; ++i)
    conditioning_set[i] = order[d_ - k + i];
  return simulate_conditional_impl(
    u_cond, conditioning_set, VinecopView(*this), qrng, num_threads, seeds);
}

//! @brief Simulates conditionally on a specified set of variables.
//!
//! @param conditioning_set The 1-based conditioning-variable indices. The
//!   first `k` columns of `u_cond` correspond to these variables in the given
//!   order. In the expanded layout, the next `k` columns contain their
//!   left-limits in the same order; in the compact layout, only left-limits of
//!   discrete variables follow.
inline Eigen::MatrixXd
Vinecop::simulate_conditional(const Eigen::MatrixXd& u_cond,
                              const std::vector<size_t>& conditioning_set,
                              const bool qrng,
                              const size_t num_threads,
                              const std::vector<int>& seeds) const
{
  auto reorientation = make_reorientation_map(conditioning_set);
  return simulate_conditional_impl(u_cond,
                                   conditioning_set,
                                   VinecopView(*this, reorientation),
                                   qrng,
                                   num_threads,
                                   seeds);
}

inline Eigen::MatrixXd
Vinecop::simulate_conditional_impl(const Eigen::MatrixXd& u_cond,
                                   const std::vector<size_t>& conditioning_set,
                                   const VinecopView& view,
                                   bool qrng,
                                   size_t num_threads,
                                   const std::vector<int>& seeds) const
{
  const size_t n = static_cast<size_t>(u_cond.rows());
  const size_t k = conditioning_set.size();
  const size_t kd = static_cast<size_t>(std::count_if(
    conditioning_set.begin(), conditioning_set.end(), [&](size_t var) {
      return var_types_[var - 1] == "d";
    }));
  const size_t n_cols = static_cast<size_t>(u_cond.cols());
  const bool expanded = n_cols == 2 * k;
  if (!expanded && (n_cols != k + kd)) {
    throw std::runtime_error(
      "u_cond must have 2 * k columns or one column per conditioning variable "
      "plus one left-limit column per discrete conditioning variable.");
  }
  if (!tools_eigen::check_if_in_unit_cube(u_cond)) {
    throw std::runtime_error("all elements of u_cond must be in (0, 1).");
  }

  auto disc_cols = tools_select::get_disc_cols(var_types_);
  Eigen::MatrixXd u_completed =
    Eigen::MatrixXd::Constant(n, d_ + get_n_discrete(), 0.5);
  size_t kd_seen = 0;
  for (size_t i = 0; i < k; ++i) {
    size_t var = conditioning_set[i] - 1;
    u_completed.col(var) = u_cond.col(i);
    if (var_types_[var] == "d") {
      size_t left_col = expanded ? k + i : k + kd_seen;
      if ((u_cond.col(left_col).array() > u_cond.col(i).array()).any()) {
        throw std::runtime_error(
          "for discrete conditioning variables, the left-limit columns of "
          "u_cond (F(x^-)) must not exceed the value columns (F(x)).");
      }
      u_completed.col(d_ + disc_cols[var]) = u_cond.col(left_col);
      ++kd_seen;
    }
  }

  Eigen::MatrixXd w =
    rosenblatt_impl(std::move(u_completed), view, num_threads, true, seeds);

  Eigen::MatrixXd u = tools_stats::simulate_uniform(n, d_, qrng, seeds);
  for (auto var : conditioning_set) {
    u.col(var - 1) = w.col(var - 1);
  }

  return inverse_rosenblatt_impl(u, view, num_threads);
}

//! @brief Evaluates the log-likelihood.
//!
//! @details The log-likelihood is defined as
//! \f[ \mathrm{loglik} = \sum_{i = 1}^n \log c(U_{1, i}, ..., U_{d, i}), \f]
//! where \f$ c \f$ is the copula density, see `Vinecop::pdf()`.
//!
//! @param u An \f$ n \times 2d \f$ or \f$ n \times (d + k) \f$ matrix of
//!   evaluation points, where \f$ k \f$ is the number of discrete variables
//!   (see `select()` or `Vinecop::pdf()`).
//! @param num_threads The number of threads to use for computations; if greater
//!   than 1, the function will be applied concurrently to `num_threads` batches
//!   of `u`.
//! @return The log-likelihood as a double.
inline double
Vinecop::loglik(const Eigen::MatrixXd& u, const size_t num_threads) const
{
  if (u.rows() < 1) {
    return this->get_loglik();
  } else {
    return pdf(u, num_threads).array().log().sum();
  }
}

//! @brief Evaluates the log-likelihood with per-observation parameters.
//!
//! Per-observation counterpart of `loglik()`: the sum over observations of
//! \f$ \log c(u_i; \theta_i) \f$, where row `i` of `parameters` is
//! \f$ \theta_i \f$. See the per-observation `pdf_full()` overload for the
//! `parameters` layout and restrictions.
inline double
Vinecop::loglik(const Eigen::MatrixXd& u,
                const Eigen::MatrixXd& parameters,
                const size_t num_threads) const
{
  return pdf(u, parameters, num_threads).array().log().sum();
}

//! @brief Evaluates the Akaike information criterion (AIC).
//!
//! @details The AIC is defined as
//! \f[ \mathrm{AIC} = -2\, \mathrm{loglik} + 2 p, \f]
//! where \f$ \mathrm{loglik} \f$ is the log-liklihood (see `Vinecop::loglik()`)
//! and \f$ p \f$ is the (effective) number of parameters of the model.
//! The AIC is a consistent model selection criterion even
//! for nonparametric models.
//!
//! @param u An \f$ n \times 2d \f$ or \f$ n \times (d + k) \f$ matrix of
//!   evaluation points, where \f$ k \f$ is the number of discrete variables
//!   (see `select()` or `Vinecop::pdf()`).
//! @param num_threads The number of threads to use for computations; if greater
//!   than 1, the function will be applied concurrently to `num_threads` batches
//!   of `u`.
//! @return The AIC as a double.
inline double
Vinecop::aic(const Eigen::MatrixXd& u, const size_t num_threads) const
{
  return -2 * this->loglik(u, num_threads) + 2 * get_npars();
}

//! @brief Evaluates the Bayesian information criterion (BIC).
//!
//! @details The BIC is defined as
//! \f[ \mathrm{BIC} = -2\, \mathrm{loglik} +  \log(n) p, \f]
//! where \f$ \mathrm{loglik} \f$ is the log-liklihood (see `Vinecop::loglik()`)
//! and \f$ p \f$ is the (effective) number of parameters of the model.
//! The BIC is a consistent model selection criterion
//! for nonparametric models.
//!
//! @param u An \f$ n \times 2d \f$ or \f$ n \times (d + k) \f$ matrix of
//!   evaluation points, where \f$ k \f$ is the number of discrete variables
//!   (see `Vinecop::select()` or `Vinecop::pdf()`).
//! @param num_threads The number of threads to use for computations; if greater
//!   than 1, the function will be applied concurrently to `num_threads` batches
//!   of `u`.
//! @return The BIC as a double.
inline double
Vinecop::bic(const Eigen::MatrixXd& u, const size_t num_threads) const
{
  return -2 * this->loglik(u, num_threads) +
         get_npars() * log(static_cast<double>(u.rows()));
}

// clang-format off
//! @brief Evaluates the modified Bayesian information criterion for vines
//! (mBICV).
//!
//! @details The mBICV is defined as
//! \f[ \mathrm{mBICV} = -2\, \mathrm{loglik} +  \log(n) p, - 2 * \sum_{t=1}^(d - 1) \{q_t \log(\psi_0^t) - (d - t - q_t) \log(1 -\psi_0^t)\},\f]
//! where \f$ \mathrm{loglik} \f$ is the log-liklihood,
//! \f$ p \f$ is the (effective) number of parameters of the model, \f$ t \f$
//! is the tree level, \f$ \psi_0 \f$ is the prior probability of having a
//! non-independence copula in the first tree, and \f$ q_t \f$ is the number of
//! non-independence copulas in tree \f$ t \f$; The vBIC is a consistent model
//! selection criterion for parametric sparse vine copula models when
//! \f$ d = o(\sqrt{n \log n})\f$.
//!
//! @param u An \f$ n \times 2d \f$ or \f$ n \times (d + k) \f$ matrix of
//!   evaluation points, where \f$ k \f$ is the number of discrete variables
//!   (see `Vinecop::select()` or `Vinecop::pdf()`).
//! @param psi0 Baseline prior probability of a non-independence copula.
//! @param num_threads The number of threads to use for computations; if greater
//!   than 1, the function will be applied concurrently to `num_threads` batches
//!   of `u`.
//! @return The mBICV as a double.
// clang-format on
inline double
Vinecop::mbicv(const Eigen::MatrixXd& u,
               const double psi0,
               const size_t num_threads) const
{

  size_t n = u.rows();
  double ll = this->loglik(u, num_threads);
  return -2.0 * ll + this->calculate_mbicv_penalty(n, psi0);
}

//! @brief Returns sum of the number of parameters for all pair copulas (see.
//! Bicop::get_npars()).
inline double
Vinecop::get_npars() const
{
  double npars = 0.0;
  for (auto& tree : pair_copulas_) {
    for (auto& pc : tree) {
      npars += pc.get_npars();
    }
  }
  return npars;
}

//! @brief Evaluates the Rosenblatt transform for a vine copula model.
//!
//! @details The Rosenblatt transform converts data from this model
//! into independent uniform variates.
//!
//! The Rosenblatt transform (Rosenblatt, 1952) \f$ U = T(V) \f$ of a random
//! vector \f$ V = (V_1,\ldots,V_d) ~ F \f$ is defined as
//! \f[ U_1= F(V_1), U_{2} = F(V_{2}|V_1), \ldots, U_d
//! =F(V_d|V_1,\ldots,V_{d-1}), \f] where \f$ F(v_k|v_1,\ldots,v_{k-1}) \f$ is
//! the conditional distribution of \f$ V_k \f$ given  \f$ V_1 \ldots, V_{k-1},
//! k = 2,\ldots,d \f$. The vector \f$ U = (U_1, \dots, U_d) \f$ then contains
//! independent standard uniform variables. The inverse operation \f[ V_1 =
//! F^{-1}(U_1), V_{2} = F^{-1}(U_2|U_1), \ldots, V_d
//! =F^{-1}(U_d|U_1,\ldots,U_{d-1}) \f] can be used to simulate from a
//! distribution. For any copula \f$ F \f$, if \f$ U\f$ is a vector of
//! independent random variables, \f$ V = T^{-1}(U) \f$ has distribution \f$ F
//! \f$.
//!
//! The formulas above assume a vine copula model with order \f$ d, \dots, 1
//! \f$. More generally, `Vinecop::rosenblatt()` returns the variables \f[
//! U_{M[d - j, j]}= F(V_{M[d - j, j]} | V_{M[d - j - 1, j - 1]}, \dots, V_{M[0,
//! 0]}), \f] where \f$ M \f$ is the structure matrix.  Similarly,
//! `Vinecop::inverse_rosenblatt()` computes \f[ V_{M[d - j, j]}= F^{-1}(U_{M[d
//! - j, j]} | U_{M[d - j - 1, j - 1]}, \dots, U_{M[0, 0]}). \f]
//!
//! If some variables have atoms, Brockwell (10.1016/j.spl.2007.02.008) proposed
//! a simple randomization scheme to ensure that output is still independent
//! uniform if the model is correct. The transformation reads \f[ U_{M[d - j,
//! j]}= W_{d - j} F(V_{M[d - j, j]} | V_{M[d - j - 1, j - 1]}, \dots, V_{M[0,
//! 0]}) + (1 - W_{d - j}) F^-(V_{M[d - j, j]} | V_{M[d - j - 1, j - 1]}, \dots,
//! V_{M[0, 0]}), \f] where \f$ F^- \f$ is the left limit of the conditional cdf
//! and \f$W_1, \dots, W_d \f$ are are independent standard uniform random
//! variables. This is used by default. If you are interested in the conditional
//! probabilities
//! \f[ F(V_{M[d - j, j]} | V_{M[d - j - 1, j - 1]}, \dots, V_{M[0, 0]}), \f]
//! set `randomize_discrete = FALSE`.
//!
//! @param u An \f$ n \times 2d \f$ or \f$ n \times (d + k) \f$ matrix of
//!   evaluation points, where \f$ k \f$ is the number of discrete variables.
//! @param num_threads The number of threads to use for computations; if greater
//!   than 1, the function will be applied concurrently to `num_threads` batches
//!   of `u`.
//! @param randomize_discrete Whether to randomize the transform for discrete
//!   variables; see Details.
//! @param seeds Seeds to scramble the quasi-random numbers; if empty (default),
//!   the random number quasi-generator is seeded randomly. Only relevant if
//!   there are discrete variables and `randomize_discrete = TRUE`.
//!
//! @return An \f$ n \times d \f$ matrix of independent uniform variates.
inline Eigen::MatrixXd
Vinecop::rosenblatt(Eigen::MatrixXd u,
                    const size_t num_threads,
                    bool randomize_discrete,
                    std::vector<int> seeds) const
{
  return rosenblatt_impl(std::move(u),
                         VinecopView(*this),
                         num_threads,
                         randomize_discrete,
                         std::move(seeds));
}

//! @brief Evaluates the Rosenblatt transform for a conditioning set.
//!
//! @details The vine is evaluated in an admissible sampling order whose tail
//! contains exactly `conditioning_set`. The model itself is not modified.
//!
//! @param u An \f$ n \times 2d \f$ or \f$ n \times (d + k) \f$ matrix of
//!   evaluation points, where \f$ k \f$ is the number of discrete variables.
//! @param conditioning_set The 1-based indices of the conditioning variables.
//! @param num_threads The number of threads to use for computations.
//! @param randomize_discrete Whether to randomize the transform for discrete
//!   variables.
//! @param seeds Seeds used for discrete randomization.
//! @return An \f$ n \times d \f$ matrix of independent uniform variates.
inline Eigen::MatrixXd
Vinecop::rosenblatt(Eigen::MatrixXd u,
                    const std::vector<size_t>& conditioning_set,
                    const size_t num_threads,
                    bool randomize_discrete,
                    std::vector<int> seeds) const
{
  auto reorientation = make_reorientation_map(conditioning_set);
  return rosenblatt_impl(std::move(u),
                         VinecopView(*this, reorientation),
                         num_threads,
                         randomize_discrete,
                         std::move(seeds));
}

inline Eigen::MatrixXd
Vinecop::rosenblatt_impl(Eigen::MatrixXd u,
                         const VinecopView& view,
                         size_t num_threads,
                         bool randomize_discrete,
                         std::vector<int> seeds) const
{
  check_data(u);
  collapse_data_inplace(u);

  size_t d = d_;
  size_t n = u.rows();

  // info about the vine structure
  const auto& structure = view.get_structure();
  size_t trunc_lvl = structure.get_trunc_lvl();
  auto order = structure.get_order();
  auto inverse_order = tools_stl::invert_permutation(order);
  auto disc_cols = tools_select::get_disc_cols(var_types_);

  // fill first row of hfunc2 matrix with evaluation points;
  // points have to be reordered to correspond to natural order
  Eigen::MatrixXd hfunc1(n, d), hfunc2(n, d), hfunc1_sub, hfunc2_sub;
  for (size_t j = 0; j < d; ++j) {
    hfunc2.col(j) = u.col(order[j] - 1);
  }
  hfunc1 = hfunc2; // just ensure data is in [0, 1]^d
  if (is_discrete()) {
    hfunc1_sub = hfunc1;
    hfunc2_sub = hfunc2;
    for (size_t j = 0; j < d_; ++j) {
      if (var_types_[order[j] - 1] == "d") {
        hfunc2_sub.col(j) = u.col(d_ + disc_cols[order[j] - 1]);
      }
    }
  }

  auto do_batch = [&](const tools_batch::Batch& b) {
    Eigen::MatrixXd u_e, u_e_sub;
    for (size_t tree = 0; tree < trunc_lvl; ++tree) {
      tools_interface::check_user_interrupt(
        static_cast<double>(n) * static_cast<double>(d) > 1e5);
      for (size_t edge = 0; edge < d - tree - 1; ++edge) {
        tools_interface::check_user_interrupt(edge % 100 == 0);
        // extract evaluation point from hfunction matrices (have been
        // computed in previous tree level)
        auto edge_copula = view.get_pair_copula(tree, edge);
        auto var_types = edge_copula.get_var_types();
        size_t m = structure.min_array(tree, edge);

        u_e.resize(b.size, 2);
        u_e.col(0) = hfunc2.block(b.begin, edge, b.size, 1);
        if (m == structure.struct_array(tree, edge, true)) {
          u_e.col(1) = hfunc2.block(b.begin, m - 1, b.size, 1);
        } else {
          u_e.col(1) = hfunc1.block(b.begin, m - 1, b.size, 1);
        }

        if ((var_types[0] == "d") || (var_types[1] == "d")) {
          u_e.conservativeResize(b.size, 4);
          u_e.col(2) = hfunc2_sub.block(b.begin, edge, b.size, 1);
          if (m == structure.struct_array(tree, edge, true)) {
            u_e.col(3) = hfunc2_sub.block(b.begin, m - 1, b.size, 1);
          } else {
            u_e.col(3) = hfunc1_sub.block(b.begin, m - 1, b.size, 1);
          }
        }

        // h-functions are only evaluated if needed in next step
        if (structure.needed_hfunc1(tree, edge)) {
          hfunc1.block(b.begin, edge, b.size, 1) = edge_copula.hfunc1(u_e);
          if (var_types[1] == "d") {
            u_e_sub = u_e;
            u_e_sub.col(1) = u_e.col(3);
            hfunc1_sub.block(b.begin, edge, b.size, 1) =
              edge_copula.hfunc1(u_e_sub);
          }
        }
        hfunc2.block(b.begin, edge, b.size, 1) = edge_copula.hfunc2(u_e);
        if (var_types[0] == "d") {
          u_e_sub = u_e;
          u_e_sub.col(0) = u_e.col(2);
          hfunc2_sub.block(b.begin, edge, b.size, 1) =
            edge_copula.hfunc2(u_e_sub);
        }
      }
    }
  };

  if (trunc_lvl > 0) {
    tools_thread::ThreadPool pool((num_threads == 1) ? 0 : num_threads);
    pool.map(do_batch, tools_batch::create_batches(n, num_threads));
    pool.join();
  }

  // go back to original order
  Eigen::MatrixXd U(n, 2 * d);
  for (size_t j = 0; j < d; j++) {
    U.col(j) = hfunc2.col(inverse_order[j]);
  }

  if (randomize_discrete && is_discrete()) {
    // fill second half of U with left-sided limits of the conditional CDF
    // (equal to conditional CDF for continuous variables)
    for (size_t j = 0; j < d; j++) {
      U.col(d + j) = var_types_[j] == "d" ? hfunc2_sub.col(inverse_order[j])
                                          : hfunc2.col(inverse_order[j]);
    }
    // randomize by weighting left and right limits with independent uniforms
    auto R =
      tools_stats::simulate_uniform(u.rows(), d, false, std::move(seeds));
    U.leftCols(d) = U.leftCols(d).array() * R.array() +
                    U.rightCols(d).array() * (1 - R.array());
  }
  return U.leftCols(d).array().min(1 - 1e-10).max(1e-10);
}

//! @brief Evaluates the inverse Rosenblatt transform.
//!
//! @details The inverse Rosenblatt transform can be used for simulation: the
//! function applied to independent uniform variates resembles simulated
//! data from the vine copula model.
//!
//! If the problem is too large, it is split recursively into halves (w.r.t.
//! \f$ n \f$, the number of observations).
//! "Too large" means that the required memory will exceed 1 GB. An
//! examplary configuration requiring less than 1 GB is \f$ n = 1000 \f$,
//! \f$ d = 200\f$.
//!
//! The Rosenblatt transform (Rosenblatt, 1952) \f$ U = T(V) \f$ of a random
//! vector \f$ V = (V_1,\ldots,V_d) ~ F \f$ is defined as
//! \f[ U_1= F(V_1), U_{2} = F(V_{2}|V_1), \ldots, U_d
//! =F(V_d|V_1,\ldots,V_{d-1}), \f] where \f$ F(v_k|v_1,\ldots,v_{k-1}) \f$ is
//! the conditional distribution of \f$ V_k \f$ given  \f$ V_1 \ldots, V_{k-1},
//! k = 2,\ldots,d \f$. The vector \f$ U = (U_1, \dots, U_d) \f$ then contains
//! independent standard uniform variables. The inverse operation \f[ V_1 =
//! F^{-1}(U_1), V_{2} = F^{-1}(U_2|U_1), \ldots, V_d
//! =F^{-1}(U_d|U_1,\ldots,U_{d-1}) \f] can be used to simulate from a
//! distribution. For any copula \f$ F \f$, if \f$ U\f$ is a vector of
//! independent random variables, \f$ V = T^{-1}(U) \f$ has distribution \f$ F
//! \f$.
//!
//! The formulas above assume a vine copula model with order \f$ d, \dots, 1
//! \f$. More generally, `Vinecop::rosenblatt()` returns the variables \f[
//! U_{M[d - j, j]}= F(V_{M[d - j, j]} | V_{M[d - j - 1, j - 1]}, \dots, V_{M[0,
//! 0]}), \f] where \f$ M \f$ is the structure matrix. Similarly,
//! `Vinecop::inverse_rosenblatt()` computes \f[ V_{M[d - j, j]}= F^{-1}(U_{M[d
//! - j, j]} | U_{M[d - j - 1, j - 1]}, \dots, U_{M[0, 0]}). \f]
//!
//! @param u An \f$ n \times 2d \f$, \f$ n \times (d + k) \f$, or
//!   \f$ n \times d \f$ matrix of independent uniform variates, where \f$ k \f$
//!   is the number of discrete variables. Only the first \f$ d \f$ columns are
//!   used.
//! @param num_threads The number of threads to use for computations; if greater
//!   than 1, the function will be applied concurrently to `num_threads` batches
//!   of `u`.
//! @return An \f$ n \times d \f$ matrix of evaluations.
inline Eigen::MatrixXd
Vinecop::inverse_rosenblatt(const Eigen::MatrixXd& u,
                            const size_t num_threads) const
{
  return inverse_rosenblatt_impl(u, VinecopView(*this), num_threads);
}

//! @brief Evaluates the inverse Rosenblatt transform for a conditioning set.
//!
//! @details The vine is evaluated in an admissible sampling order whose tail
//! contains exactly `conditioning_set`. The model itself is not modified.
//!
//! @param u An \f$ n \times 2d \f$, \f$ n \times (d + k) \f$, or
//!   \f$ n \times d \f$ matrix of independent uniform variates, where \f$ k \f$
//!   is the number of discrete variables. Only the first \f$ d \f$ columns are
//!   used.
//! @param conditioning_set The 1-based indices of the conditioning variables.
//! @param num_threads The number of threads to use for computations.
//! @return An \f$ n \times d \f$ matrix of transformed values.
inline Eigen::MatrixXd
Vinecop::inverse_rosenblatt(const Eigen::MatrixXd& u,
                            const std::vector<size_t>& conditioning_set,
                            const size_t num_threads) const
{
  auto reorientation = make_reorientation_map(conditioning_set);
  return inverse_rosenblatt_impl(
    u, VinecopView(*this, reorientation), num_threads);
}

inline Eigen::MatrixXd
Vinecop::inverse_rosenblatt_impl(const Eigen::MatrixXd& u,
                                 const VinecopView& view,
                                 size_t num_threads) const
{
  const size_t n_cols = static_cast<size_t>(u.cols());
  const size_t compact_cols = d_ + get_n_discrete();
  if ((n_cols != d_) && (n_cols != compact_cols) && (n_cols != 2 * d_)) {
    throw std::runtime_error(
      "data has wrong number of columns; expected: " + std::to_string(2 * d_) +
      ", " + std::to_string(compact_cols) + ", or " + std::to_string(d_) +
      ", actual: " + std::to_string(n_cols) + ".");
  }
  if (u.rows() < 1) {
    throw std::runtime_error("data must have at least one row");
  }
  tools_eigen::check_if_in_unit_cube(u);

  size_t n = u.rows();
  size_t d = d_;

  Eigen::MatrixXd U_vine = u.leftCols(d); // output matrix
  //                   (hfunc1 + hfunc2)      (U_vine)       (info matrices)
  size_t bytes_required =
    (size_t{ 16 } * n * d * d) + (size_t{ 8 } * n * d) + (size_t{ 16 } * d * d);
  // if the problem is too large (requires more than 1 GB memory), split
  // the data into two halves and call simulate on the reduced data.
  if ((n > 1) & (bytes_required > static_cast<size_t>(1e9))) {
    size_t n_half = n / 2;
    size_t n_left = n - n_half;
    U_vine.block(0, 0, n_half, d) =
      inverse_rosenblatt_impl(u.block(0, 0, n_half, d), view, num_threads);
    U_vine.block(n_half, 0, n_left, d) =
      inverse_rosenblatt_impl(u.block(n_half, 0, n_left, d), view, num_threads);
    return U_vine;
  }

  // info about the vine structure (in upper triangular matrix notation)
  const auto& structure = view.get_structure();
  size_t trunc_lvl = structure.get_trunc_lvl();
  auto order = structure.get_order();
  auto inverse_order = tools_stl::invert_permutation(order);

  auto do_batch = [&](const tools_batch::Batch& b) {
    // temporary storage objects for (inverse) h-functions
    TriangularArray<Eigen::VectorXd> hinv2(d + 1, trunc_lvl + 1);
    TriangularArray<Eigen::VectorXd> hfunc1(d + 1, trunc_lvl + 1);
    Eigen::MatrixXd U_e(b.size, 2);

    // initialize with independent uniforms (corresponding to natural
    // order)
    for (size_t j = 0; j < d; ++j) {
      hinv2(std::min(trunc_lvl, d - j - 1), j) =
        u.block(b.begin, order[j] - 1, b.size, 1);
    }
    hfunc1(0, d - 1) = hinv2(0, d - 1);

    // loop through variables (0 is just the initial uniform)
    for (ptrdiff_t var = d - 2; var >= 0; --var) {
      tools_interface::check_user_interrupt(
        static_cast<double>(n) * static_cast<double>(d) > 1e5);
      size_t tree_start = std::min(trunc_lvl - 1, d - var - 2);
      for (ptrdiff_t tree = tree_start; tree >= 0; --tree) {
        auto edge_copula = view.get_pair_copula(tree, var).as_continuous();

        // extract data for conditional pair
        size_t m = structure.min_array(tree, var);
        U_e.col(0) = hinv2(tree + 1, var);
        if (m == structure.struct_array(tree, var, true)) {
          U_e.col(1) = hinv2(tree, m - 1);
        } else {
          U_e.col(1) = hfunc1(tree, m - 1);
        }

        // inverse Rosenblatt transform simulates data for conditional pair
        hinv2(tree, var) = edge_copula.hinv2(U_e);

        // if required at a later stage, also calculate hfunc1
        if (var < static_cast<ptrdiff_t>(d_) - 1) {
          if (structure.needed_hfunc1(tree, var)) {
            U_e.col(0) = hinv2(tree, var);
            hfunc1(tree + 1, var) = edge_copula.hfunc1(U_e);
          }
        }
      }
    }
    // go back to original order
    for (size_t j = 0; j < d; j++) {
      U_vine.block(b.begin, j, b.size, 1) = hinv2(0, inverse_order[j]);
    }
  };

  if (trunc_lvl > 0) {
    tools_thread::ThreadPool pool((num_threads == 1) ? 0 : num_threads);
    pool.map(do_batch, tools_batch::create_batches(n, num_threads));
    pool.join();
  }

  return U_vine;
}

//! Checks if dimension d of the data matches the dimension of the vine.
inline void
Vinecop::check_data_dim(const Eigen::MatrixXd& data) const
{
  size_t d_data = data.cols();
  auto n_disc = get_n_discrete();
  size_t d_exp = d_ + n_disc;
  if ((d_data != d_exp) & (d_data != 2 * d_)) {
    std::stringstream msg;
    msg << "data has wrong number of columns; "
        << "expected: " << 2 * d_ << " or " << d_exp << ", actual: " << d_data
        << " (model contains ";
    if (n_disc == 0) {
      msg << "no discrete variables)." << std::endl;
    } else if (n_disc == 1) {
      msg << "1 discrete variable)." << std::endl;
    } else {
      msg << get_n_discrete() << " discrete variables)." << std::endl;
    }
    throw std::runtime_error(msg.str());
  }

  if (data.rows() < 1) {
    throw std::runtime_error("data must have at least one row");
  }
}

//! Checks if dimension d of the data matches the dimension of the vine.
inline void
Vinecop::check_data(const Eigen::MatrixXd& data) const
{
  check_data_dim(data);
  tools_eigen::check_if_in_unit_cube(data);
}

//! Checks if pair copulas are compatible with the R-vine structure.
inline void
Vinecop::check_pair_copulas_rvine_structure(
  const std::vector<std::vector<Bicop>>& pair_copulas) const
{
  size_t trunc_lvl = rvine_structure_.get_trunc_lvl();
  if (pair_copulas.size() > std::min(d_ - 1, trunc_lvl)) {
    std::stringstream message;
    message << "pair_copulas is too large; "
            << "expected size: < " << std::min(d_ - 1, trunc_lvl) << ", "
            << "actual size: " << pair_copulas.size() << std::endl;
    throw std::runtime_error(message.str().c_str());
  }
  for (size_t t = 0; t < pair_copulas.size(); ++t) {
    if (pair_copulas[t].size() != d_ - 1 - t) {
      std::stringstream message;
      message << "size of pair_copulas[" << t << "] "
              << "does not match dimension of matrix (" << d_ << "); "
              << "expected size: " << d_ - 1 - t << ", "
              << "actual size: " << pair_copulas[t].size() << std::endl;
      throw std::runtime_error(message.str().c_str());
    }
  }
}

inline void
Vinecop::finalize_fit(const tools_select::VinecopSelector& selector)
{
  rvine_structure_ = selector.get_rvine_structure();
  threshold_ = selector.get_threshold();
  loglik_ = selector.get_loglik();
  nobs_ = selector.get_nobs();
  pair_copulas_ = selector.get_pair_copulas();
}

//! Checks if weights are compatible with the data.
inline void
Vinecop::check_weights_size(const Eigen::VectorXd& weights,
                            const Eigen::MatrixXd& data) const
{
  if ((weights.size() > 0) && (weights.size() != data.rows())) {
    throw std::runtime_error("sizes of weights and data don't match.");
  }
}

//! Checks if data size is large enough.
inline void
Vinecop::check_enough_data(const Eigen::MatrixXd& data) const
{
  if (data.rows() == 1) {
    throw std::runtime_error("data must have more than one row");
  }
}

inline void
Vinecop::check_fitted() const
{
  if (std::isnan(loglik_)) {
    throw std::runtime_error("copula has not been fitted from data ");
  }
}

inline void
Vinecop::check_indices(const size_t tree, const size_t edge) const
{
  if (tree > d_ - 2) {
    std::stringstream message;
    message << "tree index out of bounds" << std::endl
            << "allowed: 0, ..., " << d_ - 2 << std::endl
            << "actual: " << tree << std::endl;
    throw std::runtime_error(message.str().c_str());
  }
  if (edge > d_ - tree - 2) {
    std::stringstream message;
    message << "edge index out of bounds" << std::endl
            << "allowed: 0, ..., " << d_ - tree - 2 << std::endl
            << "actual: " << edge << std::endl
            << "tree level: " << tree << std::endl;
    throw std::runtime_error(message.str().c_str());
  }
}

//! @brief Truncates the vine copula model.
//!
//! @details While model for a `d` dimensional random vector contains at most
//! `d-1` nested trees, this function extracts a sub-model based
//! on a given truncation level.
//!
//! If the model is already truncated at a level less than `trunc_lvl`,
//! the function does nothing.
//! @param trunc_lvl The truncation level.
inline void
Vinecop::truncate(size_t trunc_lvl)
{
  if (trunc_lvl < this->get_trunc_lvl()) {
    rvine_structure_.truncate(trunc_lvl);
    pair_copulas_.resize(trunc_lvl);
  }
}

//! @brief Sets all variable types to continuous.
//!
inline void
Vinecop::set_continuous_var_types()
{
  var_types_ = std::vector<std::string>(d_);
  for (auto& t : var_types_)
    t = "c";
  set_var_types_internal(var_types_);
}

//! @brief Returns the number of discrete variables.
inline int
Vinecop::get_n_discrete() const
{
  return n_discrete_;
}

inline bool
Vinecop::is_discrete() const
{
  return n_discrete_ > 0;
}

//! @brief Removes superfluous columns for continuous data.
inline Eigen::MatrixXd
Vinecop::collapse_data(const Eigen::MatrixXd& u) const
{
  if (static_cast<size_t>(u.cols()) == d_ + get_n_discrete()) {
    return u;
  }
  Eigen::MatrixXd u_new(u.rows(), d_ + get_n_discrete());
  u_new.leftCols(d_) = u.leftCols(d_);
  size_t disc_count = 0;
  for (size_t i = 0; i < d_; ++i) {
    if (var_types_[i] == "d") {
      u_new.col(d_ + disc_count++) = u.col(d_ + i);
    }
  }
  return u_new;
}

//! @brief Removes superfluous columns for continuous data, avoiding the
//! copy when the data is already in collapsed format.
inline void
Vinecop::collapse_data_inplace(Eigen::MatrixXd& u) const
{
  if (static_cast<size_t>(u.cols()) != d_ + get_n_discrete()) {
    u = collapse_data(u);
  }
}

//! @brief Summarizes the model into a string (can be used for printing).
//! @param trees A vector of tree indices to summarize; if empty, all trees.
inline std::string
Vinecop::str(const std::vector<size_t>& trees) const
{
  std::vector<size_t> trees_to_summarize;
  std::vector<size_t> all_trees(rvine_structure_.get_trunc_lvl());
  std::iota(all_trees.begin(), all_trees.end(), 0);
  if (trees.empty()) {
    trees_to_summarize = all_trees;
  } else {
    trees_to_summarize = tools_stl::intersect(all_trees, trees);
  }

  std::stringstream vinecop_str;
  vinecop_str << std::setprecision(2);
  vinecop_str << "Vinecop model with " << d_ << " variables" << std::endl;
  auto arr = rvine_structure_.get_struct_array();
  auto order = rvine_structure_.get_order();

  std::vector<std::string> trees_s;
  std::vector<std::string> edges;
  std::vector<std::string> conditioned_variables;
  std::vector<std::string> conditioning_variables;
  std::vector<std::string> var_types;
  std::vector<std::string> families;
  std::vector<std::string> rotations;
  std::vector<std::string> parameters;
  std::vector<std::string> dfs;
  std::vector<std::string> taus;
  trees_s.emplace_back("tree");
  edges.emplace_back("edge");
  conditioned_variables.emplace_back("conditioned variables");
  conditioning_variables.emplace_back("conditioning variables");
  var_types.emplace_back("var_types");
  families.emplace_back("family");
  rotations.emplace_back("rotation");
  parameters.emplace_back("parameters");
  dfs.emplace_back("df");
  taus.emplace_back("tau");

  std::stringstream params_ss;   // 2 decimal places
  std::stringstream dfs_ss;      // 1 decimal place
  std::stringstream taus_ss;     // 2 decimal places
  std::stringstream rotation_ss; // 0 decimal places
  params_ss << std::fixed << std::setprecision(2);
  dfs_ss << std::fixed << std::setprecision(1);
  taus_ss << std::fixed << std::setprecision(2);
  rotation_ss << std::fixed << std::setprecision(0);
  for (size_t t : trees_to_summarize) {
    for (size_t e = 0; e < d_ - 1 - t; ++e) {
      trees_s.push_back(std::to_string(t + 1));
      edges.push_back(std::to_string(e + 1));
      conditioned_variables.push_back(std::to_string(order[e]) + ", " +
                                      std::to_string(arr(t, e)));
      if (t > 0) {
        std::string cv;
        for (size_t cv_ = t - 1; cv_ > 0; --cv_) {
          cv += std::to_string(arr(cv_, e)) + ", ";
        }
        cv += std::to_string(arr(0, e));
        conditioning_variables.push_back(cv);
      } else {
        conditioning_variables.emplace_back("");
      }
      var_types.push_back(var_types_[order[e] - 1] + ", " +
                          var_types_[arr(t, e) - 1]);

      if (t < pair_copulas_.size()) {
        params_ss.str("");
        dfs_ss.str("");
        taus_ss.str("");
        rotation_ss.str("");
        auto bicop = pair_copulas_[t][e];
        families.push_back(bicop.get_family_name());
        if (bicop.get_family() == BicopFamily::tll) {
          params_ss << "[30x30 grid]";
          dfs_ss << bicop.get_npars();
        } else if (bicop.get_family() != BicopFamily::indep) {
          // They are concatenated on a single row with ", " as a separator
          auto bicop_params = bicop.get_parameters();
          for (long int row = 0; row < bicop_params.rows(); ++row) {
            for (long int col = 0; col < bicop_params.cols(); ++col) {
              params_ss << bicop_params(row, col);
              // Add separator if not last element
              if (col < bicop_params.cols() - 1 ||
                  row < bicop_params.rows() - 1) {
                params_ss << ", ";
              }
            }
          }
          rotation_ss << bicop.get_rotation();
          dfs_ss << bicop.get_npars();
        }

        taus_ss << bicop.get_tau();
        parameters.push_back(params_ss.str());
        dfs.push_back(dfs_ss.str());
        taus.push_back(taus_ss.str());
        rotations.push_back(rotation_ss.str());

      } else {
        families.push_back(get_family_name(BicopFamily::indep));
        rotations.emplace_back("");
        parameters.emplace_back("");
        taus.emplace_back("0.0");
        dfs.emplace_back("");
      }
    }
  }

  std::vector<std::vector<std::string>> vinecop_str_vec = {
    trees_s,
    edges,
    conditioned_variables,
    conditioning_variables,
    var_types,
    families,
    rotations,
    parameters,
    dfs,
    taus
  };

  vinecop_str << tools_stl::dataframe_to_string(vinecop_str_vec).str();
  return vinecop_str.str();
}
}
