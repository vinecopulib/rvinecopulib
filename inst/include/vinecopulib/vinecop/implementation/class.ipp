// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

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

  if (pair_copulas.size() > 0) {
    set_all_pair_copulas(pair_copulas);
  }

  if (var_types.size() > 0) {
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
    if (var_types.size() > 0) {
      d_ = var_types.size();
    } else {
      d_ = data.cols();
    }
    rvine_structure_ = RVineStructure(d_, static_cast<size_t>(0));
  }
  if (var_types.size() == 0) {
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
//! @param matrix Either an empty matrix (default) or an R-vine structure
//!     matrix, see `select()`. If empty, then it is selected as part of the
//!     fit.
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
    nobs_ = static_cast<size_t>(input["nobs_"]);
    threshold_ = static_cast<double>(input["threshold"]);
    loglik_ = static_cast<double>(input["loglik"]);
  } catch (...) {
  }
}

//! @brief Instantiates from a JSON file.
//!
//! @details The input file contains 2 attributes : `"structure"` for the vine
//! structure, which itself contains attributes `"array"` for the structure
//! triangular array and `"order"` for the order vector, and `"pair copulas"`.
//! `"pair copulas"` contains a list of attributes for the trees
//! (`"tree1"`, `"tree2"`, etc), each containing
//! a list of attributes for the edges (`"pc1"`, `"pc2"`, etc).
//! See the corresponding method of `Bicop` objects for the encoding of
//! pair-copulas.
//!
//! @param filename The name of the JSON file to read.
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

//! @brief Writes the copula object into a JSON file.
//!
//! @details The output file contains 2 attributes : `"structure"` for the vine
//! structure, which itself contains attributes `"array"` for the structure
//! triangular array and `"order"` for the order vector, and `"pair copulas"`.
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
//! @param data \f$ n \times (d + k) \f$ or \f$ n \times 2d \f$ matrix of
//!   observations, where \f$ k \f$ is the number of discrete variables.
//! @param controls The controls to the algorithm (see `FitControlsVinecop()`).
inline void
Vinecop::select(const Eigen::MatrixXd& data, const FitControlsVinecop& controls)
{
  if (controls.get_select_families()) {
    check_data(data);
    if (d_ == 1) {
      loglik_ = 0;
      nobs_ = data.rows();
      return;
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
  } else {
    fit(data, controls.get_fit_controls_bicop(), controls.get_num_threads());
  }
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
//! @param data \f$ n \times (d + k) \f$ or \f$ n \times 2d \f$ matrix of
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
    auto fit_edge = [&](size_t edge) {
      tools_interface::check_user_interrupt(edge % 5 == 0);
      // extract evaluation point from hfunction matrices (have been
      // computed in previous tree level)
      Bicop* edge_copula = &pair_copulas_[tree][edge];
      auto var_types = edge_copula->get_var_types();
      size_t m = rvine_structure_.min_array(tree, edge);

      auto u_e = Eigen::MatrixXd(n, 2), u_e_sub = Eigen::MatrixXd(n, 2);
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

      edge_copula->fit(u_e, controls);

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
}

//! @brief Automatically fits and selects a vine copula model.
//!
//! Selection of the structure is performed using the algorithm of
//! Dissmann, J. F., E. C. Brechmann, C. Czado, and D. Kurowicka (2013).
//! *Selecting and estimating regular vine copulae and application to
//! financial returns.* Computational Statistics & Data Analysis, 59 (1),
//! 52-69.
//!
//! @details When at least one variable is discrete, two types of "observations"
//! are required: the first \f$ n \times d \f$ block contains realizations of
//! \f$ F_Y(Y), F_X(X) \f$; the second \f$ n \times d \f$ block contains
//! realizations of \f$ F_Y(Y^-), F_X(X^-), ... \f$. The minus indicates a
//! left-sided limit of the cdf. For continuous variables the left limit and the
//! cdf itself coincide. For, e.g., an integer-valued variable, it holds \f$
//! F_Y(Y^-) = F_Y(Y - 1) \f$. Continuous variables in the second block can
//! be omitted.
//!
//! @param data \f$ n \times (d + k) \f$ or \f$ n \times 2d \f$ matrix of
//!   observations, where \f$ k \f$ is the number of discrete variables.
//! @param controls The controls to the algorithm (see `FitControlsVinecop()`).
inline void
Vinecop::select_all(const Eigen::MatrixXd& data,
                    const FitControlsVinecop& controls)
{
  rvine_structure_ = RVineStructure(d_, static_cast<size_t>(0));
  select(data, controls);
}

//! @brief Automatically selects all pair-copula families and fits all.
//! parameters.
//!
//!
//! @details When at least one variable is discrete, two types of "observations"
//! are required: the first \f$ n \times d \f$ block contains realizations of
//! \f$ F_Y(Y), F_X(X) \f$; the second \f$ n \times d \f$ block contains
//! realizations of \f$ F_Y(Y^-), F_X(X^-), ... \f$. The minus indicates a
//! left-sided limit of the cdf. For continuous variables the left limit and the
//! cdf itself coincide. For, e.g., an integer-valued variable, it holds \f$
//! F_Y(Y^-) = F_Y(Y - 1) \f$. Continuous variables in the second block can
//! be omitted.
//!
//! If there are missing data (i.e., NaN entries), incomplete observations are
//! discarded before fitting a pair-copula. This is done on a pair-by-pair basis
//! so that the maximal available information is used.
//!
//! @param data \f$ n \times (d + k) \f$ or \f$ n \times 2d \f$ matrix of
//!   observations, where \f$ k \f$ is the number of discrete variables.
//! @param controls The controls to the algorithm (see `FitControlsVinecop()`).
inline void
Vinecop::select_families(const Eigen::MatrixXd& data,
                         const FitControlsVinecop& controls)
{
  auto controls_trunc_lvl = controls;
  controls_trunc_lvl.set_trunc_lvl(rvine_structure_.get_trunc_lvl());
  select(data, controls);
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
//! @param natural_order Whether indices correspond to natural order.
inline TriangularArray<size_t>
Vinecop::get_struct_array(bool natural_order) const
{
  return rvine_structure_.get_struct_array(natural_order);
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
//! @param var_types A vector specifying the types of the variables,
//!   e.g., `{"c", "d"}` means first varible continuous, second discrete.
inline void
Vinecop::set_var_types(const std::vector<std::string>& var_types)
{
  check_var_types(var_types);
  set_var_types_internal(var_types);
}

//! @brief Sets all pair-copulas.
//! @param pair_copulas A vector of pair-copulas that has to be consistent with
//!     the current structure (see `Vinecop()`).
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
  for (auto t : var_types) {
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
Vinecop::set_var_types_internal(const std::vector<std::string>& var_types) const
{
  var_types_ = var_types;
  if (pair_copulas_.size() == 0) {
    return;
  }

  // set new var_types for all pair-copulas
  std::vector<std::string> natural_types(d_), pair_types(2);
  for (size_t j = 0; j < d_; ++j) {
    natural_types[j] = var_types[rvine_structure_.get_order()[j] - 1];
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
//! @param u An \f$ n \times (d + k) \f$ or \f$ n \times 2d \f$ matrix of
//!   evaluation points, where \f$ k \f$ is the number of discrete variables
//!   (see `Vinecop::select()`).
//! @param num_threads The number of threads to use for computations; if greater
//!   than 1, the function will be applied concurrently to `num_threads` batches
//!   of `u`.
//! @return A vector of length `n` containing the copula density values.
inline Eigen::VectorXd
Vinecop::pdf(Eigen::MatrixXd u, const size_t num_threads) const
{
  check_data(u);
  u = collapse_data(u);

  // info about the vine structure (reverse rows (!) for more natural indexing)
  size_t trunc_lvl = rvine_structure_.get_trunc_lvl();
  auto order = rvine_structure_.get_order();
  auto disc_cols = tools_select::get_disc_cols(var_types_);

  // initial value must be 1.0 for multiplication
  Eigen::VectorXd pdf = Eigen::VectorXd::Constant(u.rows(), 1.0);

  auto do_batch = [&](const tools_batch::Batch& b) {
    // temporary storage objects (all data must be in (0, 1))
    Eigen::MatrixXd hfunc1, hfunc2, u_e, hfunc1_sub, hfunc2_sub, u_e_sub;
    hfunc1 = Eigen::MatrixXd::Zero(b.size, d_);
    hfunc2 = Eigen::MatrixXd::Zero(b.size, d_);
    if (is_discrete()) {
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

        pdf.segment(b.begin, b.size) =
          pdf.segment(b.begin, b.size).cwiseProduct(edge_copula->pdf(u_e));

        // h-functions are only evaluated if needed in next step
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
      }
    }
  };

  if (trunc_lvl > 0) {
    tools_thread::ThreadPool pool((num_threads == 1) ? 0 : num_threads);
    pool.map(do_batch, tools_batch::create_batches(u.rows(), num_threads));
    pool.join();
  }

  return pdf;
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
//! @param u An \f$ n \times (d + k) \f$ or \f$ n \times 2d \f$ matrix of
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
             std::vector<int> seeds) const
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
  Eigen::ArrayXXd x(N, 1);
  Eigen::RowVectorXd temp(d_);
  for (size_t i = 0; i < n; i++) {
    tools_interface::check_user_interrupt(i % 1000 == 0);
    temp = u.block(i, 0, 1, d_);
    x = (u_sim.rowwise() - temp).rowwise().maxCoeff().array();
    vine_distribution(i) = static_cast<double>((x <= 0.0).count());
  }
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

//! @brief Evaluates the log-likelihood.
//!
//! @details The log-likelihood is defined as
//! \f[ \mathrm{loglik} = \sum_{i = 1}^n \log c(U_{1, i}, ..., U_{d, i}), \f]
//! where \f$ c \f$ is the copula density, see `Vinecop::pdf()`.
//!
//! @param u An \f$ n \times (d + k) \f$ or \f$ n \times 2d \f$ matrix of
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

//! @brief Evaluates the Akaike information criterion (AIC).
//!
//! @details The AIC is defined as
//! \f[ \mathrm{AIC} = -2\, \mathrm{loglik} + 2 p, \f]
//! where \f$ \mathrm{loglik} \f$ is the log-liklihood (see `Vinecop::loglik()`)
//! and \f$ p \f$ is the (effective) number of parameters of the model.
//! The AIC is a consistent model selection criterion even
//! for nonparametric models.
//!
//! @param u An \f$ n \times (d + k) \f$ or \f$ n \times 2d \f$ matrix of
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
//! @param u An \f$ n \times (d + k) \f$ or \f$ n \times 2d \f$ matrix of
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
//! @param u An \f$ n \times (d + k) \f$ or \f$ n \times 2d \f$ matrix of
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
//! \f[F(V_{M[d - j, j]} | V_{M[d - j - 1, j - 1]}, \dots, V_{M[0, 0]}), \f]
//! set `randomize_discrete = FALSE`.
//!
//! @param u An \f$ n \times d \f$ matrix of evaluation points.
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
  check_data(u);
  u = collapse_data(u);

  size_t d = d_;
  size_t n = u.rows();

  // info about the vine structure
  size_t trunc_lvl = rvine_structure_.get_trunc_lvl();
  auto order = rvine_structure_.get_order();
  auto inverse_order = tools_stl::invert_permutation(order);
  auto disc_cols = tools_select::get_disc_cols(var_types_);

  // fill first row of hfunc2 matrix with evaluation points;
  // points have to be reordered to correspond to natural order
  Eigen::MatrixXd hfunc1(n, d), hfunc2(n, d), hfunc1_sub(n, d),
    hfunc2_sub(n, d);
  for (size_t j = 0; j < d; ++j) {
    hfunc2.col(j) = u.col(order[j] - 1);
  }
  hfunc1 = hfunc2; // just ensure data is in [0, 1]^d
  if (is_discrete()) {
    hfunc1_sub = hfunc1;
    hfunc2_sub = hfunc2;
  }

  // fill first row of hfunc2 matrix with evaluation points;
  // points have to be reordered to correspond to natural order
  for (size_t j = 0; j < d_; ++j) {
    hfunc2.col(j) = u.col(order[j] - 1);
    if (var_types_[order[j] - 1] == "d") {
      hfunc2_sub.col(j) = u.col(d_ + disc_cols[order[j] - 1]);
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
        Bicop* edge_copula = &pair_copulas_[tree][edge];
        auto var_types = edge_copula->get_var_types();
        size_t m = rvine_structure_.min_array(tree, edge);

        u_e = Eigen::MatrixXd(b.size, 2);
        u_e.col(0) = hfunc2.block(b.begin, edge, b.size, 1);
        if (m == rvine_structure_.struct_array(tree, edge, true)) {
          u_e.col(1) = hfunc2.block(b.begin, m - 1, b.size, 1);
        } else {
          u_e.col(1) = hfunc1.block(b.begin, m - 1, b.size, 1);
        }

        if ((var_types[0] == "d") || (var_types[1] == "d")) {
          u_e.conservativeResize(b.size, 4);
          u_e.col(2) = hfunc2_sub.block(b.begin, edge, b.size, 1);
          if (m == rvine_structure_.struct_array(tree, edge, true)) {
            u_e.col(3) = hfunc2_sub.block(b.begin, m - 1, b.size, 1);
          } else {
            u_e.col(3) = hfunc1_sub.block(b.begin, m - 1, b.size, 1);
          }
        }

        // h-functions are only evaluated if needed in next step
        if (rvine_structure_.needed_hfunc1(tree, edge)) {
          hfunc1.block(b.begin, edge, b.size, 1) = edge_copula->hfunc1(u_e);
          if (var_types[1] == "d") {
            u_e_sub = u_e;
            u_e_sub.col(1) = u_e.col(3);
            hfunc1_sub.block(b.begin, edge, b.size, 1) =
              edge_copula->hfunc1(u_e_sub);
          }
        }
        hfunc2.block(b.begin, edge, b.size, 1) = edge_copula->hfunc2(u_e);
        if (var_types[0] == "d") {
          u_e_sub = u_e;
          u_e_sub.col(0) = u_e.col(2);
          hfunc2_sub.block(b.begin, edge, b.size, 1) =
            edge_copula->hfunc2(u_e_sub);
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
    auto R = tools_stats::simulate_uniform(u.rows(), d, false, seeds);
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
//! \f$d = 200\f$.
//!
//! The Rosenblatt transform (Rosenblatt, 1952) \f$ U = T(V) \f$ of a random
//! vector \f$ V = (V_1,\ldots,V_d) ~ F \f$ is defined as
//! \f[ U_1= F(V_1), U_{2} = F(V_{2}|V_1), \ldots, U_d
//! =F(V_d|V_1,\ldots,V_{d-1}), \f] where \f$ F(v_k|v_1,\ldots,v_{k-1}) \f$ is
//! the conditional distribution of \f$ V_k \f$ given  \f$ V_1 \ldots, V_{k-1},
//! k = 2,\ldots,d \f$. The vector \f$ U = (U_1, \dots, U_d) \f$ then contains
//! independent standard uniform variables. The inverse operation \f[V_1 =
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
//! @param u An \f$ n \times d \f$ matrix of evaluation points.
//! @param num_threads The number of threads to use for computations; if greater
//!   than 1, the function will be applied concurrently to `num_threads` batches
//!   of `u`.
//! @return An \f$ n \times d \f$ matrix of evaluations.
inline Eigen::MatrixXd
Vinecop::inverse_rosenblatt(const Eigen::MatrixXd& u,
                            const size_t num_threads) const
{
  auto var_types = get_var_types();
  set_continuous_var_types();
  check_data(u);

  size_t n = u.rows();
  size_t d = d_;

  Eigen::MatrixXd U_vine = u.leftCols(d); // output matrix
  //                   (direct + indirect)    (U_vine)       (info matrices)
  size_t bytes_required = (8 * 2 * n * d * d) + (8 * n * d) + (4 * 4 * d * d);
  // if the problem is too large (requires more than 1 GB memory), split
  // the data into two halves and call simulate on the reduced data.
  if ((n > 1) & (bytes_required > static_cast<size_t>(1e9))) {
    size_t n_half = n / 2;
    size_t n_left = n - n_half;
    U_vine.block(0, 0, n_half, d) =
      inverse_rosenblatt(u.block(0, 0, n_half, d), num_threads);
    U_vine.block(n_half, 0, n_left, d) =
      inverse_rosenblatt(u.block(n_half, 0, n_left, d), num_threads);
    return U_vine;
  }

  // info about the vine structure (in upper triangular matrix notation)
  size_t trunc_lvl = rvine_structure_.get_trunc_lvl();
  auto order = rvine_structure_.get_order();
  auto inverse_order = tools_stl::invert_permutation(order);

  auto do_batch = [&](const tools_batch::Batch& b) {
    // temporary storage objects for (inverse) h-functions
    TriangularArray<Eigen::VectorXd> hinv2(d + 1, trunc_lvl + 1);
    TriangularArray<Eigen::VectorXd> hfunc1(d + 1, trunc_lvl + 1);

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
        Bicop edge_copula = pair_copulas_[tree][var].as_continuous();

        // extract data for conditional pair
        Eigen::MatrixXd U_e(b.size, 2);
        size_t m = rvine_structure_.min_array(tree, var);
        U_e.col(0) = hinv2(tree + 1, var);
        if (m == rvine_structure_.struct_array(tree, var, true)) {
          U_e.col(1) = hinv2(tree, m - 1);
        } else {
          U_e.col(1) = hfunc1(tree, m - 1);
        }

        // inverse Rosenblatt transform simulates data for conditional pair
        hinv2(tree, var) = edge_copula.hinv2(U_e);

        // if required at later stage, also calculate hfunc2
        if (var < static_cast<ptrdiff_t>(d_) - 1) {
          if (rvine_structure_.needed_hfunc1(tree, var)) {
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

  set_var_types_internal(var_types);
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
        << "expected: " << d_exp << " or " << 2 * d_ << ", actual: " << d_data
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
//! @details The function can be const, because var_types_ is mutable.
inline void
Vinecop::set_continuous_var_types() const
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
  int n_discrete = 0;
  for (auto t : var_types_) {
    n_discrete += (t == "d");
  }
  return n_discrete;
}

inline bool
Vinecop::is_discrete() const
{
  return get_n_discrete() > 0;
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

//! @brief Summarizes the model into a string (can be used for printing).
//! @param trees A vector of tree indices to summarize; if empty, all trees.
inline std::string
Vinecop::str(const std::vector<size_t>& trees) const
{
  std::vector<size_t> trees_to_summarize;
  std::vector<size_t> all_trees(rvine_structure_.get_trunc_lvl());
  std::iota(all_trees.begin(), all_trees.end(), 0);
  if (trees.size() == 0) {
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
  trees_s.push_back("tree");
  edges.push_back("edge");
  conditioned_variables.push_back("conditioned variables");
  conditioning_variables.push_back("conditioning variables");
  var_types.push_back("var_types");
  families.push_back("family");
  rotations.push_back("rotation");
  parameters.push_back("parameters");
  dfs.push_back("df");
  taus.push_back("tau");

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
        std::string cv = "";
        for (size_t cv_ = t - 1; cv_ > 0; --cv_) {
          cv += std::to_string(arr(cv_, e)) + ", ";
        }
        cv += std::to_string(arr(0, e));
        conditioning_variables.push_back(cv);
      } else {
        conditioning_variables.push_back("");
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
        rotations.push_back("");
        parameters.push_back("");
        taus.push_back("0.0");
        dfs.push_back("");
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
