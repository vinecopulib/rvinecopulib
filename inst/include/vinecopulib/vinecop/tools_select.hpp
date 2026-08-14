// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <boost/graph/adjacency_list.hpp>
#include <vinecopulib/bicop/class.hpp>
#include <vinecopulib/misc/tools_interface.hpp>
#include <vinecopulib/vinecop/fit_controls.hpp>
#include <vinecopulib/vinecop/rvine_structure.hpp>

// to allow for (auto e : boost::edges(g)) notation
namespace std {
template<class T>
T
begin(const std::pair<T, T>& eItPair)
{
  return eItPair.first;
}

template<class T>
T
end(const std::pair<T, T>& eItPair)
{
  return eItPair.second;
}
}
namespace vinecopulib {

namespace tools_select {

double
calculate_criterion(const Eigen::MatrixXd& data,
                    const std::string& tree_criterion,
                    const Eigen::VectorXd& weights,
                    const TreeCriterionFunction& tree_criterion_function = {});

std::vector<size_t>
get_disc_cols(std::vector<std::string> var_types);

// boost::graph represenation of a vine tree
struct VertexProperties
{
  std::vector<size_t> conditioning;
  std::vector<size_t> conditioned;
  std::vector<size_t> all_indices;
  std::vector<size_t> prev_edge_indices;
  Eigen::VectorXd hfunc1;
  Eigen::VectorXd hfunc2;
  Eigen::VectorXd hfunc1_sub;
  Eigen::VectorXd hfunc2_sub;
  std::vector<std::string> var_types{ "c", "c" };
};
struct EdgeProperties
{
  std::vector<size_t> conditioning;
  std::vector<size_t> conditioned;
  std::vector<size_t> all_indices;
  Eigen::MatrixXd pc_data;
  Eigen::VectorXd hfunc1;
  Eigen::VectorXd hfunc2;
  Eigen::VectorXd hfunc1_sub;
  Eigen::VectorXd hfunc2_sub;
  std::vector<std::string> var_types{ "c", "c" };
  double weight;
  double crit;
  vinecopulib::Bicop pair_copula;
  double fit_id;
};
using VineTree = boost::adjacency_list<
  boost::vecS,
  boost::vecS,
  boost::undirectedS,
  VertexProperties,
  boost::property<boost::edge_weight_t, double, EdgeProperties>>;

using EdgeIterator = boost::graph_traits<VineTree>::edge_descriptor;
using FoundEdge = std::pair<EdgeIterator, bool>;
using WeightMap = boost::property_map<VineTree, boost::edge_weight_t>::type;

class VinecopSelector
{
public:
  VinecopSelector(const Eigen::MatrixXd& data,
                  const FitControlsVinecop& controls,
                  std::vector<std::string> var_types);

  VinecopSelector(const Eigen::MatrixXd& data,
                  const RVineStructure& vine_struct,
                  const FitControlsVinecop& controls,
                  std::vector<std::string> var_types);

  virtual ~VinecopSelector() = default;

  std::vector<std::vector<Bicop>> get_pair_copulas() const;

  RVineStructure get_rvine_structure() const;

  static std::vector<std::vector<Bicop>> make_pair_copula_store(
    size_t d,
    size_t trunc_lvl);

  void select_all_trees(const Eigen::MatrixXd& data);

  void sparse_select_all_trees(const Eigen::MatrixXd& data);

  double get_loglik() const;

  double get_threshold() const;

  size_t get_nobs() const;

  std::vector<VineTree> get_trees() const { return trees_; }
  std::vector<VineTree> get_trees_opt() const { return trees_opt_; }

protected:
  virtual void select_tree(size_t t);

  bool is_last_tree(size_t t) const;

  void finalize(size_t trunc_lvl);

  void finalize_known_structure(size_t trunc_lvl);

  void finalize_unknown_structure(size_t trunc_lvl);

  double get_mbicv_of_tree(size_t t, double loglik);

  double get_loglik_of_tree(size_t t);

  double get_npars_of_tree(size_t t);

  void set_tree_to_indep(size_t t);

  void print_pair_copulas_of_tree(size_t);

  std::vector<double> get_thresholded_crits();

  void initialize_new_fit(const Eigen::MatrixXd& data);

  void set_current_fit_as_opt(const double& loglik);

  void add_allowed_edges(VineTree& vine_tree);

  void add_allowed_edges_proximity(VineTree& vine_tree,
                                   const std::string& tree_criterion,
                                   const Eigen::VectorXd& weights,
                                   const TreeCriterionFunction& criterion_fun);

  void add_allowed_edges_structured(VineTree& vine_tree,
                                    const std::string& tree_criterion,
                                    const Eigen::VectorXd& weights,
                                    const TreeCriterionFunction& criterion_fun);

  void select_edges(VineTree& vine_tree);

  void select_edges_mst_prim(VineTree& vine_tree);

  void select_edges_mst_kruskal(VineTree& vine_tree);

  void select_edges_random(VineTree& vine_tree);

  Eigen::MatrixXd get_pc_data(size_t v0, size_t v1, const VineTree& tree);

  static const Eigen::VectorXd& get_hfunc(const VertexProperties& vertex_data,
                                          bool is_first);

  static const Eigen::VectorXd& get_hfunc_sub(
    const VertexProperties& vertex_data,
    bool is_first);

  ptrdiff_t find_common_neighbor(size_t v0, size_t v1, const VineTree& tree);

  virtual double compute_fit_id(const EdgeProperties& e);

  size_t n_;
  size_t d_;
  bool structure_unknown_{ true };
  // conditioning-aware selection: in_cond_[v] == 1 iff variable v (0-based) is
  // in the conditioning set; n_cond_ = |conditioning set|. n_cond_ == 0 means
  // ordinary (unconditional) selection and every conditioning-specific branch
  // below is skipped, leaving the default path byte-for-byte unchanged.
  std::vector<char> in_cond_;
  size_t n_cond_{ 0 };
  std::vector<std::string> var_types_;
  FitControlsVinecop controls_;
  tools_thread::ThreadPool pool_;
  std::vector<VineTree> trees_;
  RVineStructure vine_struct_;
  std::vector<std::vector<Bicop>> pair_copulas_;
  // for sparse selction
  std::vector<VineTree> trees_opt_;
  double loglik_;
  double threshold_;
  double psi0_; // initial prior probability for mbicv

  double get_next_threshold(std::vector<double>& thresholded_crits);

  // bundles the accumulators of one threshold-search pass over all trees
  struct ThresholdPass
  {
    double mbicv = 0.0;
    double mbicv_trunc = 0.0;
    double loglik = 0.0;
    double num_changed = 0.0;
    double num_total = 0.0;
    bool select_trunc_lvl = false;
    bool select_threshold = false;
  };

  ThresholdPass run_threshold_pass(bool& needs_break);

  void check_truncation_rollback(ThresholdPass& pass,
                                 size_t& t,
                                 double loglik_tree,
                                 double mbicv_tree,
                                 bool& needs_break);

  void update_optimum(const ThresholdPass& pass,
                      double& mbicv_opt,
                      bool& needs_break,
                      std::vector<double>& thresholded_crits);

  // functions for manipulation of trees ----------------
  VineTree make_base_tree(const Eigen::MatrixXd& data);

  VineTree edges_as_vertices(VineTree& prev_tree);

  void add_edge_info(VineTree& tree);

  void add_pc_info(const EdgeIterator& e, VineTree& tree);

  void remove_edge_data(VineTree& tree);

  void remove_vertex_data(VineTree& tree);

  void select_pair_copulas(VineTree& tree,
                           const VineTree& tree_opt = VineTree(),
                           bool last_tree = false);

  void fit_or_reuse_pair_copula(const EdgeIterator& e,
                                VineTree& tree,
                                const VineTree& tree_opt);

  FoundEdge find_old_fit(double fit_id, const VineTree& old_graph);

  double get_tree_loglik(const VineTree& tree);

  double get_tree_npars(const VineTree& tree);

  size_t get_num_non_indeps_of_tree(size_t t);

  std::string get_pc_index(const EdgeIterator& e, const VineTree& tree);
};
}
}

#include <vinecopulib/vinecop/implementation/tools_select.ipp>
