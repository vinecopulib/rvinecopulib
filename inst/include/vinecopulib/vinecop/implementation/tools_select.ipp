// Copyright © 2018 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/misc/tools_stl.hpp>
#include <vinecopulib/misc/tools_stats.hpp>
#include <vinecopulib/misc/tools_parallel.hpp>

#include <cmath>
#include <boost/graph/prim_minimum_spanning_tree.hpp>

#include <iostream>

namespace vinecopulib {

namespace tools_select {

using namespace tools_stl;

//! Calculate criterion for tree selection
//! @param data observations.
//! @param tree_criterion the criterion.
inline double calculate_criterion(const Eigen::Matrix<double, Eigen::Dynamic, 2>& data,
                                  std::string tree_criterion)
{
    double w = 0.0;
    Eigen::Matrix<double, Eigen::Dynamic, 2> data_no_nan =
        tools_eigen::nan_omit(data);
    double freq = static_cast<double>(data_no_nan.rows()) / static_cast<double>(data.rows());
    if (data_no_nan.rows() > 10) {
        if (tree_criterion == "tau") {
            w = tools_stats::pairwise_tau(data_no_nan);
        } else if (tree_criterion == "hoeffd") {
            w = tools_stats::pairwise_hoeffd(data_no_nan);
        } else if (tree_criterion == "rho") {
            w = tools_stats::pairwise_cor(data_no_nan);
        }
    }

    return std::fabs(w) * std::sqrt(freq);
}

//! Calculates maximal criterion for tree selection.
//! @param data observations.
//! @param tree_criterion the criterion.
inline Eigen::MatrixXd calculate_criterion_matrix(const Eigen::MatrixXd &data,
                                                  std::string tree_criterion)
{
    size_t n = data.rows();
    size_t d = data.cols();
    Eigen::MatrixXd mat(d, d);
    mat.diagonal() = Eigen::VectorXd::Constant(d, 1.0);
    Eigen::Matrix<double, Eigen::Dynamic, 2> pair_data(n, 2);
    for (size_t i = 1; i < d; ++i) {
        for (size_t j = 0; j < i; ++j) {
            pair_data.col(0) = data.col(i);
            pair_data.col(1) = data.col(j);
            mat(i, j) = calculate_criterion(pair_data, tree_criterion);
            mat(j, i) = mat(i, j);
        }
    }
    return mat;
}

// needs to be defined
inline VinecopSelector::~VinecopSelector()
{
}

inline std::vector <std::vector<Bicop>>
VinecopSelector::get_pair_copulas() const
{
    return pair_copulas_;
}

inline RVineMatrix VinecopSelector::get_rvine_matrix() const
{
    return vine_matrix_;
}

//! Initialize object for storing pair copulas
//!
//! @param d dimension of the vine copula.
//! @param truncation_level a truncation level (optional).
//! @return A nested vector such that `pc_store[t][e]` contains a Bicop.
//!     object for the pair copula corresponding to tree `t` and edge `e`.
inline std::vector <std::vector<Bicop>> VinecopSelector::make_pair_copula_store(
    size_t d,
    size_t truncation_level)
{
    if (d < 2) {
        throw std::runtime_error("the dimension should be larger than 1");
    }

    size_t n_trees = std::min(d - 1, truncation_level);
    std::vector <std::vector<Bicop>> pc_store(n_trees);
    for (size_t t = 0; t < n_trees; ++t) {
        pc_store[t].resize(d - 1 - t);
    }

    return pc_store;
}

inline void VinecopSelector::select_all_trees(const Eigen::MatrixXd &data)
{
    initialize_new_fit(data);
    for (size_t t = 0; t < d_ - 1; ++t) {
        select_tree(t);  // select pair copulas (+ structure) of tree t

        if (controls_.get_show_trace()) {
            std::stringstream tree_heading;
            std::cout << "** Tree: " << t << std::endl;
            print_pair_copulas_of_tree(t);
        }

        if (controls_.get_truncation_level() == t + 1) {
            // don't need to fit the remaining trees
            break;
        }
    }
    finalize(controls_.get_truncation_level());
}

inline void
VinecopSelector::sparse_select_all_trees(const Eigen::MatrixXd &data)
{
    // family set must be reset after each iteration of the threshold search
    auto family_set = controls_.get_family_set();

    std::vector<double> thresholded_crits;
    if (controls_.get_select_threshold()) {
        // initialize threshold with maximum pairwise |tau| (all pairs get
        // thresholded)
        auto tree_crit = controls_.get_tree_criterion();
        auto pairwise_crits = calculate_criterion_matrix(data, tree_crit);
        for (size_t i = 1; i < d_; ++i) {
            for (size_t j = 0; j < i; ++j) {
                thresholded_crits.push_back(pairwise_crits(i, j));
            }
        }
        // this is suboptimal for fixed structures, because several
        // iterations have to be run before first non-thresholded copula
        // appears.
    }

    double mbicv_opt = 0.0;
    bool needs_break = false;
    while (!needs_break) {
        // restore family set in case previous threshold iteration also
        // truncated the model
        controls_.set_family_set(family_set);
        controls_.set_truncation_level(std::numeric_limits<size_t>::max());
        initialize_new_fit(data);

        // decrease the threshold
        if (controls_.get_select_threshold()) {
            controls_.set_threshold(get_next_threshold(thresholded_crits));
            if (controls_.get_show_trace()) {
                std::cout <<
                    "***** threshold: " <<
                    controls_.get_threshold() <<
                    std::endl;
            }
        }

        // helper variables for checking whether an optimum was found
        double mbicv = 0.0;
        double mbicv_trunc = 0.0;
        
        for (size_t t = 0; t < d_ - 1; ++t) {
            if (controls_.get_truncation_level() < t) {
                break;  // don't need to fit the remaining trees
            }

            // select pair copulas (and possibly tree structure)
            select_tree(t);

            // update fit statistic
            mbicv_trunc += get_mbicv_of_tree(t);

            // print trace for this tree level
            if (controls_.get_show_trace()) {
                std::cout << "** Tree: " << t;
                if (controls_.get_select_truncation_level()) {
                    std::cout << ", mbicv: " << mbicv_trunc;
                }
                std::cout << std::endl;
                // print fitted pair-copulas for this tree
                print_pair_copulas_of_tree(t);
            }
            
            // mbicv comparison
            if (controls_.get_select_truncation_level() & (mbicv_trunc >= mbicv)) {
                // mbicv did not improve, truncate
                controls_.set_truncation_level(t);
                set_tree_to_indep(t);
                set_current_fit_as_opt();
                if (!controls_.get_select_threshold()) {
                    // if threshold is fixed, no need to go further
                    needs_break = true;
                }
            } else {
                mbicv = mbicv_trunc;
            }
        }

        if (controls_.get_show_trace()) {
            std::cout << "--> mbicv = " << mbicv << std::endl << std::endl;
        }


        // check whether mbicv-optimal model has been found
        if (mbicv == 0.0) {
            set_current_fit_as_opt();
            if (!controls_.get_select_threshold()) {
                // threshold is fixed, optimal truncation level has been found
                needs_break = true;
            }
        } else if (mbicv >= mbicv_opt) {
            // old model is optimal
            needs_break = true;
        } else {
            // optimum hasn't been found
            set_current_fit_as_opt();
            mbicv_opt = mbicv;
            // while loop is only for threshold selection
            needs_break = needs_break | !controls_.get_select_threshold();
            // threshold is too close to 0
            needs_break = needs_break | (controls_.get_threshold() < 0.01);
            // prepare for possible next iteration
            thresholded_crits = get_thresholded_crits();
        }
    }
    trees_ = trees_opt_;
    finalize(controls_.get_truncation_level());
}

inline void VinecopSelector::set_tree_to_indep(size_t t)
{
    // trees_[0] is base tree, see make_base_tree()
    for (auto e : boost::edges(trees_[t + 1])) {
        trees_[t + 1][e].pair_copula = Bicop();
    }
}

// extracts the current threshold value
inline double VinecopSelector::get_threshold() const
{
    return threshold_;
}


//! chooses threshold for next iteration such that at a proportion of at
//! least 2.5% of the previously thresholded pairs become non-thresholded.
inline double VinecopSelector::get_next_threshold(
    std::vector<double> &thresholded_crits)
{
    if (thresholded_crits.size() == 0) {
        return 0.0;
    }
    // sort in descending order
    std::sort(thresholded_crits.begin(), thresholded_crits.end());
    std::reverse(thresholded_crits.begin(), thresholded_crits.end());
    // pick threshold that changes at least alpha*100 % of the pair-copulas
    double alpha = 0.05;
    size_t m = thresholded_crits.size();
    return thresholded_crits[std::ceil(static_cast<double>(m) * alpha) - 1];
}


inline FamilySelector::FamilySelector(const Eigen::MatrixXd &data,
                                      const RVineMatrix &vine_matrix,
                                      const FitControlsVinecop &controls)
{
    n_ = data.rows();
    d_ = data.cols();
    trees_.resize(1);
    controls_ = controls;
    vine_matrix_ = vine_matrix;
    threshold_ = controls.get_threshold();
    psi0_ = controls.get_psi0();
}

inline StructureSelector::StructureSelector(const Eigen::MatrixXd &data,
                                            const FitControlsVinecop &controls)
{
    n_ = data.rows();
    d_ = data.cols();
    trees_.resize(1);
    controls_ = controls;
    threshold_ = controls.get_threshold();
    psi0_ = controls.get_psi0();
}

//! Add edges allowed by the proximity condition
//!
//! Also calculates the edge weight (e.g., 1-|tau| for tree_criterion =
//! "itau").
//!
//! @param vine_tree tree of a vine.
inline void StructureSelector::add_allowed_edges(VineTree &vine_tree)
{
    std::string tree_criterion = controls_.get_tree_criterion();
    double threshold = controls_.get_threshold();

    std::mutex m;
    auto add_edge = [&] (size_t v0) {
        tools_interface::check_user_interrupt(v0 % 50 == 0);
        for (size_t v1 = 0; v1 < v0; ++v1) {
            // check proximity condition: common neighbor in previous tree
            // (-1 means 'no common neighbor')
            if (find_common_neighbor(v0, v1, vine_tree) > -1) {
                auto pc_data = get_pc_data(v0, v1, vine_tree);
                double crit = calculate_criterion(pc_data, tree_criterion);
                double w = 1.0 - static_cast<double>(crit >= threshold) * crit;
                {
                    std::lock_guard<std::mutex> lk(m);
                    auto e = boost::add_edge(v0, v1, w, vine_tree).first;
                    vine_tree[e].weight = w;
                    vine_tree[e].crit = crit;
                }
            }
        }
    };
    
    tools_parallel::map_on_pool(add_edge, 
                                boost::vertices(vine_tree), 
                                controls_.get_num_threads());
}

inline void StructureSelector::finalize(size_t trunc_lvl)
{
    using namespace tools_stl;
    pair_copulas_ = make_pair_copula_store(d_, trunc_lvl);
    Eigen::Matrix <size_t, Eigen::Dynamic, Eigen::Dynamic> mat(d_, d_);
    mat.fill(0);
    std::vector <size_t> ning_set;
    // fill matrix column by column
    for (size_t col = 0; col < d_ - 1; ++col) {
        tools_interface::check_user_interrupt();
        // matrix above trunc_lvl will be filled more efficiently later
        size_t t = 
            std::max(std::min(trunc_lvl, d_ - 1 - col), static_cast<size_t>(1));
        // start with highest tree in this column
        for (auto e : boost::edges(trees_[t])) {
            // find an edge that contains a leaf
            size_t v0 = boost::source(e, trees_[t]);
            size_t v1 = boost::target(e, trees_[t]);
            size_t min_deg = std::min(boost::out_degree(v0, trees_[t]),
                                      boost::out_degree(v1, trees_[t]));
            if (min_deg > 1) {
                continue;  // not a leaf
            }
            // find position of leaf in the edge
            ptrdiff_t pos = (boost::out_degree(v1, trees_[t]) == 1);
            if (pos == 1) {
                trees_[t][e].pair_copula.flip();
            }
            // fill diagonal entry with leaf index
            mat(d_ - 1 - col, col) = trees_[t][e].conditioned[pos];
            // entry in row t-1 is other index of the edge
            mat(t - 1, col) = trees_[t][e].conditioned[std::abs(1 - pos)];

            // assign fitted pair copula to appropriate entry, see
            // `Vinecop::get_pair_copula()`.
            if (trunc_lvl > 0) {
                pair_copulas_[t - 1][col] = trees_[t][e].pair_copula;
            }

            // initialize running set with full conditioning set of this edge
            ning_set = trees_[t][e].conditioning;

            // remove edge (must not be reused in another column!)
            boost::remove_edge(v0, v1, trees_[t]);
            break;
        }

        // fill column bottom to top
        for (size_t k = 1; k < t; ++k) {
            auto check_set = cat(mat(d_ - 1 - col, col), ning_set);
            for (auto e : boost::edges(trees_[t - k])) {
                // search for an edge in lower tree that shares all
                // indices in the conditioning set + diagonal entry
                if (!is_same_set(trees_[t - k][e].all_indices, check_set)) {
                    continue;
                }
                // found suitable edge ->
                // next matrix entry is conditioned variable of new edge
                // that's not equal to the diagonal entry of this column
                auto e_new = trees_[t - k][e];
                ptrdiff_t pos = (mat(d_ - 1 - col, col) == e_new.conditioned[1]);
                if (pos == 1) {
                    e_new.pair_copula.flip();
                }
                mat(t - k - 1, col) = e_new.conditioned[std::abs(1 - pos)];

                // assign fitted pair copula to appropriate entry, see
                // Vinecop::get_pair_copula().
                pair_copulas_[t - 1 - k][col] = e_new.pair_copula;

                // start over with conditioned set of next edge
                ning_set = e_new.conditioning;

                // remove edge (must not be reused in another column!)
                size_t v0 = boost::source(e, trees_[t - k]);
                size_t v1 = boost::target(e, trees_[t - k]);
                boost::remove_edge(v0, v1, trees_[t - k]);
                break;
            }
        }
    }

    // The last column contains a single element which must be different
    // from all other diagonal elements. Based on the properties of an
    // R-vine matrix, this must be the element next to it.
    mat(0, d_ - 1) = mat(0, d_ - 2);

    // change to user-facing format
    // (variable index starting at 1 instead of 0)
    for (size_t i = 0; i < d_; ++i) {
        for (size_t j = 0; j < d_ - i; ++j) {
            mat(i, j) += 1;
        }
    }
    
    // fill missing entries in case vine was truncated
    RVineMatrix::complete_matrix(mat, trunc_lvl, controls_.get_num_threads());

    // return as RVineMatrix
    vine_matrix_ = RVineMatrix(mat, false);
}

inline bool FamilySelector::belongs_to_structure(size_t v0, size_t v1,
                                                 const VineTree &vine_tree)
{
    // -1 means no common neighbor in previous tree
    if (find_common_neighbor(v0, v1, vine_tree) > -1) {
        std::vector <size_t> conditioning;
        std::vector <size_t> conditioned(2);
        if (vine_tree[v0].conditioned.size() == 1) {
            // first tree
            conditioned = cat(vine_tree[v0].conditioned,
                              vine_tree[v1].conditioned);
            conditioning = std::vector<size_t>(0);
        } else {
            // compute new conditioned/conditioning sets
            conditioned = set_sym_diff(vine_tree[v0].all_indices,
                                       vine_tree[v1].all_indices);
            conditioning = intersect(vine_tree[v0].all_indices,
                                     vine_tree[v1].all_indices);
        }

        // to convert from vinecop to rvine_matrix indices
        auto add_one = [](std::vector <size_t> &v) {
            std::for_each(v.begin(), v.end(), [](size_t &d) { d += 1; });
        };
        add_one(conditioned);
        add_one(conditioning);

        // check whether the edge belongs to the structure
        return vine_matrix_.belongs_to_structure(conditioned, conditioning);
    }

    // there was no common neighbor
    return false;
}

//! Add edges allowed by vine matrix structure
//!
//! @param vine_tree tree of a vine.
inline void FamilySelector::add_allowed_edges(VineTree &vine_tree)
{
    double w = 1.0;
    std::string tree_criterion = controls_.get_tree_criterion();
    for (auto v0 : boost::vertices(vine_tree)) {
        tools_interface::check_user_interrupt(v0 % 10000 == 0);
        for (auto v1 : boost::vertices(vine_tree)) {
            if (v0 == v1)
                continue;
            // check whether edege (v0, v1) belongs to the structure
            // given in rvine_matrix_
            if (belongs_to_structure(v0, v1, vine_tree)) {
                Eigen::MatrixXd pc_data;
                EdgeIterator e;
                pc_data = get_pc_data(v0, v1, vine_tree);
                e = boost::add_edge(v0, v1, w, vine_tree).first;
                double crit = calculate_criterion(pc_data, tree_criterion);
                vine_tree[e].weight = w;
                vine_tree[e].crit = crit;
            }
        }
    }
}

inline void FamilySelector::finalize(size_t trunc_lvl)
{
    pair_copulas_ = make_pair_copula_store(d_, trunc_lvl);
    for (size_t tree = 0; tree < pair_copulas_.size(); tree++) {
        for (auto e : boost::edges(trees_[tree + 1])) {
            // check in which column of the matrix the pair-copula e is
            size_t edge = find_column_in_matrix(trees_[tree + 1][e].conditioned);
            // trees_[0] is base tree, vine copula starts at trees_[1]
            pair_copulas_[tree][edge] = trees_[tree + 1][e].pair_copula;
            edge++;
        }
    }
}

inline size_t FamilySelector::find_column_in_matrix(
    const std::vector<size_t>& conditioned)
{
    Eigen::Matrix<size_t, Eigen::Dynamic, 1> inds =
        vine_matrix_.get_order().reverse();
    std::vector<size_t> vinds(inds.data(), inds.data() + inds.rows());
    return tools_stl::find_position(conditioned[0] + 1, vinds);
}


//! Extract pair copula pseudo-observations from h-functions
//!
//! @param v0,v1 vertex indices.
//! @param tree a vine tree.
//! @return The pseudo-observations for the pair coula, extracted from
//!     the h-functions calculated in the previous tree.
inline Eigen::MatrixXd VinecopSelector::get_pc_data(size_t v0, size_t v1,
                                                    const VineTree &tree)
{
    Eigen::MatrixXd pc_data(tree[v0].hfunc1.size(), 2);
    size_t ei_common = find_common_neighbor(v0, v1, tree);
    if (find_position(ei_common, tree[v0].prev_edge_indices) == 0) {
        pc_data.col(0) = tree[v0].hfunc1;
    } else {
        pc_data.col(0) = tree[v0].hfunc2;
    }
    if (find_position(ei_common, tree[v1].prev_edge_indices) == 0) {
        pc_data.col(1) = tree[v1].hfunc1;
    } else {
        pc_data.col(1) = tree[v1].hfunc2;
    }

    return pc_data;
}

//! Select and fit next tree of the vine
//!
//! The next tree is found the following way:
//!     1. Edges of the previous tree become edges in the new tree.
//!     2. All edges allowed by the proximity condition are added to the
//!        new graph.
//!     3. Collapse the new graph to a maximum spanning tree for edge
//!        weight.
//!     4. Populate edges with conditioned/conditioning sets and pseudo-
//!        observations.
//!     5. Fit and select a copula model for each edge.
//!
//! @param prev_tree tree T_{k}.
//! @param controls the controls for fitting a vine copula
//!     (see FitControlsVinecop).
//! @param tree_opt the current optimal tree (used only for sparse
//!     selection).
inline void VinecopSelector::select_tree(size_t t)
{
    auto new_tree = edges_as_vertices(trees_[t]);
    remove_edge_data(trees_[t]); // no longer needed
    add_allowed_edges(new_tree);
    if (boost::num_vertices(new_tree) > 2) {
        // has no effect in FamilySelector
        min_spanning_tree(new_tree);
    }
    add_edge_info(new_tree);       // for pc estimation and next tree
    remove_vertex_data(new_tree);  // no longer needed
    if (controls_.get_selection_criterion() == "mbicv") {
        // adjust prior probability to tree level
        controls_.set_psi0(std::pow(psi0_, t + 1));
    }
    if (trees_opt_.size() > t + 1) {
        select_pair_copulas(new_tree, trees_opt_[t + 1]);
    } else {
        select_pair_copulas(new_tree);
    }
    // make sure there is space for new tree
    trees_.resize(t + 2);
    trees_[t + 1] = new_tree;
}

inline double VinecopSelector::get_mbicv_of_tree(size_t t)
{
    double loglik = get_loglik_of_tree(t);
    double npars = get_npars_of_tree(t);
    size_t non_indeps = get_num_non_indeps_of_tree(t);
    size_t indeps = d_ - t - 1 - non_indeps;
    double psi0 = std::pow(psi0_, t + 1);
    double log_prior = 
        static_cast<double>(non_indeps) * std::log(psi0) +
        static_cast<double>(indeps) * std::log(1.0 - psi0);
    return -2 * loglik + std::log(n_) * npars - 2 * log_prior;    
}


//! calculates the log-likelihood of a tree.
inline double VinecopSelector::get_loglik_of_tree(size_t t)
{
    double ll = 0.0;
    // trees_[0] is base tree, see make_base_tree()
    for (const auto &e : boost::edges(trees_[t + 1])) {
        ll += trees_[t + 1][e].loglik;
    }
    return ll;
}

//! calculates the numbers of parameters of a tree.
inline double VinecopSelector::get_npars_of_tree(size_t t)
{
    double npars = 0.0;
    // trees_[0] is base tree, see make_base_tree()
    for (const auto &e : boost::edges(trees_[t + 1])) {
        npars += trees_[t + 1][e].npars;
    }
    return npars;
}

//! calculates the numbers of independence copulas in a tree.
inline size_t VinecopSelector::get_num_non_indeps_of_tree(size_t t)
{
    size_t num_non_indeps = 0;
    // trees_[0] is base tree, see make_base_tree()
    for (const auto &e : boost::edges(trees_[t + 1])) {
        num_non_indeps += static_cast<double>(trees_[t + 1][e].npars > 0);
    }
    return num_non_indeps;
}


//! Print indices, family, and parameters for each pair-copula
//! @param tree a vine tree.
inline void VinecopSelector::print_pair_copulas_of_tree(size_t t)
{
    // trees_[0] is the base tree, see make_base_tree()
    for (auto e : boost::edges(trees_[t + 1])) {
        std::cout << 
            get_pc_index(e, trees_[t + 1]) << " <-> " <<
            trees_[t + 1][e].pair_copula.str() << std::endl;
    }
}

//! extracts all criterion values that got thresholded to zero.
inline std::vector<double> VinecopSelector::get_thresholded_crits()
{
    std::vector<double> crits;
    for (size_t t = 1; t < trees_.size(); ++t) {
        for (auto e : boost::edges(trees_[t])) {
            if (trees_[t][e].crit < controls_.get_threshold()) {
                crits.push_back(trees_[t][e].crit);
            }
        }
    }

    return crits;
}

inline void VinecopSelector::initialize_new_fit(const Eigen::MatrixXd &data)
{
    trees_[0] = make_base_tree(data);
}

inline void VinecopSelector::set_current_fit_as_opt()
{
    threshold_ = controls_.get_threshold();
    trees_opt_ = trees_;
}

//! Create base tree of the vine
//!
//!  The base tree is a star on d + 1 variables, where the conditioned
//!  set of each edge consists of a single number. When building the next
//!  tree, the edges become vertices. Because the base graph was a star
//!  all edges are allowed by the proximity condition, and the edges will
//!  have a conditioned set consisting of the two vertex indices. This
//!  will be the first actual tree of the vine.
//!
//!  @param data nxd matrix of copula data.
//!  @return A VineTree object containing the base graph.
inline VineTree VinecopSelector::make_base_tree(const Eigen::MatrixXd &data)
{
    size_t d = data.cols();
    VineTree base_tree(d);
    // a star connects the root node (d) with all other nodes
    for (size_t target = 0; target < d; ++target) {
        tools_interface::check_user_interrupt(target % 10000 == 0);
        // add edge and extract edge iterator
        auto e = add_edge(d, target, base_tree).first;

        // inititialize hfunc1 with actual data for variable "target"
        base_tree[e].hfunc1 = data.col(boost::target(e, base_tree));
        // identify edge with variable "target" and initialize sets
        base_tree[e].conditioned.reserve(2);
        base_tree[e].conditioned.push_back(boost::target(e, base_tree));
        base_tree[e].conditioning.reserve(d - 2);
        base_tree[e].all_indices = base_tree[e].conditioned;
    }

    return base_tree;
}

//! Convert edge set into vertex set of a new graph
//!
//! Further information about the structure is passed along:
//!     - conditioned/conditioning set,
//!     - indices of vertices connected by the edge in the previous tree.
//!
//! @param tree T_{k}.
//! @return A edge-less graph of vertices, each representing one edge of the
//!     previous tree.
inline VineTree VinecopSelector::edges_as_vertices(const VineTree &prev_tree)
{
    // start with full graph
    size_t d = num_edges(prev_tree);
    VineTree new_tree(d);

    // cut & paste information from previous tree
    int i = 0;
    for (auto e : boost::edges(prev_tree)) {
        new_tree[i].hfunc1 = prev_tree[e].hfunc1;
        new_tree[i].hfunc2 = prev_tree[e].hfunc2;
        new_tree[i].conditioned = prev_tree[e].conditioned;
        new_tree[i].conditioning = prev_tree[e].conditioning;
        new_tree[i].all_indices = prev_tree[e].all_indices;
        new_tree[i].prev_edge_indices.reserve(2);
        new_tree[i].prev_edge_indices.push_back(boost::source(e, prev_tree));
        new_tree[i].prev_edge_indices.push_back(boost::target(e, prev_tree));
        ++i;
    }

    return new_tree;
}

// Find common neighbor in previous tree
//
// @param v0,v1 vertices in the tree.
// @param tree the current tree.
// @return Gives the index of the vertex in the previous tree that was
//     shared by e0, e1, the edge representations of v0, v1.
inline ptrdiff_t VinecopSelector::find_common_neighbor(size_t v0, size_t v1,
                                                       const VineTree &tree)
{
    auto ei0 = tree[v0].prev_edge_indices;
    auto ei1 = tree[v1].prev_edge_indices;
    auto ei_common = intersect(ei0, ei1);

    if (ei_common.size() == 0) {
        return -1;
    } else {
        return ei_common[0];
    }
}

//! Collapse a graph to the minimum spanning tree
//!
//! @param graph the input graph.
//! @return the input graph with all non-MST edges removed.
inline void VinecopSelector::min_spanning_tree(VineTree &graph)
{
    size_t d = num_vertices(graph);
    std::vector <size_t> targets(d);
    prim_minimum_spanning_tree(graph, targets.data());
    for (size_t v1 = 0; v1 < d; ++v1) {
        for (size_t v2 = 0; v2 < v1; ++v2) {
            if ((v2 != targets[v1]) & (v1 != targets[v2])) {
                boost::remove_edge(v1, v2, graph);
            }
        }
    }
}

//! Add conditioned info and data for each edge
//!
//! See, e.g., Czado (2010), "Pair-copula constructions of multivariate
//! copulas", url: https://mediatum.ub.tum.de/doc/1079253/file.pdf
//! @param tree a vine tree.
inline void VinecopSelector::add_edge_info(VineTree &tree)
{
    for (auto e : boost::edges(tree)) {
        auto v0 = boost::source(e, tree);
        auto v1 = boost::target(e, tree);
        tree[e].pc_data = get_pc_data(v0, v1, tree);
        tree[e].conditioned = set_sym_diff(tree[v0].all_indices,
                                           tree[v1].all_indices);
        tree[e].conditioning = intersect(tree[v0].all_indices,
                                         tree[v1].all_indices);
        tree[e].all_indices = cat(tree[e].conditioned, tree[e].conditioning);
    }
}

//! Remove data (hfunc1/hfunc2/pc_data) from all edges of a vine tree
//! @param tree a vine tree.
inline void VinecopSelector::remove_edge_data(VineTree &tree)
{
    for (auto e : boost::edges(tree)) {
        tree[e].hfunc1 = Eigen::VectorXd();
        tree[e].hfunc2 = Eigen::VectorXd();
        tree[e].pc_data = Eigen::Matrix<double, Eigen::Dynamic, 2>(0, 2);
    }
}

//! Remove data (hfunc1/hfunc2/pc_data) from all vertices of a vine tree
//! @param tree a vine tree.
inline void VinecopSelector::remove_vertex_data(VineTree &tree)
{
    for (auto v : boost::vertices(tree)) {
        tree[v].hfunc1 = Eigen::VectorXd();
        tree[v].hfunc2 = Eigen::VectorXd();
    }
}

//! Fit and select a pair copula for each edges
//! @param tree a vine tree preprocessed with add_edge_info().
//! @param tree_opt the current optimal tree (used only for sparse
//!     selection).
inline void VinecopSelector::select_pair_copulas(VineTree &tree,
                                                 const VineTree &tree_opt)
{
    auto select_pc = [&](EdgeIterator e) -> void {
        tools_interface::check_user_interrupt();
        bool is_thresholded = (tree[e].crit < controls_.get_threshold());
        bool used_old_fit = false;

        if (controls_.needs_sparse_select()) {
            // the formula is quite arbitrary, but sufficient for
            // identifying situations where fits can be re-used
            tree[e].fit_id =
                (tree[e].pc_data.col(0) - 2 * tree[e].pc_data.col(1)).sum();
            tree[e].fit_id += 5.0 * static_cast<double>(is_thresholded);
            if (boost::num_edges(tree_opt) > 0) {
                auto old_fit = find_old_fit(tree[e].fit_id, tree_opt);
                if (old_fit.second) {  // indicates if match was found
                    // data and thresholding status haven't changed,
                    // we can use old fit
                    used_old_fit = true;
                    tree[e].pair_copula = tree_opt[old_fit.first].pair_copula;
                }
            }
        }

        if (!used_old_fit) {
            if (is_thresholded) {
                tree[e].pair_copula = vinecopulib::Bicop();
            } else {
                tree[e].pair_copula.select(tree[e].pc_data, controls_);
            }
        }

        tree[e].hfunc1 = tree[e].pair_copula.hfunc1(tree[e].pc_data);
        tree[e].hfunc2 = tree[e].pair_copula.hfunc2(tree[e].pc_data);
        if (controls_.needs_sparse_select()) {
            tree[e].loglik = tree[e].pair_copula.loglik(tree[e].pc_data);
            tree[e].npars = tree[e].pair_copula.calculate_npars();
        }
    };
    
    // make sure that Bicop.select() doesn't spawn new threads
    size_t num_threads = controls_.get_num_threads();
    controls_.set_num_threads(1);
    tools_parallel::map_on_pool(select_pc, 
                                boost::edges(tree), 
                                num_threads);
    controls_.set_num_threads(num_threads);
}

//! finds the fitted pair-copula from the previous iteration.
inline FoundEdge VinecopSelector::find_old_fit(double fit_id,
                                               const VineTree &old_graph)
{
    auto edge = boost::edge(0, 1, old_graph).first;
    bool fit_with_same_id = false;
    for (auto e : boost::edges(old_graph)) {
        if (fit_id == old_graph[e].fit_id) {
            fit_with_same_id = true;
            edge = e;
        }
    }
    return std::make_pair(edge, fit_with_same_id);
}

//! Get edge index for the vine (like 1, 2; 3)
//! @param e a descriptor for the edge.
//! @param tree a vine tree.
inline std::string VinecopSelector::get_pc_index(const EdgeIterator &e,
                                                 const VineTree &tree)
{
    std::stringstream index;
    // add 1 everywhere for user-facing representation (boost::graph
    // starts at 0)
    index <<
          tree[e].conditioned[0] + 1 <<
          "," <<
          tree[e].conditioned[1] + 1;
    if (tree[e].conditioning.size() > 0) {
        index << " | ";
        for (unsigned int i = 0; i < tree[e].conditioning.size(); ++i) {
            index << tree[e].conditioning[i] + 1;
            if (i < tree[e].conditioning.size() - 1)
                index << ",";
        }
    }

    return index.str().c_str();
}

}

}
