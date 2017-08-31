// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/vinecop/tools_select.hpp>
#include <vinecopulib/misc/tools_stl.hpp>
#include <vinecopulib/misc/tools_stats.hpp>
#include <vinecopulib/vinecop/class.hpp>

#include <cmath>
#include <boost/graph/prim_minimum_spanning_tree.hpp>

#include <iostream>
namespace vinecopulib {
    
namespace tools_select {

    using namespace tools_stl;
    
    //! Calculate criterion for tree selection
    //! @param data observations.
    //! @param tree_criterion the criterion.
    double calculate_criterion(Eigen::Matrix<double, Eigen::Dynamic, 2> data,
                               std::string tree_criterion)
    {
        double w = 0.0;
        if (tree_criterion == "tau") {
            w = std::fabs(tools_stats::pairwise_tau(data));
        } else if (tree_criterion == "hoeffd") {
            // scale to [0,1]
            w = (30 * tools_stats::pairwise_hoeffd(data) + 0.5) / 1.5;
        } else if (tree_criterion == "rho") {
            w = std::fabs(tools_stats::pairwise_cor(data));
        }
        
        return w;
    }
    
    //! Calculates maximal criterion for tree selection.
    //! @param data observations.
    //! @param tree_criterion the criterion.
    Eigen::MatrixXd calculate_criterion_matrix(const Eigen::MatrixXd& data, 
                                               std::string tree_criterion)
    {
        Eigen::MatrixXd w = tools_stats::dependence_matrix(data, tree_criterion);
        if (tree_criterion == "hoeffd") {
            // hoeefd needs to be scaled to [0,1]
            w = (w.array() + 0.5) / 1.5;
        }
    
        return w.array().abs();
    }
    
    //! calculates the Generalized Information Criterion.
    double calculate_gic(double loglik, double npars, int n)
    {
        double log_npars = std::log(npars);
        if (npars == 0.0) {
            log_npars = 0.0;
        }
        return -2 * loglik + std::log(std::log((double) n)) * log_npars * npars;
    }
    
    // needs to be defined
    VinecopSelector::~VinecopSelector() {}
    
    std::vector<std::vector<Bicop>> VinecopSelector::get_pair_copulas() const
    {
            return pair_copulas_;
    }
    
    RVineMatrix VinecopSelector::get_rvine_matrix() const
    {
        return vine_matrix_;
    }
    
    //! truncates the vine copula at the current tree level (tracked by 
    //! `trees_fitted_` data member).
    void VinecopSelector::truncate() 
    {
        controls_.set_family_set({BicopFamily::indep});
    }
    
    void VinecopSelector::select_all_trees(const Eigen::MatrixXd& data)
    {
        initialize_new_fit(data);
        for (size_t t = 0; t < d_ - 1; ++t) {
            select_tree(t);  // select pair copulas (+ structure) of tree t
        
            if (controls_.get_show_trace()) {
                std::stringstream tree_heading;
                tree_heading << "** Tree: " << t << std::endl;
                tools_interface::print(tree_heading.str().c_str());
                print_pair_copulas_of_tree(t);
            }
            
            if (controls_.get_truncation_level() == t + 1) {
                truncate();  // only allow for Independence copula from here on
            }
        }
        finalize();
    }
    
    
    void VinecopSelector::sparse_select_all_trees(const Eigen::MatrixXd& data)
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
        
        double gic_opt = 0.0;
        bool needs_break = false;
        std::stringstream msg;
        while (!needs_break) {     
            // restore family set in case previous threshold iteration also 
            // truncated the model
            controls_.set_family_set(family_set);  
            initialize_new_fit(data);
            
            // decrease the threshold
            if (controls_.get_select_threshold()) {
                controls_.set_threshold(get_next_threshold(thresholded_crits));
                if (controls_.get_show_trace()) {
                    msg << 
                        "***** threshold: " << 
                        controls_.get_threshold() <<
                        std::endl;
                    tools_interface::print(msg.str().c_str());
                    msg.str("");  // clear stream
                }
            }
            
            // helper variables for checking whether an optimum was found
            double loglik = 0.0;
            double npars = 0.0;
            double gic = 0.0;     
            double gic_trunc = 0.0;   
            
            for (size_t t = 0; t < d_ - 1; ++t) {
                if (controls_.get_show_trace()) {
                    msg << "** Tree: " << t;
                }
                
                // select pair copulas (and possibly tree structure)
                select_tree(t);
                                                
                // update fit statistics
                loglik += get_loglik_of_tree(t);
                npars += get_npars_of_tree(t);
                gic_trunc = calculate_gic(loglik, npars, n_);
                
                if (controls_.get_select_truncation_level()) {
                    if (gic_trunc >= gic) {
                        // gic did not improve, set this and all remaining trees
                        // to independence
                        set_tree_to_indep(t);
                        truncate();
                    } else {
                        gic = gic_trunc; 
                    }
                } else {
                    gic = gic_trunc; 
                }
                
                if (controls_.get_show_trace()) {
                    if (controls_.get_select_truncation_level()) {
                        msg << ", GIC: " << gic_trunc;
                    }
                    msg << std::endl;
                }
            
                if (controls_.get_truncation_level() == t) {
                    truncate();
                }
                
                // print trace for this tree level
                if (controls_.get_show_trace()) {
                    tools_interface::print(msg.str().c_str());
                    msg.str("");  // clear stream
                    // print fitted pair-copulas for this tree
                    if (controls_.get_show_trace()) {
                        tools_interface::print(msg.str().c_str());
                        msg.str("");  // clear stream
                        print_pair_copulas_of_tree(t);
                    }
                }    
            }
            
            if (controls_.get_show_trace()) {
                msg << "--> GIC = " << gic << std::endl << std::endl;  
                tools_interface::print(msg.str().c_str());
                msg.str("");  // clear stream
            }
            
            // prepare for possible next iteration
            thresholded_crits = get_thresholded_crits();        
            set_current_fit_as_opt();
    
            if (gic == 0.0) {
                // stil independence, threshold needs to be reduced further
            } else if (gic >= gic_opt) {
                // new model is optimal
                needs_break = true;
            } else {
                // optimum hasn't been found
                gic_opt = gic;
                // while loop is only for threshold selection
                needs_break = needs_break | !controls_.get_select_threshold();
                // threshold is too close to 0
                needs_break = needs_break | (controls_.get_threshold() < 0.01);
            }
        }
        
        finalize();
    }
    
    //! chooses threshold for next iteration such that at a proportion of at 
    //! least 2.5% of the previously thresholded pairs become non-thresholded.
    double VinecopSelector::get_next_threshold(
        std::vector<double>& thresholded_crits)
    {
        if (thresholded_crits.size() == 0) {
            return 0.0;
        }
        // sort in descending order
        std::sort(thresholded_crits.begin(), thresholded_crits.end());
        std::reverse(thresholded_crits.begin(), thresholded_crits.end());
        // pick threshold that changes at least alpha*100 % of the pair-copulas
        double alpha = 0.025;
        size_t m = thresholded_crits.size();
        return thresholded_crits[std::ceil(m * alpha) - 1];
    }
    
    
    FamilySelector::FamilySelector(const Eigen::MatrixXd& data,
                                         const RVineMatrix& vine_matrix,
                                         const FitControlsVinecop& controls)
    {
        n_ = data.rows();
        d_ = data.cols();
        trees_ = std::vector<VineTree>(d_);
        trees_opt_ = trees_;  // for sparse selection
        controls_ = controls;
        pair_copulas_ = Vinecop::make_pair_copula_store(d_);
        vine_matrix_ = vine_matrix;
    }

    StructureSelector::StructureSelector(const Eigen::MatrixXd& data,
                                         const FitControlsVinecop& controls)
    {
        n_ = data.rows();
        d_ = data.cols();
        trees_ = std::vector<VineTree>(d_);
        trees_opt_ = trees_;  // for sparse selection
        controls_ = controls;
        pair_copulas_ = Vinecop::make_pair_copula_store(d_);
    }

    //! Add edges allowed by the proximity condition
    //!
    //! Also calculates the edge weight (e.g., 1-|tau| for tree_criterion =
    //! "itau").
    //!
    //! @param vine_tree tree of a vine.
    void StructureSelector::add_allowed_edges(VineTree& vine_tree)
    {
        std::string tree_criterion = controls_.get_tree_criterion();
        double threshold = controls_.get_threshold();
        for (auto v0 : boost::vertices(vine_tree)) {
            tools_interface::check_user_interrupt(v0 % 10000 == 0);
            for (size_t v1 = 0; v1 < v0; ++v1) {
                // check proximity condition: common neighbor in previous tree
                // (-1 means 'no common neighbor')
                if (find_common_neighbor(v0, v1, vine_tree) > -1) {
                    auto pc_data = get_pc_data(v0, v1, vine_tree);
                    double crit = calculate_criterion(pc_data, tree_criterion);
                    double w = 1.0 - (double)(crit >= threshold) * crit;
                    auto e = boost::add_edge(v0, v1, w, vine_tree).first;
                    vine_tree[e].weight = w;
                    vine_tree[e].crit = crit;
                }
            }
        }
    }
    
    void StructureSelector::finalize()
    {
        using namespace tools_stl;
        Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> mat(d_, d_);
        mat.fill(0);
    
        for (size_t col = 0; col < d_ - 1; ++col) {
            tools_interface::check_user_interrupt();
            size_t t = d_ - 1 - col;
            // start with highest tree in this column and fill first two
            // entries by conditioned set
            auto e0 = *boost::edges(trees_[t]).first;
            mat(t, col) = trees_[t][e0].conditioned[0];
            mat(t - 1, col) = trees_[t][e0].conditioned[1];
    
            // assign fitted pair copula to appropriate entry, see
            // vinecopulib::Vinecop::get_pair_copula().
            pair_copulas_[t - 1][col] = trees_[t][e0].pair_copula;
    
            // initialize running set with full conditioing set of this edge
            auto ning_set = trees_[t][e0].conditioning;
    
            // iteratively search for an edge in lower tree that shares all 
            // indices in the conditioned set + diagonal entry
            for (size_t k = 1; k < t; ++k) {
                auto reduced_set = cat(mat(t, col), ning_set);
                for (auto e : boost::edges(trees_[t - k])) {
                    if (is_same_set(trees_[t - k][e].all_indices, reduced_set)) {
                        // next matrix entry is conditioned variable of new edge
                        // that's not equal to the diagonal entry of this column
                        auto e_new = trees_[t - k][e];
                        ptrdiff_t pos = 
                            find_position(mat(t, col), e_new.conditioned);
                        mat(t - k - 1, col) = 
                            e_new.conditioned[std::abs(1 - pos)];
                        if (pos == 1) {
                            e_new.pair_copula.flip();
                        }
                        // assign fitted pair copula to appropriate entry, see
                        // vinecopulib::Vinecop::get_pair_copula().
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
        }
    
        // The last column contains a single element which must be different
        // from all other diagonal elements. Based on the properties of an
        // R-vine matrix, this must be the element next to it.
        mat(0, d_ - 1) = mat(0, d_ - 2);
    
        // change to user-facing format
        // (variable index starting at 1 instead of 0)
        auto new_mat = mat;
        for (size_t i = 0; i < d_; ++i) {
            for (size_t j = 0; j < d_ - i; ++j) {
                new_mat(i, j) += 1;    
            }
        }
        vine_matrix_ = RVineMatrix(new_mat);
    }

    bool FamilySelector::belongs_to_structure(size_t v0, size_t v1,
                                              const VineTree& vine_tree) {
        // -1 means no common neighbor in previous tree
        if (find_common_neighbor(v0, v1, vine_tree) > -1) {
            std::vector<size_t> conditioning;
            std::vector<size_t> conditioned(2);
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
            auto add_one = [] (std::vector<size_t>& v) {
                std::for_each(v.begin(), v.end(), [](size_t& d) { d+=1;}); };
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
    void FamilySelector::add_allowed_edges(VineTree& vine_tree)
    {
        double w = 1.0;
        std::string tree_criterion = controls_.get_tree_criterion();
        for (auto v0 : boost::vertices(vine_tree)) {
            tools_interface::check_user_interrupt(v0 % 10000 == 0);
            for (auto v1 : boost::vertices(vine_tree)) {
                if (v0 == v1) continue;
                // check whether edege (v0, v1) belongs to the structure
                // given in rvine_matrix_
                if (belongs_to_structure(v0, v1, vine_tree)) {
                    Eigen::MatrixXd pc_data;
                    boost::graph_traits<VineTree>::edge_descriptor e;
                    pc_data = get_pc_data(v0, v1, vine_tree);
                    e = boost::add_edge(v0, v1, w, vine_tree).first;
                    double crit = calculate_criterion(pc_data, tree_criterion);
                    vine_tree[e].weight = w;
                    vine_tree[e].crit = crit;
                }
            }
        }
    }
    
    void FamilySelector::finalize()
    {
        for (size_t tree = 0; tree < d_ - 1; tree++) {
            int edge = 0;
            // trees_[0] is base tree, vine copula starts at trees_[1]
            for (auto e : boost::edges(trees_[tree + 1])) {
                pair_copulas_[tree][edge] = trees_[tree + 1][e].pair_copula;
                edge++;
            }
        }
    }


    //! Extract pair copula pseudo-observations from h-functions
    //!
    //! @param v0,v1 vertex indices.
    //! @param tree a vine tree.
    //! @return The pseudo-observations for the pair coula, extracted from
    //!     the h-functions calculated in the previous tree.
    Eigen::MatrixXd VinecopSelector::get_pc_data(size_t v0, size_t v1,
                                                 const VineTree& tree)
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
    void VinecopSelector::select_tree(size_t t)
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
        select_pair_copulas(new_tree, trees_opt_[t + 1]);
    
        trees_[t + 1] = new_tree;
    }
        
    //! calculates the log-likelihood of a tree.
    double VinecopSelector::get_loglik_of_tree(size_t t)
    {
        double ll = 0.0;
        // trees_[0] is base tree, see make_base_tree()
        for (const auto& e : boost::edges(trees_[t + 1])) {
            ll += trees_[t + 1][e].loglik;
        }
        return ll;
    }
    
    //! calculates the numbers of parameters of a tree.
    double VinecopSelector::get_npars_of_tree(size_t t)
    {
        double npars = 0.0;
        // trees_[0] is base tree, see make_base_tree()
        for (const auto& e : boost::edges(trees_[t + 1])) {
            npars += trees_[t + 1][e].npars;
        }
        return npars;
    }
    
    void VinecopSelector::set_tree_to_indep(size_t t)
    {
        // trees_[0] is base tree, see make_base_tree()
        for (auto e : boost::edges(trees_[t + 1])) {
            trees_[t + 1][e].pair_copula = Bicop();
        }
    }
    
    //! Print indices, family, and parameters for each pair-copula
    //! @param tree a vine tree.
    void VinecopSelector::print_pair_copulas_of_tree(size_t t)
    {
        // trees_[0] is the base tree, see make_base_tree()
        for (auto e : boost::edges(trees_[t + 1])) {
            std::stringstream pc_info;
            pc_info << get_pc_index(e, trees_[t + 1]) << " <-> " <<
                        trees_[t + 1][e].pair_copula.str() << std::endl;
            vinecopulib::tools_interface::print(pc_info.str().c_str());
        }
    }
    
    //! extracts all criterion values that got thresholded to zero.
    std::vector<double> VinecopSelector::get_thresholded_crits()
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
    
    void VinecopSelector::set_current_fit_as_opt()
    {
        trees_opt_ = trees_;
    }
    
    void VinecopSelector::initialize_new_fit(const Eigen::MatrixXd& data)
    {
        trees_[0] = make_base_tree(data);
        trees_fitted_ = 0;
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
    VineTree VinecopSelector::make_base_tree(const Eigen::MatrixXd& data)
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
    VineTree VinecopSelector::edges_as_vertices(const VineTree& prev_tree)
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
    ptrdiff_t VinecopSelector::find_common_neighbor(size_t v0, size_t v1,
                                                      const VineTree& tree)
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
    void VinecopSelector::min_spanning_tree(VineTree &graph)
    {
        size_t d =  num_vertices(graph);
        std::vector<size_t> targets(d);
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
    void VinecopSelector::add_edge_info(VineTree& tree)
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
    void VinecopSelector::remove_edge_data(VineTree& tree)
    {
        for (auto e : boost::edges(tree)) {
            tree[e].hfunc1 = Eigen::VectorXd();
            tree[e].hfunc2 = Eigen::VectorXd();
            tree[e].pc_data = Eigen::Matrix<double, Eigen::Dynamic, 2>(0, 2);
        }
    }
    
    //! Remove data (hfunc1/hfunc2/pc_data) from all vertices of a vine tree
    //! @param tree a vine tree.
    void VinecopSelector::remove_vertex_data(VineTree& tree)
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
    void VinecopSelector::select_pair_copulas(VineTree& tree,
                                                const VineTree& tree_opt)
    {
        size_t count = 0;
        for (auto e : boost::edges(tree)) {
            bool is_thresholded = (tree[e].crit < controls_.get_threshold());
            bool used_old_fit = false;
            
            if (controls_.needs_sparse_select()) {
                // the formula is quite arbitrary, but sufficient for 
                // identifying situations where fits can be re-used
                tree[e].fit_id = 
                    tree[e].pc_data(0, 0) - 2 * tree[e].pc_data(0, 1); 
                tree[e].fit_id += 5.0 * (double) is_thresholded;
                if (boost::num_edges(tree_opt) > 0) {
                    auto old_fit = find_old_fit(tree[e].fit_id, tree_opt);
                    if (old_fit.second)  {  // indicates if match was found
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
                tree[e].npars  = tree[e].pair_copula.calculate_npars();
            }
        }
    }
    
    //! finds the fitted pair-copula from the previous iteration.
    FoundEdge VinecopSelector::find_old_fit(double fit_id,
                                              const VineTree& old_graph) 
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
    std::string VinecopSelector::get_pc_index(
            boost::graph_traits<VineTree>::edge_descriptor e,
            VineTree& tree
    )
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
