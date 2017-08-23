// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <cstddef>
#include <boost/graph/adjacency_list.hpp>
#include <vinecopulib/bicop/class.hpp>
#include <vinecopulib/vinecop/fit_controls.hpp>
#include <vinecopulib/vinecop/rvine_matrix.hpp>

// to allow for (auto e : boost::edges(g)) notation
namespace std
{
    template <class T>
    T begin(const std::pair<T,T>& eItPair) { return eItPair.first; }

    template <class T>
    T end(const std::pair<T,T>& eItPair) { return eItPair.second; }
}
namespace vinecopulib {

namespace tools_select {
    
    double calculate_criterion(Eigen::Matrix<double, Eigen::Dynamic, 2> data, 
                               std::string tree_criterion);
    Eigen::MatrixXd calculate_criterion_matrix(const Eigen::MatrixXd& data, 
                                               std::string tree_criterion);
    double calculate_gic(double loglik, double npars, int n);

    // boost::graph represenation of a vine tree
    struct VertexProperties {
        std::vector<size_t> conditioning;
        std::vector<size_t> conditioned;
        std::vector<size_t> all_indices;
        std::vector<size_t> prev_edge_indices;
        Eigen::VectorXd hfunc1;
        Eigen::VectorXd hfunc2;
    };
    struct EdgeProperties {
        std::vector<size_t> conditioning;
        std::vector<size_t> conditioned;
        std::vector<size_t> all_indices;
        Eigen::Matrix<double, Eigen::Dynamic, 2> pc_data;
        Eigen::VectorXd hfunc1;
        Eigen::VectorXd hfunc2;
        double weight;
        double crit;
        vinecopulib::Bicop pair_copula;
        double loglik;
        double npars;
        double fit_id;
    };
    typedef boost::adjacency_list <
            boost::vecS,
            boost::vecS,
            boost::undirectedS,
            VertexProperties,
            boost::property<boost::edge_weight_t, double, EdgeProperties>
    > VineTree;

    typedef std::pair<
            boost::graph_traits<VineTree>::edge_descriptor,
            bool
    > FoundEdge;
    
    class VinecopSelector 
    {
    public:
        virtual ~VinecopSelector() = 0;
    
        std::vector<std::vector<Bicop>> get_pair_copulas() const;
        RVineMatrix get_rvine_matrix() const;
        
        void select_all_trees(const Eigen::MatrixXd& data);
        void sparse_select_all_trees(const Eigen::MatrixXd& data);
    
    protected:
        void truncate();
        //! selects tree of vine copula, assuming all previous trees have 
        //! been fit.
        void select_tree(size_t t);
        virtual void finalize() = 0;

        double get_loglik_of_tree(size_t t);
        double get_npars_of_tree(size_t t);
        void set_tree_to_indep(size_t t);
        void print_pair_copulas_of_tree(size_t);
        std::vector<double> get_thresholded_crits();
        void set_current_fit_as_opt();
        void initialize_new_fit(const Eigen::MatrixXd& data);

        virtual void add_allowed_edges(VineTree& tree) = 0;
        Eigen::MatrixXd get_pc_data(size_t v0, size_t v1, const VineTree& tree);
        ptrdiff_t find_common_neighbor(size_t v0, size_t v1,
                                       const VineTree& tree);

    
        size_t n_;
        size_t d_;
        size_t trees_fitted_;
        FitControlsVinecop controls_;
        RVineMatrix vine_matrix_;
        std::vector<std::vector<Bicop>> pair_copulas_;
        std::vector<VineTree> trees_;
        std::vector<VineTree> trees_opt_;  // for sparse selction
    
    private:
        double get_next_threshold(std::vector<double>& thresholded_crits);
        // functions for manipulation of trees ----------------
        VineTree make_base_tree(const Eigen::MatrixXd& data);
        VineTree edges_as_vertices(const VineTree& prev_tree);
        void min_spanning_tree(VineTree &tree);
        void add_edge_info(VineTree& tree);
        void remove_edge_data(VineTree& tree);
        void remove_vertex_data(VineTree& tree);
        void select_pair_copulas(VineTree& tree, const VineTree& tree_opt);
        FoundEdge find_old_fit(double fit_id, const VineTree& old_graph);
        double get_tree_loglik(const VineTree& tree);
        double get_tree_npars(const VineTree& tree);
        std::string get_pc_index(boost::graph_traits<VineTree>::edge_descriptor e,
                                 VineTree& tree);
    };
    
    class StructureSelector : public VinecopSelector
    {
    public:
        StructureSelector(const Eigen::MatrixXd& data, 
                          const FitControlsVinecop& controls);
        ~StructureSelector() {}
        
    private:
        void add_allowed_edges(VineTree& tree);
        void finalize();
    };

    class FamilySelector : public VinecopSelector
    {
    public:
        FamilySelector(const Eigen::MatrixXd& data,
                       const RVineMatrix& vine_matrix,
                       const FitControlsVinecop& controls);
        ~FamilySelector() {}

    private:
        RVineMatrix vine_matrix_;
        void add_allowed_edges(VineTree& tree);
        bool belongs_to_structure(size_t v0, size_t v1, 
                                  const VineTree& vine_tree);
        void finalize();
    };
    
}

}
