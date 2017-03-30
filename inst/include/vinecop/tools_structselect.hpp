// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include <cstddef>
#include <boost/graph/adjacency_list.hpp>
#include "bicop/class.hpp"
#include "fit_controls.hpp"

// to allow for (auto e : boost::edges(g)) notation
namespace std
{
    template <class T>
    T begin(const std::pair<T,T>& eItPair) { return eItPair.first; }

    template <class T>
    T end(const std::pair<T,T>& eItPair) { return eItPair.second; }
}

namespace tools_structselect {
    // boost::graph represenation of a vine tree ----------------
    struct VertexProperties {
        std::vector<size_t> conditioning;
        std::vector<size_t> conditioned;
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
        vinecopulib::Bicop pair_copula;
    };
    typedef boost::adjacency_list <
        boost::vecS,
        boost::vecS,
        boost::undirectedS,
        VertexProperties,
        boost::property<boost::edge_weight_t, double, EdgeProperties>
    > VineTree;


    // functions for manipulation of trees ----------------
    VineTree make_base_tree(const Eigen::MatrixXd& data);
    VineTree select_next_tree(VineTree& prev_tree,
                              vinecopulib::FitControlsVinecop& controls);
    VineTree edges_as_vertices(const VineTree& prev_tree);
    void add_allowed_edges(VineTree& tree, std::string tree_criterion,
                           double threshold);
    double get_tree_criterion(Eigen::Matrix<double, Eigen::Dynamic, 2> data,
                              std::string tree_criterion, double threshold);
    ptrdiff_t find_common_neighbor(size_t v0, size_t v1, const VineTree& tree);
    Eigen::MatrixXd get_pc_data(size_t v0, size_t v1, const VineTree& tree);
    void min_spanning_tree(VineTree &tree);
    void add_edge_info(VineTree& tree);
    void remove_edge_data(VineTree& tree);
    void remove_vertex_data(VineTree& tree);
    void select_pair_copulas(VineTree& tree,
                             vinecopulib::FitControlsVinecop& controls);

    void print_pair_copulas(VineTree& tree);
    std::string get_pc_index(boost::graph_traits<VineTree>::edge_descriptor e,
                             VineTree& tree);

}
