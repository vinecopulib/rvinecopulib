// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#include <vinecopulib/vinecop/tools_structselect.hpp>
#include <vinecopulib/misc/tools_stl.hpp>
#include <vinecopulib/misc/tools_stats.hpp>

#include <cmath>
#include <boost/graph/prim_minimum_spanning_tree.hpp>

namespace vinecopulib {
    
namespace tools_structselect {

    using namespace tools_stl;

    //! Create base tree of the vine
    //!
    //!  The base tree is a star on d + 1 variables, where the conditioning
    //!  set of each edge consists of a single number. When building the next
    //!  tree, the edges become vertices. Because the base graph was a star
    //!  all edges are allowed by the proximity condition, and the edges will
    //!  have a conditioning set consisting of the two vertex indices. This
    //!  will be the first actual tree of the vine.
    //!
    //!  @param data nxd matrix of copula data.
    //!  @return A VineTree object containing the base graph.
    VineTree make_base_tree(const Eigen::MatrixXd& data)
    {
        size_t d = data.cols();
        VineTree base_tree(d);
        // a star connects the root node (d) with all other nodes
        for (size_t target = 0; target < d; ++target) {
            // add edge and extract edge iterator
            auto e = add_edge(d, target, base_tree).first;

            // inititialize hfunc1 with actual data for variable "target"
            base_tree[e].hfunc1 = data.col(boost::target(e, base_tree));
            // identify edge with variable "target" and initialize sets
            base_tree[e].conditioning.reserve(2);
            base_tree[e].conditioning.push_back(boost::target(e, base_tree));
            base_tree[e].conditioned.reserve(d - 2);
        }

        return base_tree;
    }

    //! Select and fit next tree of the vine
    //!
    //! The next tree is found the following way:
    //!     1. Edges of the previous tree become edges in the new tree.
    //!     2. All edges allowed by the proximity condition are added to the new
    //!        graph.
    //!     3. Collapse the new graph to a maximum spanning tree for edge 
    //!        weight.
    //!     4. Populate edges with conditioning/conditioned sets and pseudo-
    //!        observations.
    //!     5. Fit and select a copula model for each edge.
    //!
    //! @param prev_tree tree T_{k}.
    //! @param controls the controls for fitting a vine copula (see FitControlsVinecop)
    //! @return tree T_{k+1}.
    VineTree select_next_tree(VineTree& prev_tree,
                              vinecopulib::FitControlsVinecop& controls)
    {
        auto new_tree = edges_as_vertices(prev_tree);
        remove_edge_data(prev_tree); // no longer needed
        add_allowed_edges(new_tree, controls.get_tree_criterion(),
                          controls.get_threshold());
        if (boost::num_vertices(new_tree) > 2) {
            min_spanning_tree(new_tree);
        }
        add_edge_info(new_tree);  // for pc estimation and next tree
        remove_vertex_data(new_tree);  // no longer needed
        select_pair_copulas(new_tree, controls);

        return new_tree;
    }

    //! Convert edge set into vertex set of a new graph
    //!
    //! Further information about the structure is passed along:
    //!     - conditioning/conditioned set,
    //!     - indices of vertices connected by the edge in the previous tree.
    //!
    //! @param tree T_{k}.
    //! @return A edge-less graph of vertices, each representing one edge of the
    //!     previous tree.
    VineTree edges_as_vertices(const VineTree& prev_tree)
    {
        // start with full graph
        size_t d = num_edges(prev_tree);
        VineTree new_tree(d);

        // cut & paste information from previous tree
        int i = 0;
        for (auto e : boost::edges(prev_tree)) {
            new_tree[i].hfunc1 = prev_tree[e].hfunc1;
            new_tree[i].hfunc2 = prev_tree[e].hfunc2;
            new_tree[i].conditioning = prev_tree[e].conditioning;
            new_tree[i].conditioned = prev_tree[e].conditioned;
            new_tree[i].prev_edge_indices.reserve(2);
            new_tree[i].prev_edge_indices.push_back(boost::source(e, prev_tree));
            new_tree[i].prev_edge_indices.push_back(boost::target(e, prev_tree));
            ++i;
        }

        return new_tree;
    }

    //! Add edges allowed by the proximity condition
    //!
    //! Also calculates the edge weight (e.g., 1-|tau| for tree_criterion =
    //! "itau").
    //!
    //! @param vine_tree tree of a vine.
    //! @param tree_criterion the criterion for selecting the maximum spanning
    //!     tree ("tau", "hoeffd" and "rho" implemented so far).
    //! @param threshold for thresholded vines.
    void add_allowed_edges(VineTree& vine_tree, std::string tree_criterion,
                           double threshold)
    {
        for (auto v0 : boost::vertices(vine_tree)) {
            for (size_t v1 = 0; v1 < v0; ++v1) {
                // check proximity condition: common neighbor in previous tree
                // (-1 means 'no common neighbor')
                if (find_common_neighbor(v0, v1, vine_tree) > -1) {
                    auto pc_data = get_pc_data(v0, v1, vine_tree);
                    auto w = get_tree_criterion(pc_data, tree_criterion, threshold);
                    auto e = boost::add_edge(v0, v1, w, vine_tree).first;
                    vine_tree[e].weight = w;
                }
            }
        }
    }

    double get_tree_criterion(Eigen::Matrix<double, Eigen::Dynamic, 2> data,
                              std::string tree_criterion, double threshold)
    {
        double w;
        if (tree_criterion == "tau") {
            w = 1.0 - std::fabs(tools_stats::pairwise_ktau(data));
        } else if (tree_criterion == "hoeffd") {
            // scale to [0,1]
            w = 1.0 - (30*tools_stats::pairwise_hoeffd(data)+0.5)/1.5;
        } else if (tree_criterion == "rho") {
            w = 1.0 - std::fabs(tools_stats::pairwise_cor(data));
        } else {
            throw std::runtime_error("tree criterion not implemented");
        }
        if (w > 1 - threshold) {
            w = 1.0;
        }

        return w;
    }

    // Find common neighbor in previous tree
    //
    // @param v0,v1 vertices in the tree.
    // @param tree the current tree.
    // @return Gives the index of the vertex in the previous tree that was
    //     shared by e0, e1, the edge representations of v0, v1.
    ptrdiff_t find_common_neighbor(size_t v0, size_t v1, const VineTree& tree)
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

    // Extract pair copula pseudo-observations from h-functions
    //
    // @param v0,v1 vertex indices.
    // @param tree a vine tree.
    // @return The pseudo-observations for the pair coula, extracted from
    //     the h-functions calculated in the previous tree.
    Eigen::MatrixXd get_pc_data(size_t v0, size_t v1, const VineTree& tree)
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

    //! Collapse a graph to the minimum spanning tree
    //!
    //! @param graph the input graph.
    //! @return the input graph with all non-MST edges removed.
    void min_spanning_tree(VineTree &graph)
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

    //! Add conditioning info and data for each edge
    //!
    //! See, e.g., Czado (2010), "Pair-copula constructions of multivariate
    //! copulas", url: https://mediatum.ub.tum.de/doc/1079253/file.pdf
    //! @param tree a vine tree.
    void add_edge_info(VineTree& tree)
    {
        for (auto e : boost::edges(tree)) {
            auto v0 = boost::source(e, tree);
            auto v1 = boost::target(e, tree);
            tree[e].pc_data = get_pc_data(v0, v1, tree);

            auto v0_indices = cat(tree[v0].conditioning, tree[v0].conditioned);
            auto v1_indices = cat(tree[v1].conditioning, tree[v1].conditioned);

            auto test = intersect(v0_indices, v1_indices);
            auto d01 = set_diff(v0_indices, v1_indices);
            auto d10 = set_diff(v1_indices, v0_indices);

            tree[e].conditioning = cat(d01, d10);
            tree[e].conditioned = intersect(v0_indices, v1_indices);
            tree[e].all_indices = cat(tree[e].conditioning, tree[e].conditioned);
        }
    }

    //! Remove data (hfunc1/hfunc2/pc_data) from all edges of a vine tree
    //! @param tree a vine tree.
    void remove_edge_data(VineTree& tree)
    {
        for (auto e : boost::edges(tree)) {
            tree[e].hfunc1 = Eigen::VectorXd();
            tree[e].hfunc2 = Eigen::VectorXd();
            tree[e].pc_data = Eigen::Matrix<double, Eigen::Dynamic, 2>(0, 2);
        }
    }

    //! Remove data (hfunc1/hfunc2/pc_data) from all vertices of a vine tree
    //! @param tree a vine tree.
    void remove_vertex_data(VineTree& tree)
    {
        for (auto v : boost::vertices(tree)) {
            tree[v].hfunc1 = Eigen::VectorXd();
            tree[v].hfunc2 = Eigen::VectorXd();
        }
    }

    //! Fit and select a pair copula for each edges
    //! @param tree a vine tree preprocessed with add_edge_info().
    //! @param controls the controls for fitting a vine copula (see FitControlsVinecop)
    void select_pair_copulas(VineTree& tree,
                             vinecopulib::FitControlsVinecop& controls)
    {
        // Controls for selecting the pair copulas
        vinecopulib::FitControlsBicop controls_bicop(controls.get_fit_controls_bicop());

        for (auto e : boost::edges(tree)) {
            if (tree[e].weight > 1 - controls.get_threshold()) {
                tree[e].pair_copula = vinecopulib::Bicop();
            } else {
                tree[e].pair_copula = vinecopulib::Bicop(tree[e].pc_data,
                                                         controls_bicop);
            }

            tree[e].hfunc1 = tree[e].pair_copula.hfunc1(tree[e].pc_data);
            tree[e].hfunc2 = tree[e].pair_copula.hfunc2(tree[e].pc_data);
        }
    }

    //! Print indices, family, and parameters for each pair-copula
    //! @param tree a vine tree.
    void print_pair_copulas(VineTree& tree)
    {
        for (auto e : boost::edges(tree)) {
            std::stringstream pc_info;
            pc_info << get_pc_index(e, tree) << " <-> " <<
                        tree[e].pair_copula.str() << std::endl;
            vinecopulib::tools_interface::print(pc_info.str().c_str());
        }
    }

    //! Get edge index for the vine (like 1, 2; 3)
    //! @param e a descriptor for the edge.
    //! @param tree a vine tree.
    std::string get_pc_index(
            boost::graph_traits<VineTree>::edge_descriptor e,
            VineTree& tree
    )
    {
        std::stringstream index;
        // add 1 everywhere for user-facing representation (boost::graph starts
        // at 0)
        index <<
              tree[e].conditioning[0] + 1 <<
              "," <<
              tree[e].conditioning[1] + 1;
        if (tree[e].conditioned.size() > 0) {
            index << " ; ";
            for (unsigned int i = 0; i < tree[e].conditioned.size(); ++i) {
                index << tree[e].conditioned[i] + 1;
                if (i < tree[e].conditioned.size() - 1)
                    index << ",";
            }
        }

        return index.str().c_str();
    }

}

}
