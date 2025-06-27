// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <vinecopulib/misc/tools_stl.hpp>
#include <vinecopulib/misc/triangular_array.hpp>

namespace vinecopulib {

class RVineStructure; // forward declaration

/**
 * @brief Represents a list-of-trees decomposition of an R-vine structure.
 * Each tree consists of edges defined_ by conditioned_ pairs and_ conditioning
 * sets.
 */
class RVineTrees
{
public:
  using Edge = std::tuple<size_t, size_t, std::set<size_t>>;
  using Tree = std::vector<Edge>;

  struct AugmentedTree
  {
    // std::map<size_t, Edge> nodes; // Not used in this implementation
    std::vector<std::tuple<size_t, size_t, Edge>> edges;
    std::map<size_t, size_t> degrees;
  };

  // default constructor allows empty RVineTrees
  RVineTrees()
    : d_(0)
    , trunc_lvl_(0)
    , trees_() {};

  explicit RVineTrees(const std::vector<size_t>& order,
                      const TriangularArray<size_t>& struct_array);
  explicit RVineTrees(const std::vector<RVineTrees>& vines);

  const std::vector<Tree>& get_trees() const { return trees_; }

  std::tuple<std::vector<size_t>, TriangularArray<size_t>> to_struct_array(
    bool fill_missing = false) const;

  size_t get_dim() const { return d_; }
  size_t get_trunc_lvl() const { return trunc_lvl_; }

private:
  size_t d_;
  size_t trunc_lvl_;
  std::vector<Tree> trees_;

  std::map<std::pair<size_t, std::set<size_t>>, size_t> build_lookup(
    const std::vector<std::tuple<size_t, size_t, Edge>>& edges) const;

  void check_missing_vars(const std::vector<Edge>& edges, size_t d) const;
  std::vector<AugmentedTree> trees_to_augmented() const;
  void fill_missing_edges(
    std::vector<AugmentedTree>& augmented) const;
};

}

#include <vinecopulib/vinecop/implementation/rvine_trees.ipp>
