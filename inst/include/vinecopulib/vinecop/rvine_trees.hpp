// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <algorithm>
#include <functional>
#include <map>
#include <tuple>
#include <vector>
#include <vinecopulib/bicop/class.hpp>
#include <vinecopulib/misc/tools_stl.hpp>
#include <vinecopulib/misc/triangular_array.hpp>

namespace vinecopulib {

//! @brief A list-of-trees decomposition of an R-vine structure.
//!
//! @details The vine is stored as one edge list per tree. Each edge holds its
//! conditioned pair `(a, b)`, its conditioning set `C`, and (optionally) the
//! fitted pair-copula. The class converts to and from the `(order,
//! struct_array)` representation used by `RVineStructure`, validating the
//! proximity condition on the way back. It is the shared primitive behind
//! `Vinecop::reorient()` and the structure finalisation during selection.
//!
//! The orientation convention is that an edge's pair-copula has its first
//! argument aligned with the conditioned variable `a`; a copula is therefore
//! flipped iff the variable placed on the matrix diagonal differs from `a`.
class RVineTrees
{
public:
  //! @brief A single edge: conditioned pair `(a, b)`, conditioning set `C`,
  //! and the associated pair-copula (independence by default).
  struct Edge
  {
    size_t a{ 0 };           //!< first conditioned variable (arg1)
    size_t b{ 0 };           //!< second conditioned variable
    std::vector<size_t> C{}; //!< conditioning set, kept sorted ascending
    std::vector<size_t> all_indices{}; //!< sorted {a} ∪ {b} ∪ C (precomputed)
    Bicop pair_copula{}; //!< pair-copula (independence by default)

    Edge() = default;
    Edge(size_t a_arg, size_t b_arg, std::vector<size_t> c_arg)
      : a(a_arg)
      , b(b_arg)
      , C(std::move(c_arg))
    {
      finalize();
    }
    Edge(size_t a_arg, size_t b_arg, std::vector<size_t> c_arg, Bicop pc)
      : a(a_arg)
      , b(b_arg)
      , C(std::move(c_arg))
      , pair_copula(std::move(pc))
    {
      finalize();
    }

  private:
    //! Canonicalizes `C` (sort + unique) and precomputes `all_indices`.
    void finalize()
    {
      std::sort(C.begin(), C.end());
      C.erase(std::unique(C.begin(), C.end()), C.end());
      all_indices = C;
      all_indices.reserve(C.size() + 2);
      tools_stl::insert_sorted(all_indices, a);
      tools_stl::insert_sorted(all_indices, b);
    }
  };
  using Tree = std::vector<Edge>;

  //! @brief Chooses which variable sits on the diagonal of a given column.
  //!
  //! @details Called with the column index and, for each leaf edge of the
  //! current top tree (in iteration order), the variables that may go on the
  //! diagonal (its leaf endpoints, 1 or 2 of them). Must return one of them.
  using DiagonalPolicy =
    std::function<size_t(size_t col,
                         const std::vector<std::vector<size_t>>& leaf_edges)>;

  //! @brief The result of a copula-aware `to_struct_array()`.
  struct Decomposition
  {
    std::vector<size_t> order;
    TriangularArray<size_t> struct_array;
    std::vector<std::vector<Bicop>> pair_copulas; //!< indexed `[tree][edge]`
  };

  //! @brief One edge of the augmented (line-graph) view: its two incident node
  //! ids plus a non-owning pointer to the underlying `Edge` (in `trees_`).
  struct AugmentedEdge
  {
    size_t node1{ 0 };
    size_t node2{ 0 };
    const Edge* edge{ nullptr };
    bool consumed{ false };
  };

  //! @brief The augmented (line-graph) view of a tree used during conversion.
  //! `degrees` is indexed by node id (contiguous), avoiding a `std::map`.
  struct AugmentedTree
  {
    std::vector<AugmentedEdge> edges;
    std::vector<size_t> degrees;
  };

  RVineTrees()
    : d_(0)
    , trunc_lvl_(0)
    , trees_()
  {
  }

  explicit RVineTrees(const std::vector<size_t>& order,
                      const TriangularArray<size_t>& struct_array);
  RVineTrees(const std::vector<size_t>& order,
             const TriangularArray<size_t>& struct_array,
             const std::vector<std::vector<Bicop>>& pair_copulas);
  RVineTrees(size_t d, std::vector<Tree> trees);

  const std::vector<Tree>& get_trees() const { return trees_; }

  //! @brief The default diagonal policy: the `conditioned[0]` (first) leaf
  //! endpoint of the first leaf edge. It is flip-free — each edge stores its
  //! pair copula aligned to `conditioned[0]`, so placing that endpoint on the
  //! diagonal needs no flip — and it is the convention shared by `select()`
  //! finalization and the `RVineStructure` round-trip.
  static DiagonalPolicy default_diagonal_policy();

  //! @brief Converts back to matrix form, choosing diagonals with
  //! `diagonal_policy` (default: `default_diagonal_policy()`) and
  //! carrying/flipping the pair-copulas.
  Decomposition to_struct_array(
    const DiagonalPolicy& diagonal_policy = default_diagonal_policy()) const;

  size_t get_dim() const { return d_; }
  size_t get_trunc_lvl() const { return trunc_lvl_; }

private:
  size_t d_;
  size_t trunc_lvl_;
  std::vector<Tree> trees_;

  std::map<std::pair<size_t, std::vector<size_t>>, size_t> build_lookup(
    const std::vector<AugmentedEdge>& edges) const;
  void check_missing_vars(const std::vector<Edge>& edges, size_t d) const;
  std::vector<AugmentedTree> trees_to_augmented() const;
  Decomposition peel(const DiagonalPolicy& diagonal_policy,
                     bool carry_copulas) const;
};

}

#include <vinecopulib/vinecop/implementation/rvine_trees.ipp>
