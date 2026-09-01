// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

namespace vinecopulib {

//! @brief Builds the tree list from an `(order, struct_array)` representation.
//! @param order The variable order (diagonal of the R-vine matrix).
//! @param struct_array The structure array, labelled consistently with `order`.
inline RVineTrees::RVineTrees(const std::vector<size_t>& order,
                              const TriangularArray<size_t>& struct_array)
  : RVineTrees(order, struct_array, {})
{
}

//! @brief Builds the tree list, attaching pair-copulas to each edge.
//! @param order The variable order (diagonal of the R-vine matrix).
//! @param struct_array The structure array, labelled consistently with `order`.
//! @param pair_copulas The pair-copulas, indexed `[tree][edge]`; each is stored
//!   with its first argument aligned to the diagonal variable `order[edge]`. If
//!   empty, every edge is independence; otherwise it must cover every tree of
//!   `struct_array`.
inline RVineTrees::RVineTrees(
  const std::vector<size_t>& order,
  const TriangularArray<size_t>& struct_array,
  const std::vector<std::vector<Bicop>>& pair_copulas)
{
  if (order.size() != struct_array.get_dim()) {
    throw std::runtime_error(
      "Order and structure array dimensions do not match.");
  }
  d_ = order.size();
  trunc_lvl_ = struct_array.get_trunc_lvl();
  const bool with_copulas = !pair_copulas.empty();
  if (with_copulas && (pair_copulas.size() < trunc_lvl_)) {
    throw std::runtime_error(
      "pair_copulas covers " + std::to_string(pair_copulas.size()) +
      " trees, the structure array " + std::to_string(trunc_lvl_) +
      "; pass an empty store to mean independence throughout.");
  }
  trees_.resize(trunc_lvl_);
  for (size_t t = 0; t < trunc_lvl_; ++t) {
    if (with_copulas && (pair_copulas[t].size() != d_ - 1 - t)) {
      throw std::runtime_error("pair_copulas[" + std::to_string(t) + "] has " +
                               std::to_string(pair_copulas[t].size()) +
                               " entries, expected " +
                               std::to_string(d_ - 1 - t) + ".");
    }
    for (size_t e = 0; e < d_ - 1 - t; ++e) {
      std::vector<size_t> conditioning;
      conditioning.reserve(t);
      for (size_t k = 0; k < t; ++k)
        conditioning.push_back(struct_array(k, e));
      trees_[t].push_back(
        with_copulas
          ? Edge(order[e], struct_array(t, e), conditioning, pair_copulas[t][e])
          : Edge(order[e], struct_array(t, e), conditioning));
    }
  }
}

//! @brief Builds directly from a list of trees (edges carry their pair-copula).
//! @param d The dimension of the vine.
//! @param trees One edge list per tree; `trees[t]` must have `d - 1 - t` edges
//!   for a valid vine (validated on conversion, not here).
inline RVineTrees::RVineTrees(size_t d, std::vector<Tree> trees)
  : d_(d)
  , trunc_lvl_(trees.size())
  , trees_(std::move(trees))
{
}

//! @brief The default diagonal policy: the `conditioned[0]` (first) leaf
//! endpoint of the first leaf edge — flip-free, and shared by `select()`
//! finalization and the `RVineStructure` round-trip.
inline RVineTrees::DiagonalPolicy
RVineTrees::default_diagonal_policy()
{
  return
    [](size_t, const std::vector<std::vector<size_t>>& leaf_edges) noexcept {
      return leaf_edges[0][0];
    };
}

//! @brief Converts back to matrix form, choosing diagonals with a policy and
//! carrying (flipping) the pair-copulas.
//! @param diagonal_policy Chooses which leaf variable sits on each column's
//!   diagonal (see `DiagonalPolicy`).
inline RVineTrees::Decomposition
RVineTrees::to_struct_array(const DiagonalPolicy& diagonal_policy) const
{
  check_tree_sizes();
  return peel(diagonal_policy, true);
}

inline RVineTrees::Decomposition
RVineTrees::to_struct_array_map(const DiagonalPolicy& diagonal_policy) const
{
  check_tree_sizes();
  return peel(diagonal_policy, false);
}

inline void
RVineTrees::check_tree_sizes() const
{
  for (size_t t = 0; t < trunc_lvl_; ++t) {
    if (trees_[t].size() != d_ - 1 - t) {
      throw std::runtime_error(
        "Tree " + std::to_string(t) +
        " does not have the correct number of edges. Expected " +
        std::to_string(d_ - 1 - t) + " edges, but got " +
        std::to_string(trees_[t].size()) + ".");
    }
  }
}

//! @brief Shared leaf-peeling: fills the R-vine matrix column by column.
inline RVineTrees::Decomposition
RVineTrees::peel(const DiagonalPolicy& diagonal_policy,
                 bool carry_copulas) const
{
  std::vector<AugmentedTree> augmented = trees_to_augmented();

  std::vector<size_t> order(d_);
  TriangularArray<size_t> struct_array(d_, trunc_lvl_);
  TriangularArray<PairCopulaLocation> pair_copula_locations(d_, trunc_lvl_);
  std::vector<std::vector<Bicop>> pair_copulas;
  if (carry_copulas) {
    pair_copulas.resize(trunc_lvl_);
    for (size_t t = 0; t < trunc_lvl_; ++t)
      pair_copulas[t].resize(d_ - 1 - t);
  }

  for (size_t col = 0; col < d_ - 1; ++col) {
    size_t t =
      std::max(std::min(trunc_lvl_, d_ - 1 - col), static_cast<size_t>(1));
    AugmentedTree& tree = augmented[t - 1];

    // gather the leaf edges of the top tree in iteration order; each yields
    // its leaf endpoint(s) as diagonal options (endpoint `a` before `b`)
    std::vector<std::vector<size_t>> leaf_edges;
    for (const auto& e : tree.edges) {
      if (e.consumed)
        continue;
      const Edge& ed = *e.edge;
      std::vector<size_t> options;
      if (tree.degrees[e.node1] == 1)
        options.push_back(ed.a);
      if (tree.degrees[e.node2] == 1)
        options.push_back(ed.b);
      if (!options.empty())
        leaf_edges.push_back(std::move(options));
    }
    if (leaf_edges.empty()) {
      throw std::runtime_error("RVineTrees: no leaf found while peeling tree " +
                               std::to_string(t - 1) + ".");
    }
    size_t diag = diagonal_policy(col, leaf_edges);
    order[col] = diag;

    // locate the (unique) leaf edge carrying `diag`, consume it, and seed the
    // running conditioning set for the descent
    std::vector<size_t> check_set;
    bool found = false;
    for (auto& e : tree.edges) {
      if (e.consumed)
        continue;
      const Edge& ed = *e.edge;
      bool as_a = (tree.degrees[e.node1] == 1) && (ed.a == diag);
      bool as_b = (tree.degrees[e.node2] == 1) && (ed.b == diag);
      if (!as_a && !as_b)
        continue;
      struct_array(t - 1, col) = as_a ? ed.b : ed.a;
      if (carry_copulas) {
        pair_copulas[t - 1][col] = ed.pair_copula;
        if (diag != ed.a)
          pair_copulas[t - 1][col].flip();
      }
      pair_copula_locations(t - 1,
                            col) = { t - 1, e.source_edge, diag != ed.a };
      check_set = ed.C;
      tree.degrees[e.node1]--;
      tree.degrees[e.node2]--;
      e.consumed = true;
      found = true;
      break;
    }
    if (!found) {
      throw std::runtime_error(
        "RVineTrees: diagonal variable is not an available leaf.");
    }

    // descend the column: at each lower tree, match the edge whose variables
    // are exactly `{diag} ∪ check_set`, carrying/flipping its pair-copula
    for (size_t k = 1; k < t; ++k) {
      size_t tree_idx = t - 1 - k;
      AugmentedTree& ltree = augmented[tree_idx];
      tools_stl::insert_sorted(check_set, diag);
      bool matched = false;
      for (auto& e : ltree.edges) {
        if (e.consumed)
          continue;
        const Edge& ed = *e.edge;
        if (ed.all_indices != check_set)
          continue;
        struct_array(tree_idx, col) = (ed.a == diag) ? ed.b : ed.a;
        if (carry_copulas) {
          pair_copulas[tree_idx][col] = ed.pair_copula;
          if (diag != ed.a)
            pair_copulas[tree_idx][col].flip();
        }
        pair_copula_locations(tree_idx,
                              col) = { tree_idx, e.source_edge, diag != ed.a };
        check_set = ed.C;
        ltree.degrees[e.node1]--;
        ltree.degrees[e.node2]--;
        e.consumed = true;
        matched = true;
        break;
      }
      if (!matched) {
        throw std::runtime_error(
          "RVineTrees: proximity condition violated while peeling.");
      }
    }
  }

  order[d_ - 1] = struct_array(0, d_ - 2);
  return Decomposition{ std::move(order),
                        std::move(struct_array),
                        std::move(pair_copulas),
                        std::move(pair_copula_locations) };
}

//! @brief Builds the map from `(variable, conditioning ∪ partner)` to the edge
//! index in the previous tree (i.e. the incident node in the line graph).
inline std::map<std::pair<size_t, std::vector<size_t>>, size_t>
RVineTrees::build_lookup(const std::vector<AugmentedEdge>& edges) const
{
  std::map<std::pair<size_t, std::vector<size_t>>, size_t> lookup;
  for (size_t i = 0; i < edges.size(); ++i) {
    const Edge& ed = *edges[i].edge;
    std::vector<size_t> key = ed.C; // already sorted ascending
    tools_stl::insert_sorted(key, ed.b);
    lookup[{ ed.a, key }] = i;
    key = ed.C;
    tools_stl::insert_sorted(key, ed.a);
    lookup[{ ed.b, key }] = i;
  }
  return lookup;
}

//! @brief Checks that tree 0 mentions every variable `1, ..., d`.
inline void
RVineTrees::check_missing_vars(const std::vector<RVineTrees::Edge>& edges,
                               size_t d) const
{
  std::vector<bool> seen(d + 1, false); // index 0 unused (1-based variables)
  for (const auto& ed : edges) {
    if (ed.a >= 1 && ed.a <= d)
      seen[ed.a] = true;
    if (ed.b >= 1 && ed.b <= d)
      seen[ed.b] = true;
  }
  std::vector<size_t> missing;
  for (size_t i = 1; i <= d; ++i)
    if (!seen[i])
      missing.push_back(i);
  if (!missing.empty()) {
    std::string problem =
      "Tree 0 is missing " + std::to_string(missing.size()) + " variables: ";
    for (auto v : missing)
      problem += std::to_string(v) + " ";
    throw std::runtime_error(problem);
  }
}

//! @brief Builds the line-graph (augmented) view of every tree, validating the
//! proximity condition against the tree below.
inline std::vector<RVineTrees::AugmentedTree>
RVineTrees::trees_to_augmented() const
{
  std::vector<AugmentedTree> augmented(trunc_lvl_);
  for (size_t t = 0; t < trunc_lvl_; ++t) {
    const Tree& edges = trees_[t];
    AugmentedTree atree;
    atree.edges.reserve(edges.size());
    if (t == 0) {
      check_missing_vars(edges, d_);
      for (size_t i = 0; i < edges.size(); ++i) {
        const auto& ed = edges[i];
        atree.edges.push_back({ ed.a, ed.b, i, &ed, false });
      }
      atree.degrees.assign(d_ + 1, 0); // 1-based labels; index 0 unused
    } else {
      auto lookup = build_lookup(augmented[t - 1].edges);
      for (size_t i = 0; i < edges.size(); ++i) {
        const Edge& ed = edges[i];
        auto it1 = lookup.find({ ed.a, ed.C });
        auto it2 = lookup.find({ ed.b, ed.C });
        if (it1 == lookup.end() || it2 == lookup.end()) {
          std::string problem =
            "Proximity condition violated in tree " + std::to_string(t) +
            " for edge " + std::to_string(i) + " (" + std::to_string(ed.a) +
            ", " + std::to_string(ed.b) + "; ";
          for (auto v : ed.C)
            problem += std::to_string(v) + " ";
          problem += ")";
          throw std::runtime_error(problem);
        }
        atree.edges.push_back({ it1->second, it2->second, i, &ed, false });
      }
      atree.degrees.assign(augmented[t - 1].edges.size(), 0);
    }
    for (const auto& e : atree.edges) {
      atree.degrees[e.node1]++;
      atree.degrees[e.node2]++;
    }
    augmented[t] = std::move(atree);
  }
  return augmented;
}

}
