// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

namespace vinecopulib {

inline RVineTrees::RVineTrees(
  const std::vector<size_t>& order,
  const TriangularArray<size_t>& struct_array) // TODO: natural order
{
  if (order.size() != struct_array.get_dim()) {
    throw std::runtime_error(
      "Order and structure array dimensions do not match.");
  }
  d_ = order.size();
  trunc_lvl_ = struct_array.get_trunc_lvl();
  trees_.clear();
  trees_.resize(trunc_lvl_);

  for (size_t t = 0; t < trunc_lvl_; ++t) {
    Tree tree;
    for (size_t e = 0; e < d_ - 1 - t; ++e) {
      size_t a = order[e];
      size_t b = static_cast<size_t>(struct_array(t, e));
      std::set<size_t> C;
      for (size_t k = 0; k < t; ++k) {
        C.insert(struct_array(k, e));
      }
      tree.push_back(Edge(a, b, C));
    }
    trees_[t] = tree;
  }
}

inline std::tuple<std::vector<size_t>, TriangularArray<size_t>>
RVineTrees::to_struct_array(bool fill_missing) const
{
  if (!fill_missing) {
    for (size_t t = 0; t < trunc_lvl_; ++t) {
      if (trees_[t].size() != d_ - 1 - t) {
        throw std::runtime_error(
          "Tree " + std::to_string(t) +
          " does not have the correct number of edges. " + "Expected " +
          std::to_string(d_ - 1 - t) + " edges, but got " +
          std::to_string(trees_[t].size()) + ".");
      }
    }
  }

  std::vector<AugmentedTree> augmented = trees_to_augmented();
  if (fill_missing) {
    fill_missing_edges(augmented);
  }

  TriangularArray<size_t> struct_array(d_, trunc_lvl_);
  std::vector<size_t> order(d_);
  for (size_t col = 0; col < d_ - 1; ++col) {
    size_t t =
      std::max(std::min(trunc_lvl_, d_ - 1 - col), static_cast<size_t>(1));
    AugmentedTree& tree = augmented[t - 1];

    Edge edge;
    size_t a, b, diag, offdiag, node1, node2;
    std::set<size_t> check_set, C, all_indices;
    for (typename std::vector<std::tuple<size_t, size_t, Edge>>::iterator it =
           tree.edges.begin();
         it != tree.edges.end();
         ++it) {
      node1 = std::get<0>(*it);
      node2 = std::get<1>(*it);

      if (tree.degrees[node1] == 1 || tree.degrees[node2] == 1) {
        edge = std::get<2>(*it);
        a = std::get<0>(edge);
        b = std::get<1>(edge);
        check_set = std::get<2>(edge);
        diag = (tree.degrees[node1] == 1) ? a : b;
        offdiag = (diag == a) ? b : a;
        struct_array(t - 1, col) = offdiag;
        order[col] = diag;
        tree.degrees[node1]--;
        tree.degrees[node2]--;
        tree.edges.erase(it);
        break;
      }
    }

    for (size_t k = 1; k < t; ++k) {
      size_t tree_idx = t - 1 - k;
      AugmentedTree& ltree = augmented[tree_idx];
      check_set.insert(order[col]);

      for (typename std::vector<std::tuple<size_t, size_t, Edge>>::iterator it =
             ltree.edges.begin();
           it != ltree.edges.end();
           ++it) {
        node1 = std::get<0>(*it);
        node2 = std::get<1>(*it);
        edge = std::get<2>(*it);
        a = std::get<0>(edge);
        b = std::get<1>(edge);
        C = std::get<2>(edge);

        all_indices = C;
        all_indices.insert(a);
        all_indices.insert(b);

        if (all_indices == check_set) {
          size_t next_var = (a == order[col]) ? b : a;
          check_set = C;
          struct_array(tree_idx, col) = next_var;
          ltree.degrees[node1]--;
          ltree.degrees[node2]--;
          ltree.edges.erase(it);
          break;
        }
      }
    }
  }

  // The last column contains a single element which must be different
  // from all other diagonal elements. Based on the properties of an
  // R-vine matrix, this must be the element next to it.
  order[d_ - 1] = struct_array(0, d_ - 2);

  return std::make_tuple(order, struct_array);
}

inline std::map<std::pair<size_t, std::set<size_t>>, size_t>
RVineTrees::build_lookup(
  const std::vector<std::tuple<size_t, size_t, Edge>>& edges) const
{
  size_t idx = 0;
  std::map<std::pair<size_t, std::set<size_t>>, size_t> lookup;
  for (std::size_t i = 0; i < edges.size(); ++i) {
    const Edge& edge = std::get<2>(edges[i]);
    size_t a = std::get<0>(edge);
    size_t b = std::get<1>(edge);
    std::set<size_t> C = std::get<2>(edge);

    std::set<size_t> Ca = C;
    Ca.insert(b);
    lookup[{ a, Ca }] = idx;
    std::set<size_t> Cb = C;
    Cb.insert(a);
    lookup[{ b, Cb }] = idx;
    ++idx;
  }
  return lookup;
}

inline void
RVineTrees::check_missing_vars(const std::vector<RVineTrees::Edge>& edges,
                               size_t d) const
{
  std::vector<bool> seen(d + 1,
                         false); // index 0 is unused for 1-based indexing

  for (const auto& edge : edges) {
    size_t a = std::get<0>(edge);
    size_t b = std::get<1>(edge);
    if (a >= 1 && a <= d)
      seen[a] = true;
    if (b >= 1 && b <= d)
      seen[b] = true;
  }

  std::vector<size_t> missing;
  for (size_t i = 1; i <= d; ++i) {
    if (!seen[i])
      missing.push_back(i);
  }
  if (!missing.empty()) {
    std::string problem =
      "Tree 0 is missing " + std::to_string(missing.size()) + " variables: ";
    for (auto v : missing) {
      problem += std::to_string(v) + " ";
    }
    throw std::runtime_error(problem);
  }
}

inline std::vector<RVineTrees::AugmentedTree>
RVineTrees::trees_to_augmented() const
{
  std::vector<AugmentedTree> augmented;
  augmented.resize(trunc_lvl_);
  for (std::size_t t = 0; t < trunc_lvl_; ++t) {
    const Tree& edges = trees_[t];
    AugmentedTree atree;

    if (t == 0) {
      check_missing_vars(edges, d_);

      for (std::size_t i = 0; i < edges.size(); ++i) {
        size_t a = std::get<0>(edges[i]);
        size_t b = std::get<1>(edges[i]);
        std::set<size_t> C = std::get<2>(edges[i]);
        atree.edges.push_back(std::make_tuple(a, b, Edge(a, b, C)));
      }
    } else {
      auto nodes = augmented[t - 1].edges;
      auto lookup = build_lookup(nodes);

      for (std::size_t i = 0; i < edges.size(); ++i) {
        size_t a = std::get<0>(edges[i]);
        size_t b = std::get<1>(edges[i]);
        std::set<size_t> C = std::get<2>(edges[i]);

        auto it1 = lookup.find({ a, C });
        auto it2 = lookup.find({ b, C });
        if (it1 == lookup.end() || it2 == lookup.end()) {
          std::string problem = "Proximity condition violated in tree " +
                                std::to_string(t) + " for edge " +
                                std::to_string(i) + " (" + std::to_string(a) +
                                ", " + std::to_string(b) + "; ";
          for (auto v : C) {
            problem += std::to_string(v) + " ";
          }
          problem += ")";
          throw std::runtime_error(problem);
        }
        size_t node1 = it1->second;
        size_t node2 = it2->second;
        atree.edges.push_back(std::make_tuple(node1, node2, Edge(a, b, C)));
      }
    }

    for (std::size_t i = 0; i < atree.edges.size(); ++i) {
      size_t node1 = std::get<0>(atree.edges[i]);
      size_t node2 = std::get<1>(atree.edges[i]);
      atree.degrees[node1]++;
      atree.degrees[node2]++;
    }
    augmented[t] = atree;
  }
  return augmented;
}

inline RVineTrees::RVineTrees(const std::vector<RVineTrees>& vines)
{
  if (vines.empty()) {
    throw std::runtime_error("Cannot merge empty list of RVineTrees.");
  }

  // Infer d_ from unique variable indices in all Tree 0 edges
  std::set<size_t> unique_vars;
  trunc_lvl_ = 0;
  for (const auto& vine : vines) {
    if (vine.trees_.empty()) {
      throw std::runtime_error(
        "Each RVineTrees must contain at least one tree.");
    }
    trunc_lvl_ = std::max(trunc_lvl_, vine.get_trunc_lvl());
    for (const auto& edge : vine.trees_[0]) {
      unique_vars.insert(std::get<0>(edge));
      unique_vars.insert(std::get<1>(edge));
    }
  }
  d_ = unique_vars.size();

  for (size_t i = 1; i <= d_; ++i) {
    if (!unique_vars.count(i)) {
      throw std::runtime_error("First trees contain " + std::to_string(d_) +
                               "variables, but " + std::to_string(i) +
                               " is missing.");
    }
  }

  trees_.clear();
  trees_.resize(trunc_lvl_);

  for (size_t t = 0; t < trunc_lvl_; ++t) {
    for (const auto& vine : vines) {
      if (t < vine.get_trunc_lvl()) {
        trees_[t].insert(
          trees_[t].end(), vine.trees_[t].begin(), vine.trees_[t].end());
      }
    }
  }
}

inline void
RVineTrees::fill_missing_edges(
  std::vector<RVineTrees::AugmentedTree>& augmented) const
{
  const size_t d = d_;

  for (size_t t = 0; t < augmented.size(); ++t) {
    auto& tree = augmented[t];
    const size_t target_edges = d - 1 - t;

    if (tree.edges.size() == target_edges) {
      continue; // already full
    }

    // Reconstruct nodes
    std::map<size_t, Edge> nodes;
    std::map<size_t, std::set<size_t>> node_to_vars;

    if (t == 0) {
      // Collect unique variables
      std::set<size_t> var_set;
      for (const auto& edge : tree.edges) {
        var_set.insert(std::get<0>(edge));
        var_set.insert(std::get<1>(edge));
      }

      std::vector<size_t> var_list(var_set.begin(), var_set.end());
      std::map<size_t, size_t> var_to_node;
      for (size_t idx = 0; idx < var_list.size(); ++idx) {
        nodes[idx] = Edge(var_list[idx], 0, {}); // info is variable only
        node_to_vars[idx] = { var_list[idx] };
        var_to_node[var_list[idx]] = idx;
      }

      // Remap edges using node indices
      std::vector<std::tuple<size_t, size_t, Edge>> new_edges;
      for (const auto& edge : tree.edges) {
        size_t var1 = std::get<0>(edge);
        size_t var2 = std::get<1>(edge);
        size_t node1 = var_to_node.at(var1);
        size_t node2 = var_to_node.at(var2);
        new_edges.emplace_back(node1, node2, Edge(var1, var2, {}));
      }
      tree.edges = std::move(new_edges);
    } else {
      const auto& prev_edges = augmented[t - 1].edges;
      for (size_t idx = 0; idx < prev_edges.size(); ++idx) {
        const Edge& e = std::get<2>(prev_edges[idx]);
        nodes[idx] = e;
        std::set<size_t> vars = { std::get<0>(e), std::get<1>(e) };
        vars.insert(std::get<2>(e).begin(), std::get<2>(e).end());
        node_to_vars[idx] = std::move(vars);
      }
    }

    // Union-find structure to detect cycles
    tools_stl::UnionFind<size_t> uf;
    // Set to keep track of already connected nodes
    std::set<std::pair<size_t, size_t>> connected;

    // Add nodes to union-find
    std::vector<size_t> node_list;
    for (const auto& kv : nodes) {
      size_t node = kv.first;
      node_list.push_back(node);
      uf.add(node);
    }

    // Add edges to union-find and connected set
    for (const auto& e : tree.edges) {
      size_t n1 = std::get<0>(e), n2 = std::get<1>(e);
      if (!uf.unite(n1, n2)) {
        throw std::runtime_error("Tree " + std::to_string(t) +
                                 " contains a cycle.");
      }
      connected.insert({ n1, n2 });
      connected.insert({ n2, n1 });
    }

    // Propose candidate edges
    std::vector<std::tuple<size_t, size_t, Edge>> candidates;
    for (size_t i = 0; i < node_list.size(); ++i) {
      for (size_t j = i + 1; j < node_list.size(); ++j) {
        size_t n1 = node_list[i], n2 = node_list[j];
        if (connected.count({ n1, n2 }))
          continue;

        const auto& vars1 = node_to_vars[n1];
        const auto& vars2 = node_to_vars[n2];

        std::set<size_t> shared = tools_stl::intersect(vars1, vars2);

        if (shared.size() == vars1.size() - 1) {
          std::set<size_t> diff1 = tools_stl::set_diff(vars1, shared),
                           diff2 = tools_stl::set_diff(vars2, shared);
          size_t a = *diff1.begin();
          size_t b = *diff2.begin();
          candidates.emplace_back(n1, n2, Edge(a, b, shared));
        }
      }
    }

    // Add missing edges greedily
    while (tree.edges.size() < target_edges && !candidates.empty()) {
      auto candidate = std::move(candidates.back());
      candidates.pop_back();
      size_t n1 = std::get<0>(candidate);
      size_t n2 = std::get<1>(candidate);
      const Edge& e = std::get<2>(candidate);
      if (uf.find(n1) == uf.find(n2))
        continue;
      tree.edges.emplace_back(n1, n2, e);
      // std::cout << "Added edge (" << std::get<0>(e) << ", " << std::get<1>(e)
      //           << "; ";
      // for (const auto& v : std::get<2>(e)) {
      //   std::cout << v << " ";
      // }
      // std::cout << ")" << std::endl;
      uf.unite(n1, n2);
    }

    if (tree.edges.size() < target_edges) {
      throw std::runtime_error("Could not fill tree " + std::to_string(t));
    }

    // Recompute degrees
    tree.degrees.clear();
    for (const auto& e : tree.edges) {
      tree.degrees[std::get<0>(e)]++;
      tree.degrees[std::get<1>(e)]++;
    }
  }
}

} // namespace vinecopulib
