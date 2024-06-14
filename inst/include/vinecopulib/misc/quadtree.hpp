#include <vector>
#include <random>
#include <unordered_map>
#include <array>

namespace vinecopulib {

struct Point {
  double x, y;
};

// A set of points with O(1) insert, remove, and uniform sampling;
// inspired by https://gist.github.com/matovitch/1c262d2d870024ca0b81
class PointSet
{
public:
  size_t size() { return vector_.size(); }

  void insert(const Point& p)
  {
    if (ptr_set_.find(&p) == ptr_set_.end()) {  // point not already in the set
      vector_.push_back(p);
      ptr_set_.insert(&(*vector_.crbegin()));
    }
  }

  void remove(const Point& p)
  {
    auto el = ptr_set_.find(&p);
    if (el != ptr_set_.end()) {
      // first replace element by last entry of vector_, then pop last entry
      if (*el != &(*vector_.crbegin())) {
        *(const_cast<Point*>(*el)) = *vector_.rbegin();
      }
      vector_.pop_back();
    }
  }

  const Point& sample() const
  {
    auto n = vector_.size();
    if (n == 0) {
      throw std::runtime_error("Cannot sample from empty set.");
    }
    std::mt19937 rng(std::random_device{}());
    auto distribution = std::uniform_int_distribution<std::size_t>(0, n - 1);
    return vector_[distribution(rng) % vector_.size()];
  }

private:
  std::unordered_set<const Point*> ptr_set_;
  std::vector<Point> vector_;
};


struct BoundingBox {
  double x_min, y_min, x_mid, y_mid, x_max, y_max;

  BoundingBox(double x_min, double y_min, double x_max, double y_max)
    : x_min(x_min), y_min(y_min), x_max(x_max), y_max(y_max),
      x_mid(0.5 * (x_min + x_max)), y_mid(0.5 * (y_min + y_max))
  {}

  bool contains(const BoundingBox& other) const {
    return other.x_min >= x_min && other.x_max <= x_max &&
      other.y_min >= y_min && other.y_max <= y_max;
  }

  bool intersects(const BoundingBox& other) const {
    return !(other.x_min > x_max || other.x_max < x_min ||
             other.y_min > y_max || other.y_max < y_min);
  }
};

// Fixed depth quadtree class
class QuadTree {
private:
  struct Node {
    BoundingBox boundary;
    PointSet points;

    std::array<Node*, 4> children;
    Node(const BoundingBox& boundary) : boundary(boundary) {}
  };

  void construct_children(Node* node, uint16_t depth, int& node_idx) {
    if (depth >= depth_) return;

    const auto& b = node->boundary;
    nodes_.emplace_back(BoundingBox{b.x_min, b.y_min, b.x_mid, b.y_mid});
    nodes_.emplace_back(BoundingBox{b.x_mid, b.y_min, b.x_max, b.y_mid});
    nodes_.emplace_back(BoundingBox{b.x_min, b.y_mid, b.x_mid, b.y_max});
    nodes_.emplace_back(BoundingBox{b.x_mid, b.y_mid, b.x_max, b.y_max});

    for (int i = 0; i < 4; ++i) {
      node->children[i] = &nodes_[node_idx++];
      auto c = node->children[i];
    }

    for (auto& child : node->children) {
      construct_children(child, depth + 1, node_idx);
    }
  }

  Node* find_child_with_point(const Node* node, const Point& point)
  {
    if (point.x <= node->boundary.x_mid && point.y <= node->boundary.y_mid) {
      return node->children[0];
    } else if (point.y <= node->boundary.y_mid) {
      return node->children[1];
    } else if (point.x <= node->boundary.x_mid) {
      return node->children[2];
    } else {
      return node->children[3];
    }
  }

  std::vector<Node> nodes_;
  std::vector<Node*> terminal_nodes_;
  uint16_t depth_;

public:
  QuadTree(const BoundingBox& boundary, int depth)
    : depth_(depth)
  {
    // geometric series sum formula 4^0 + 4^1 + ... + 4^depth
    size_t num_nodes = (std::pow(4, depth + 1) - 1) / 3;
    nodes_.reserve(num_nodes);
    terminal_nodes_.reserve(num_nodes);

    nodes_.emplace_back(boundary);
    int node_idx = 1;
    construct_children(&nodes_[0], 0, node_idx);
  }

  void insert(const Point& point)
  {
    Node* node = &nodes_[0]; // start at root
    for (int d = 0; d <= depth_; ++d) {
      node->points.insert(point);
      node = find_child_with_point(node, point);
    }
  }

  void remove(const Point& point) {
    Node* node = &nodes_[0]; // start at root
    for (int d = 0; d <= depth_; ++d) {
      node->points.remove(point);
      node = find_child_with_point(node, point);
    }
  }

  // Helper function to recursively calculate points within the range
  void find_terminal_nodes(Node* node, const BoundingBox& range, int depth)
  {
    if (depth == 0) {
      terminal_nodes_.clear(); // starting new search
    }

    // If the node is outside the range, no points are within the range
    if (!range.intersects(node->boundary)) {
      return;
    }

    // If the node is fully contained in the range, all its points are
    // If at deepest level, return the node even though it only intersects
    if (range.contains(node->boundary) || (depth == depth_)) {
      terminal_nodes_.push_back(node);
      return;
    }

    // Otherwise, check each child node
    for (auto& child : node->children) {
      find_terminal_nodes(child, range, depth + 1);
    }

  };

  Point sample(const BoundingBox& range) {
    // Find nodes at highest possible level fully contained in the range,
    // starting at root
    find_terminal_nodes(&nodes_[0], range, 0);

    int total_points_in_range = 0;
    for (auto node : terminal_nodes_) {
      auto c = node;
      total_points_in_range += node->points.size();
    }
    if (total_points_in_range == 0) {
      throw std::runtime_error("No points in the specified range.");
    }

    // Select a child node based on the counts
    std::uniform_int_distribution<int> dist(0, total_points_in_range - 1);
    std::mt19937 rng(std::random_device{}());
    int r = dist(rng);

    Node* sampled_node;
    for (auto node : terminal_nodes_) {
      if (r < node->points.size()) {
        sampled_node = node;
        break;
      }
      r -= node->points.size();
    }

    return sampled_node->points.sample();
  }
};

} // namespace vinecopulib
