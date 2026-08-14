// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <Eigen/Dense>
#include <cmath>
#include <limits>
#include <vector>

namespace vinecopulib {

//! @brief Tools for reparameterizing bounded parameters as unconstrained
//! optimization variables.
//!
//! The optimizer works in an unconstrained space `eta`; a `Transform` maps it
//! to a natural, bound-constrained parameter `theta` and back. The kind of
//! transform is inferred purely from the parameter bounds, so the statistical
//! model never has to know about it.
namespace tools_transforms {

//! @brief The kind of per-coordinate reparameterization, inferred from bounds.
enum class TransformKind
{
  Identity,      //!< unbounded `(-inf, inf)`: `theta = eta`
  SoftplusLower, //!< lower bounded `[a, inf)`: softplus
  SoftplusUpper, //!< upper bounded `(-inf, b]`: mirrored softplus
  Logistic,      //!< interval bounded `[a, b]`: scaled logistic
  Fixed          //!< degenerate `a == b`: coordinate held constant
};

//! @brief A smooth bijection between an unconstrained variable `eta` and a
//! bounded natural parameter `theta`.
//!
//! The kind is inferred from the (finite/infinite) bounds; the struct owns the
//! forward map, its inverse, and the derivative `dtheta/deta` used for the
//! chain rule. All maps are numerically stable and keep `theta` inside
//! `[a, b]`.
struct Transform
{
  TransformKind kind = TransformKind::Identity;
  double a = -std::numeric_limits<double>::infinity(); //!< lower bound
  double b = std::numeric_limits<double>::infinity();  //!< upper bound
  double eta_max = 30.0; //!< clamp on `|eta|` for bounded kinds

  //! @brief forward map; the result lies in `[a, b]`.
  double to_theta(double eta) const;

  //! @brief inverse map; `theta` is clamped strictly inside `[a, b]` first, so
  //! the result is always finite.
  double to_eta(double theta) const;

  //! @brief derivative `dtheta/deta` evaluated at `eta`.
  double dtheta_deta(double eta) const;

  //! @brief clamps `eta` to the transform's usable range (a no-op for
  //! `Identity`; one-sided for softplus; `[-eta_max, eta_max]` for `Logistic`).
  double clamp(double eta) const;
};

//! @brief A per-coordinate collection of `Transform`s for a parameter vector.
class ParameterTransform
{
public:
  ParameterTransform() = default;

  //! @brief builds one `Transform` per coordinate from the bounds.
  //! @param lower_bounds lower bounds (finite or `-inf`).
  //! @param upper_bounds upper bounds (finite or `+inf`).
  //! @param eta_max clamp on `|eta|` for bounded coordinates.
  static ParameterTransform from_bounds(const Eigen::VectorXd& lower_bounds,
                                        const Eigen::VectorXd& upper_bounds,
                                        double eta_max = 30.0);

  //! @brief maps an unconstrained vector to natural parameters.
  Eigen::VectorXd to_theta(const Eigen::VectorXd& eta) const;

  //! @brief maps natural parameters to an unconstrained vector.
  Eigen::VectorXd to_eta(const Eigen::VectorXd& theta) const;

  //! @brief per-coordinate derivative `dtheta/deta` (the Jacobian diagonal).
  Eigen::VectorXd dtheta_deta(const Eigen::VectorXd& eta) const;

  //! @brief clamps each coordinate of `eta` to its transform's usable range.
  Eigen::VectorXd clamp_eta(const Eigen::VectorXd& eta) const;

  //! @brief the number of coordinates.
  size_t size() const { return transforms_.size(); }

private:
  std::vector<Transform> transforms_;
};

} // namespace tools_transforms

} // namespace vinecopulib

#include <vinecopulib/misc/implementation/tools_transforms.ipp>
