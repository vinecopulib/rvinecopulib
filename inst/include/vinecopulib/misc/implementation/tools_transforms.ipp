// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <algorithm>

namespace vinecopulib {

namespace tools_transforms {

//! @brief numerically stable logistic sigmoid `1 / (1 + exp(-x))`.
inline double
sigmoid(double x)
{
  return x >= 0.0 ? 1.0 / (1.0 + std::exp(-x))
                  : std::exp(x) / (1.0 + std::exp(x));
}

//! @brief numerically stable softplus `log(1 + exp(x))`.
inline double
softplus(double x)
{
  return std::max(x, 0.0) + std::log1p(std::exp(-std::fabs(x)));
}

//! @brief inverse softplus, `log(exp(y) - 1)` for `y >= 0`.
inline double
softplus_inv(double y)
{
  y = std::max(y, std::numeric_limits<double>::min());
  return y + std::log1p(-std::exp(-y));
}

//! @brief forward map from the unconstrained `eta` to the natural `theta`.
inline double
Transform::to_theta(double eta) const
{
  switch (kind) {
    case TransformKind::SoftplusLower:
      return a + softplus(eta);
    case TransformKind::SoftplusUpper:
      return b - softplus(-eta);
    case TransformKind::Logistic: {
      double theta = a + (b - a) * sigmoid(eta);
      return std::min(std::max(theta, a), b); // guard floating-point rounding
    }
    case TransformKind::Fixed:
      return a;
    case TransformKind::Identity:
    default:
      return eta;
  }
}

//! @brief inverse map from the natural `theta` to the unconstrained `eta`.
inline double
Transform::to_eta(double theta) const
{
  switch (kind) {
    case TransformKind::SoftplusLower:
      return std::max(softplus_inv(theta - a), -eta_max);
    case TransformKind::SoftplusUpper:
      return std::min(-softplus_inv(b - theta), eta_max);
    case TransformKind::Logistic: {
      double eps_z = sigmoid(-eta_max); // keeps eta within [-eta_max, eta_max]
      double p = (theta - a) / (b - a);
      p = std::min(std::max(p, eps_z), 1.0 - eps_z);
      double eta = std::log(p) - std::log1p(-p);
      return std::min(std::max(eta, -eta_max), eta_max);
    }
    case TransformKind::Fixed:
      return 0.0;
    case TransformKind::Identity:
    default:
      return theta;
  }
}

//! @brief derivative `dtheta/deta`.
inline double
Transform::dtheta_deta(double eta) const
{
  switch (kind) {
    case TransformKind::SoftplusLower:
      return sigmoid(eta);
    case TransformKind::SoftplusUpper:
      return sigmoid(-eta);
    case TransformKind::Logistic:
      return (b - a) * sigmoid(eta) * sigmoid(-eta);
    case TransformKind::Fixed:
      return 0.0;
    case TransformKind::Identity:
    default:
      return 1.0;
  }
}

//! @brief clamps `eta` to the transform's usable range.
inline double
Transform::clamp(double eta) const
{
  switch (kind) {
    case TransformKind::SoftplusLower:
      return std::max(eta, -eta_max); // unbounded above
    case TransformKind::SoftplusUpper:
      return std::min(eta, eta_max); // unbounded below
    case TransformKind::Logistic:
      return std::min(std::max(eta, -eta_max), eta_max);
    case TransformKind::Fixed:
      return 0.0;
    case TransformKind::Identity:
    default:
      return eta; // unbounded both ways
  }
}

//! @brief builds one `Transform` per coordinate from the bounds.
inline ParameterTransform
ParameterTransform::from_bounds(const Eigen::VectorXd& lower_bounds,
                                const Eigen::VectorXd& upper_bounds,
                                double eta_max)
{
  ParameterTransform pt;
  const Eigen::Index n = lower_bounds.size();
  pt.transforms_.resize(static_cast<size_t>(n));
  for (Eigen::Index k = 0; k < n; ++k) {
    Transform t;
    t.eta_max = eta_max;
    t.a = lower_bounds(k);
    t.b = upper_bounds(k);
    const bool lo_finite = std::isfinite(t.a);
    const bool hi_finite = std::isfinite(t.b);
    if (lo_finite && hi_finite) {
      t.kind =
        (t.b - t.a > 0.0) ? TransformKind::Logistic : TransformKind::Fixed;
    } else if (lo_finite) {
      t.kind = TransformKind::SoftplusLower;
    } else if (hi_finite) {
      t.kind = TransformKind::SoftplusUpper;
    } else {
      t.kind = TransformKind::Identity;
    }
    pt.transforms_[static_cast<size_t>(k)] = t;
  }
  return pt;
}

inline Eigen::VectorXd
ParameterTransform::to_theta(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd theta(eta.size());
  for (Eigen::Index k = 0; k < eta.size(); ++k) {
    theta(k) = transforms_[static_cast<size_t>(k)].to_theta(eta(k));
  }
  return theta;
}

inline Eigen::VectorXd
ParameterTransform::to_eta(const Eigen::VectorXd& theta) const
{
  Eigen::VectorXd eta(theta.size());
  for (Eigen::Index k = 0; k < theta.size(); ++k) {
    eta(k) = transforms_[static_cast<size_t>(k)].to_eta(theta(k));
  }
  return eta;
}

inline Eigen::VectorXd
ParameterTransform::dtheta_deta(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd d(eta.size());
  for (Eigen::Index k = 0; k < eta.size(); ++k) {
    d(k) = transforms_[static_cast<size_t>(k)].dtheta_deta(eta(k));
  }
  return d;
}

inline Eigen::VectorXd
ParameterTransform::clamp_eta(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd out(eta.size());
  for (Eigen::Index k = 0; k < eta.size(); ++k) {
    out(k) = transforms_[static_cast<size_t>(k)].clamp(eta(k));
  }
  return out;
}

} // namespace tools_transforms

} // namespace vinecopulib
