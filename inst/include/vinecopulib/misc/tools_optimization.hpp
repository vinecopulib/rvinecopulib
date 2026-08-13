// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <Eigen/Dense>
#include <functional>

namespace vinecopulib {

//! @brief Utilities for numerical optimization.
//!
//! A small optimizer used for parametric maximum likelihood: Brent's
//! derivative-free bracketing search for one-dimensional problems, and BFGS
//! for higher dimensions. For BFGS, bound constraints are handled by
//! optimizing over an unconstrained space and mapping back with
//! tools_transforms, so the objective always sees natural parameters.
namespace tools_optimization {

//! @brief Controls for the BFGS optimizer.
struct Controls
{
  size_t maxeval = 1000; //!< maximal number of objective evaluations per run
  double gtol = 1e-6;    //!< convergence: max-norm of the (eta-)gradient
  double ftol = 1e-9;    //!< convergence: relative change in the objective
  double xtol = 1e-8;    //!< convergence: relative change in eta
  double eta_max = 30.0; //!< clamp on `|eta|` for bounded coordinates
};

//! @brief The objective in the natural parameter space.
//!
//! Returns the value `f(theta)`. If `grad` has non-zero size on entry (equal
//! to `theta`), it must be filled with the gradient `df/dtheta`; the optimizer
//! passes an empty `grad` when only the value is needed (1-d Brent search and
//! backtracked line-search trials) and a sized `grad` where the gradient is
//! consumed.
using Objective =
  std::function<double(const Eigen::VectorXd& theta, Eigen::VectorXd& grad)>;

//! @brief A box-constrained maximizer: Brent in 1-d, BFGS with automatic
//! bound transforms otherwise.
class Optimizer
{
public:
  Optimizer() = default;

  explicit Optimizer(Controls controls);

  //! @brief maximizes `objective` over the box `[lower_bounds, upper_bounds]`.
  //!
  //! @param initial_parameters starting values, in natural coordinates.
  //! @param lower_bounds lower bounds (finite or `-inf`).
  //! @param upper_bounds upper bounds (finite or `+inf`).
  //! @param objective the objective (value and gradient), see Objective.
  //! @return the maximizing parameters, in natural coordinates.
  Eigen::VectorXd optimize(const Eigen::VectorXd& initial_parameters,
                           const Eigen::VectorXd& lower_bounds,
                           const Eigen::VectorXd& upper_bounds,
                           const Objective& objective);

  //! @brief the number of objective evaluations (cumulative over calls).
  size_t get_objective_calls() const;

  //! @brief the objective value at the maximum found by the last call.
  double get_objective_max() const;

private:
  void check_parameters_size(const Eigen::VectorXd& initial_parameters,
                             const Eigen::VectorXd& lower_bounds,
                             const Eigen::VectorXd& upper_bounds) const;

  Controls controls_{};
  size_t objective_calls_{ 0 };
  double objective_max_{ 0 };
};

} // namespace tools_optimization

} // namespace vinecopulib

#include <vinecopulib/misc/implementation/tools_optimization.ipp>
