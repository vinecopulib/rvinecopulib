// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <boost/math/tools/minima.hpp>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vinecopulib/misc/tools_transforms.hpp>

namespace vinecopulib {

namespace tools_optimization {

inline Optimizer::Optimizer(Controls controls)
  : controls_(controls)
{
}

//! @brief maximizes `objective` over a box.
//!
//! One-dimensional problems with finite bounds use Brent's derivative-free
//! bracketing search (it needs only ~10 value evaluations on the narrowed
//! intervals used for fitting, which a gradient method cannot beat). Higher
//! dimensions use unconstrained BFGS: the problem is transformed via
//! tools_transforms (identity / softplus / logistic inferred from the bounds)
//! and minimized in `eta`-space with a backtracking-Armijo line search; the
//! chain rule maps the natural-space gradient supplied by `objective` to
//! `eta`-space. The best point seen is always returned, so a bad line-search
//! step cannot regress the result.
inline Eigen::VectorXd
Optimizer::optimize(const Eigen::VectorXd& initial_parameters,
                    const Eigen::VectorXd& lower_bounds,
                    const Eigen::VectorXd& upper_bounds,
                    const Objective& objective)
{
  check_parameters_size(initial_parameters, lower_bounds, upper_bounds);

  const Eigen::Index n = initial_parameters.size();

  if ((n == 1) && std::isfinite(lower_bounds(0)) &&
      std::isfinite(upper_bounds(0))) {
    // 1-d: Brent's method on the (slightly shrunk) bracket; value-only.
    const double eps = 1e-6;
    Eigen::VectorXd no_grad; // empty: the objective skips the gradient
    auto f = [&](double x) {
      ++objective_calls_;
      return -objective(Eigen::VectorXd::Constant(1, x), no_grad);
    };
    auto result = boost::math::tools::brent_find_minima(
      f, lower_bounds(0) + eps, upper_bounds(0) - eps, 20);
    objective_max_ = -result.second;
    return Eigen::VectorXd::Constant(1, result.first);
  }
  const auto transform = tools_transforms::ParameterTransform::from_bounds(
    lower_bounds, upper_bounds, controls_.eta_max);

  size_t n_eval =
    0; // per-call evaluation budget (objective_calls_ is cumulative)

  // Evaluates m(eta) = -f(theta(eta)) and, when `need_grad`, its eta-gradient.
  // Returns false if the objective throws or is not finite at `eta`.
  auto eval = [&](const Eigen::VectorXd& eta,
                  bool need_grad,
                  double& m,
                  Eigen::VectorXd& grad_eta) -> bool {
    Eigen::VectorXd theta = transform.to_theta(eta);
    Eigen::VectorXd grad_theta =
      need_grad ? Eigen::VectorXd::Zero(n) : Eigen::VectorXd();
    ++n_eval;
    ++objective_calls_;
    double f;
    try {
      f = objective(theta, grad_theta);
    } catch (...) {
      return false;
    }
    if (!std::isfinite(f)) {
      return false;
    }
    m = -f;
    if (need_grad) {
      if (grad_theta.allFinite()) {
        // chain rule: d(-f)/deta = -(dtheta/deta) .* df/dtheta
        grad_eta =
          -(transform.dtheta_deta(eta).array() * grad_theta.array()).matrix();
      } else {
        grad_eta = Eigen::VectorXd::Zero(n);
      }
    }
    return true;
  };

  Eigen::VectorXd eta = transform.clamp_eta(transform.to_eta(
    initial_parameters.cwiseMax(lower_bounds).cwiseMin(upper_bounds)));

  double m;
  Eigen::VectorXd grad(n);
  if (!eval(eta, true, m, grad)) {
    // objective unusable at the start; report the (clamped) initial point.
    objective_max_ = -std::numeric_limits<double>::infinity();
    return transform.to_theta(eta);
  }

  Eigen::VectorXd best_eta = eta;
  double best_m = m;

  const Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(n, n);
  // initial inverse Hessian, scaled so the first step has magnitude ~1 in
  // eta-space; the raw gradient scales with the number of observations, and
  // an O(n) first step would only burn line-search halvings. Overwritten by
  // the Nocedal-Wright rescale after the first accepted step.
  Eigen::MatrixXd H = Id / std::max(1.0, grad.cwiseAbs().maxCoeff());
  const double eps_c = 1e-10;
  const double c1 = 1e-4;
  const double alpha_min = 1e-12;
  bool did_scale = false;

  while (n_eval < controls_.maxeval) {
    // search direction; reset to steepest descent if not a descent direction.
    Eigen::VectorXd p = -(H * grad);
    double gp = grad.dot(p);
    if (!p.allFinite() || gp >= 0.0) {
      H = Id;
      p = -grad;
      gp = grad.dot(p);
    }
    if (gp >= 0.0) {
      break; // gradient is (numerically) zero: stationary point
    }

    // backtracking-Armijo line search. The first trial (usually accepted)
    // also computes the gradient, so the common case needs no separate
    // re-evaluation at the accepted point; backtracked trials are value-only.
    double alpha = 1.0;
    Eigen::VectorXd eta_new;
    Eigen::VectorXd grad_new(n);
    double m_new = m;
    bool ls_ok = false;
    bool have_grad = false;
    bool first_trial = true;
    while (n_eval < controls_.maxeval) {
      eta_new = transform.clamp_eta(eta + alpha * p);
      double m_trial = m;
      Eigen::VectorXd grad_trial;
      if (eval(eta_new, first_trial, m_trial, grad_trial) &&
          m_trial <= m + c1 * alpha * gp) {
        m_new = m_trial;
        if (first_trial) {
          grad_new = grad_trial;
          have_grad = true;
        }
        ls_ok = true;
        break;
      }
      first_trial = false;
      alpha *= 0.5;
      if (alpha < alpha_min) {
        break;
      }
    }
    if (!ls_ok) {
      break; // line search failed: stop and return the best point seen
    }

    // gradient at the accepted iterate (only needed after backtracking).
    if (!have_grad) {
      double m_unused;
      if (!eval(eta_new, true, m_unused, grad_new)) {
        break;
      }
    }

    const Eigen::VectorXd s = eta_new - eta;
    const Eigen::VectorXd y = grad_new - grad;
    const double sy = s.dot(y);

    // Nocedal & Wright (6.20): scale the initial inverse Hessian once.
    if (!did_scale) {
      const double yy = y.dot(y);
      if (sy > 0.0 && yy > 0.0) {
        H = (sy / yy) * Id;
      }
      did_scale = true;
    }

    // BFGS inverse-Hessian update, skipped when the curvature condition fails
    // (keeps H positive definite under noisy / non-smooth gradients).
    if (sy > eps_c * s.norm() * y.norm()) {
      const double rho = 1.0 / sy;
      const Eigen::MatrixXd V = Id - rho * (s * y.transpose());
      H = V * H * V.transpose() + rho * (s * s.transpose());
    }

    if (m_new < best_m) {
      best_m = m_new;
      best_eta = eta_new;
    }

    // convergence tests (all relative / scale invariant).
    const double gnorm = grad_new.cwiseAbs().maxCoeff();
    const double dstep = s.cwiseAbs().maxCoeff();
    const bool conv_g = gnorm <= controls_.gtol * (1.0 + std::fabs(m_new));
    const bool conv_x =
      dstep <= controls_.xtol * (1.0 + eta_new.cwiseAbs().maxCoeff());
    const bool conv_f =
      std::fabs(m_new - m) <= controls_.ftol * (1.0 + std::fabs(m_new));

    eta = eta_new;
    m = m_new;
    grad = grad_new;

    if (conv_g || conv_x || conv_f) {
      break;
    }
  }

  objective_max_ = -best_m;
  return transform.to_theta(best_eta);
}

//! @brief the number of objective evaluations (cumulative over calls).
inline size_t
Optimizer::get_objective_calls() const
{
  return objective_calls_;
}

//! @brief the objective value at the maximum found by the last call.
inline double
Optimizer::get_objective_max() const
{
  return objective_max_;
}

//! @brief checks that the sizes of the parameters and bounds match.
inline void
Optimizer::check_parameters_size(const Eigen::VectorXd& initial_parameters,
                                 const Eigen::VectorXd& lower_bounds,
                                 const Eigen::VectorXd& upper_bounds) const
{
  if (initial_parameters.size() != upper_bounds.size()) {
    throw std::runtime_error(
      "initial parameters and and bounds must have same size.");
  }
  if (lower_bounds.size() != upper_bounds.size()) {
    throw std::runtime_error("lower and upper bounds must have same size.");
  }
  if (lower_bounds.size() < 1) {
    throw std::runtime_error("n_parameters should be larger than 0.");
  }
}

} // namespace tools_optimization

} // namespace vinecopulib
