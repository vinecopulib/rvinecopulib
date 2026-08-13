// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/misc/tools_optimization.hpp>
#include <vinecopulib/misc/tools_stats.hpp>
#include <vinecopulib/misc/tools_stl.hpp>
#include <wdm/eigen.hpp>

namespace vinecopulib {
inline Eigen::MatrixXd
ParBicop::get_parameters() const
{
  return parameters_;
}

inline Eigen::MatrixXd
ParBicop::get_parameters_lower_bounds() const
{
  return parameters_lower_bounds_;
}

inline Eigen::MatrixXd
ParBicop::get_parameters_upper_bounds() const
{
  return parameters_upper_bounds_;
}

inline void
ParBicop::set_parameters(const Eigen::MatrixXd& parameters)
{
  check_parameters(parameters);
  parameters_ = parameters;
}

inline void
ParBicop::flip()
{
  // Most parametric families can be flipped by changing the rotation.
  // This is done in Bicop::flip() directly. All other families need to
  // override this method.
}

// calculate number of parameters
inline double
ParBicop::get_npars() const
{
  // indepence copula has no parameters
  if (family_ == BicopFamily::indep) {
    return 0.0;
  }
  // otherwise, return length of parameter vector
  return static_cast<double>(parameters_.size());
}

inline void
ParBicop::set_npars(const double&)
{
  // does nothing
}

// fit
inline void
ParBicop::fit(const Eigen::MatrixXd& data,
              std::string method,
              double,
              size_t,
              const Eigen::VectorXd& weights)
{
  // for independence copula we don't have to do anything
  if (family_ == BicopFamily::indep) {
    set_loglik(0.0);
    return;
  }

  check_fit_method(method);
  double tau = wdm::wdm(data, "tau", weights)(0, 1);

  // for method itau and one-parameter families we don't need to optimize
  int npars = static_cast<int>(get_npars()) - (method == "itau");
  if (npars == 0) {
    set_parameters(tau_to_parameters(tau));
    set_loglik(loglik(data, weights));
    return;
  }

  // Set bounds and starting values
  auto lb = get_parameters_lower_bounds();
  auto ub = get_parameters_upper_bounds();
  adjust_parameters_bounds(lb, ub, tau, method);
  auto initial_parameters = get_start_parameters(winsorize_tau(tau));

  // find (pseudo-) mle. the optimizer works in an unconstrained space and only
  // ever sees natural parameters; it needs both the value and the gradient.
  // continuous data uses the analytic (or finite-difference) score leaves via
  // `logpdf_deriv_raw`; discrete data (whose likelihood uses h-function
  // differences) falls back to a finite difference of the total log-likelihood.
  const bool continuous = (var_types_ == std::vector<std::string>{ "c", "c" });
  tools_optimization::Objective objective;
  if (method == "mle") {
    objective = [&data, &weights, this, continuous, lb, ub](
                  const Eigen::VectorXd& pars, Eigen::VectorXd& grad) {
      this->set_parameters(pars);
      double val = this->loglik(data, weights);
      if (grad.size() > 0) {
        if (continuous) {
          Eigen::ArrayXd w = Eigen::ArrayXd::Ones(data.rows());
          if (weights.size() > 0)
            w = weights.array();
          for (Eigen::Index k = 0; k < pars.size(); ++k) {
            Eigen::ArrayXd dw =
              this
                ->logpdf_deriv_raw(
                  data, pars.transpose(), "par" + std::to_string(k + 1))
                .array() *
              w;
            grad(k) = dw.isFinite().select(dw, 0.0).sum();
          }
        } else {
          grad = this->fd_grad(
            [&data, &weights, this](const Eigen::VectorXd& p) {
              this->set_parameters(p);
              return this->loglik(data, weights);
            },
            pars,
            lb,
            ub);
        }
      }
      return val;
    };
  } else {
    // profile likelihood: the first parameter is held fixed, only the second
    // is optimized.
    set_parameters(initial_parameters);
    initial_parameters(0) = initial_parameters(1);
    initial_parameters.conservativeResize(1);
    objective = [&data, &weights, this, continuous, lb, ub](
                  const Eigen::VectorXd& pars, Eigen::VectorXd& grad) {
      Eigen::VectorXd newpars(2);
      newpars(0) = this->get_parameters()(0);
      newpars(1) = pars(0);
      this->set_parameters(newpars);
      double val = this->loglik(data, weights);
      if (grad.size() > 0) {
        if (continuous) {
          Eigen::ArrayXd w = Eigen::ArrayXd::Ones(data.rows());
          if (weights.size() > 0)
            w = weights.array();
          Eigen::ArrayXd dw =
            this->logpdf_deriv_raw(data, newpars.transpose(), "par2").array() *
            w;
          grad(0) = dw.isFinite().select(dw, 0.0).sum();
        } else {
          grad = this->fd_grad(
            [&data, &weights, this](const Eigen::VectorXd& p) {
              Eigen::VectorXd np(2);
              np(0) = this->get_parameters()(0);
              np(1) = p(0);
              this->set_parameters(np);
              return this->loglik(data, weights);
            },
            pars,
            lb,
            ub);
        }
      }
      return val;
    };
  }

  tools_optimization::Optimizer optimizer;
  auto newpars = optimizer.optimize(initial_parameters, lb, ub, objective);

  // check if fit is reasonable, otherwise increase search interval
  // and refit
  if (tools_stl::is_member(family_, bicop_families::one_par) &&
      (optimizer.get_objective_max() < -0.1)) {
    newpars = optimizer.optimize(initial_parameters,
                                 get_parameters_lower_bounds(),
                                 get_parameters_upper_bounds(),
                                 objective);
  }

  // finalize fitted model
  if (method == "itau") {
    // only the second parameter has been optimized
    newpars.conservativeResize(2);
    std::swap(newpars(0), newpars(1));
    newpars(0) = get_parameters()(0);
  }

  set_parameters(newpars);
  set_loglik(optimizer.get_objective_max());
}

//! ensures that starting values are sufficiently separated from bounds
//! @param tau Kendall's tau
inline double
ParBicop::winsorize_tau(double tau) const
{
  double sign = 1.0;
  if (tau < 0) {
    sign = -1.0;
  }
  if (std::abs(tau) < 0.01) {
    tau = 0.01 * sign;
  } else if (std::abs(tau) > 0.9) {
    tau = 0.9 * sign;
  }
  return tau;
}

//! adjusts parameter bounds for better search intervals.
inline void
ParBicop::adjust_parameters_bounds(Eigen::MatrixXd& lb,
                                   Eigen::MatrixXd& ub,
                                   const double& tau,
                                   const std::string& method)
{
  if (method == "itau") {
    // for pseudo mle, we can ignore the first parameter
    lb(0) = lb(1);
    ub(0) = ub(1);
    lb.conservativeResize(1, 1);
    ub.conservativeResize(1, 1);
    if (family_ == BicopFamily::student) {
      // the df parameter doesn't need to be estimated as accurately
      ub(0) = 15;
    }
  }

  // refine search interval for Brent algorithm
  double eps = (var_types_ == std::vector<std::string>{ "c", "c" }) ? 0.1 : 0.6;
  if (tools_stl::is_member(family_, bicop_families::one_par)) {
    auto lb2 = lb;
    auto ub2 = ub;
    if (tools_stl::is_member(family_, bicop_families::rotationless)) {
      lb = tau_to_parameters(std::max(tau - eps, -0.99));
      ub = tau_to_parameters(std::min(tau + eps, 0.99));
    } else {
      lb = tau_to_parameters(std::max(std::fabs(tau) - eps, 1e-10));
      ub = tau_to_parameters(std::min(std::fabs(tau) + eps, 0.95));
    }
    // make sure that parameter bounds are respected
    lb = lb2.cwiseMax(lb);
    ub = ub2.cwiseMin(ub);
  }

  if (family_ == BicopFamily::tawn) {
    Eigen::VectorXd lb2(3), ub2(3);
    lb2 << 0.3, 0.3, 1.5;
    ub2 << 1, 1, 7;
    lb = lb2;
    ub = ub2;
  }
}

//! Sanity checks
//! @{
inline void
ParBicop::check_parameters(const Eigen::MatrixXd& parameters)
{
  check_parameters_size(parameters);
  check_parameters_lower(parameters);
  check_parameters_upper(parameters);
}

inline void
ParBicop::check_parameters_size(const Eigen::MatrixXd& parameters)
{
  // Validate rows and cols unconditionally: a same-size-but-transposed shape
  // (e.g. 1 x p vs p x 1) would otherwise slip through and reach the
  // coefficient-wise bound comparisons in check_parameters_lower/_upper with
  // mismatched shapes (an out-of-bounds read under NDEBUG).
  if (parameters.rows() != parameters_.rows()) {
    std::stringstream message;
    message << "parameters have has wrong number of rows "
            << "for " << get_family_name() << " copula; "
            << "expected: " << parameters_.rows() << ", "
            << "actual: " << parameters.rows() << std::endl;
    throw std::runtime_error(message.str().c_str());
  }
  if (parameters.cols() != parameters_.cols()) {
    std::stringstream message;
    message << "parameters have wrong number of columns "
            << "for " << get_family_name() << " copula; "
            << "expected: " << parameters_.cols() << ", "
            << "actual: " << parameters.cols() << std::endl;
    throw std::runtime_error(message.str().c_str());
  }
}

inline void
ParBicop::check_parameters_lower(const Eigen::MatrixXd& parameters)
{
  if (parameters_lower_bounds_.size() > 0) {
    std::stringstream message;
    if ((parameters.array() < parameters_lower_bounds_.array()).any()) {
      message << "parameters exceed lower bound "
              << "for " << get_family_name() << " copula; " << std::endl
              << "bound:" << std::endl
              << parameters_lower_bounds_ << std::endl
              << "actual:" << std::endl
              << parameters << std::endl;
      throw std::runtime_error(message.str().c_str());
    }
  }
}

inline void
ParBicop::check_parameters_upper(const Eigen::MatrixXd& parameters)
{
  if (parameters_upper_bounds_.size() > 0) {
    std::stringstream message;
    if ((parameters.array() > parameters_upper_bounds_.array()).any()) {
      message << "parameters exceed upper bound "
              << "for " << get_family_name() << " copula; " << std::endl
              << "bound:" << std::endl
              << parameters_upper_bounds_ << std::endl
              << "actual:" << std::endl
              << parameters << std::endl;
      throw std::runtime_error(message.str().c_str());
    }
  }
}

inline void
ParBicop::check_fit_method(const std::string& method)
{
  if (!tools_stl::is_member(method, { "itau", "mle" })) {
    throw std::runtime_error("Method not implemented.");
  }

  if (method == "itau") {
    if (!tools_stl::is_member(family_, bicop_families::itau)) {
      throw std::runtime_error("itau method is not available for this family.");
    }
  }
}

//! @}

//! @name Finite-difference derivative fallbacks
//!
//! Central differences of the value leaves w.r.t. a single component; the
//! second-order versions difference the (possibly analytic) first-derivative
//! leaves. Steps are clipped to the parameter bounds / the unit interval,
//! with the effective step kept in the denominator.
//! @{

//! differentiates `f` w.r.t. component `comp` (0-based parameter index, `-1`
//! for the first argument, `-2` for the second) by central differences.
inline Eigen::VectorXd
ParBicop::fd_deriv(
  const std::function<Eigen::VectorXd(const Eigen::MatrixXd&,
                                      const Eigen::MatrixXd&)>& f,
  const Eigen::MatrixXd& u,
  const Eigen::MatrixXd& parameters,
  int comp)
{
  if (comp >= 0) {
    Eigen::MatrixXd par_plus = parameters;
    Eigen::MatrixXd par_minus = parameters;
    double lb = parameters_lower_bounds_(comp);
    double ub = parameters_upper_bounds_(comp);
    for (Eigen::Index i = 0; i < parameters.rows(); ++i) {
      double par = parameters(i, comp);
      double eps = 1e-4 * std::max(1.0, std::fabs(par));
      par_plus(i, comp) = std::min(par + eps, ub);
      par_minus(i, comp) = std::max(par - eps, lb);
    }
    Eigen::VectorXd diff = f(u, par_plus) - f(u, par_minus);
    if (parameters.rows() == 1) {
      return diff / (par_plus(0, comp) - par_minus(0, comp));
    }
    Eigen::ArrayXd eps = (par_plus.col(comp) - par_minus.col(comp)).array();
    return (diff.array() / eps).matrix();
  }

  // argument derivative; stay strictly inside the unit interval
  Eigen::Index col = (comp == -1) ? 0 : 1;
  Eigen::MatrixXd u_plus = u;
  Eigen::MatrixXd u_minus = u;
  u_plus.col(col) = (u.col(col).array() + 1e-5).min(1 - 1e-10);
  u_minus.col(col) = (u.col(col).array() - 1e-5).max(1e-10);
  Eigen::ArrayXd eps = u_plus.col(col).array() - u_minus.col(col).array();
  return ((f(u_plus, parameters) - f(u_minus, parameters)).array() / eps)
    .matrix();
}

//! central finite differences of a scalar objective `f` w.r.t. each
//! optimization variable, with steps clipped to `[lb, ub]`.
inline Eigen::VectorXd
ParBicop::fd_grad(const std::function<double(const Eigen::VectorXd&)>& f,
                  const Eigen::VectorXd& x,
                  const Eigen::VectorXd& lb,
                  const Eigen::VectorXd& ub)
{
  Eigen::VectorXd grad(x.size());
  for (Eigen::Index k = 0; k < x.size(); ++k) {
    double eps = 1e-4 * std::max(1.0, std::fabs(x(k)));
    double xp = std::min(x(k) + eps, ub(k));
    double xm = std::max(x(k) - eps, lb(k));
    Eigen::VectorXd xpv = x;
    Eigen::VectorXd xmv = x;
    xpv(k) = xp;
    xmv(k) = xm;
    double denom = xp - xm;
    grad(k) = denom > 0.0 ? (f(xpv) - f(xmv)) / denom : 0.0;
  }
  return grad;
}

inline Eigen::VectorXd
ParBicop::pdf_deriv_raw(const Eigen::MatrixXd& u,
                        const Eigen::MatrixXd& parameters,
                        const std::string& deriv)
{
  auto comps = tools_deriv::parse_components(deriv);
  auto f = [this](const Eigen::MatrixXd& uu, const Eigen::MatrixXd& pp) {
    return pdf_raw(uu, pp);
  };
  return fd_deriv(f, u, parameters, comps[0]);
}

inline Eigen::VectorXd
ParBicop::pdf_deriv2_raw(const Eigen::MatrixXd& u,
                         const Eigen::MatrixXd& parameters,
                         const std::string& deriv)
{
  // difference the first-derivative leaf w.r.t. the second component; this
  // uses the analytic first derivative when the family provides one
  auto comps = tools_deriv::parse_components(deriv);
  auto first = tools_deriv::comp_to_string(comps[0]);
  auto f = [this, first](const Eigen::MatrixXd& uu, const Eigen::MatrixXd& pp) {
    return pdf_deriv_raw(uu, pp, first);
  };
  return fd_deriv(f, u, parameters, comps[1]);
}

inline Eigen::VectorXd
ParBicop::hfunc1_deriv_raw(const Eigen::MatrixXd& u,
                           const Eigen::MatrixXd& parameters,
                           const std::string& deriv)
{
  auto comps = tools_deriv::parse_components(deriv);
  auto f = [this](const Eigen::MatrixXd& uu, const Eigen::MatrixXd& pp) {
    return hfunc1_raw(uu, pp);
  };
  return fd_deriv(f, u, parameters, comps[0]);
}

inline Eigen::VectorXd
ParBicop::hfunc1_deriv2_raw(const Eigen::MatrixXd& u,
                            const Eigen::MatrixXd& parameters,
                            const std::string& deriv)
{
  auto comps = tools_deriv::parse_components(deriv);
  auto first = tools_deriv::comp_to_string(comps[0]);
  auto f = [this, first](const Eigen::MatrixXd& uu, const Eigen::MatrixXd& pp) {
    return hfunc1_deriv_raw(uu, pp, first);
  };
  return fd_deriv(f, u, parameters, comps[1]);
}

inline Eigen::VectorXd
ParBicop::hfunc2_deriv_raw(const Eigen::MatrixXd& u,
                           const Eigen::MatrixXd& parameters,
                           const std::string& deriv)
{
  auto comps = tools_deriv::parse_components(deriv);
  auto f = [this](const Eigen::MatrixXd& uu, const Eigen::MatrixXd& pp) {
    return hfunc2_raw(uu, pp);
  };
  return fd_deriv(f, u, parameters, comps[0]);
}

inline Eigen::VectorXd
ParBicop::hfunc2_deriv2_raw(const Eigen::MatrixXd& u,
                            const Eigen::MatrixXd& parameters,
                            const std::string& deriv)
{
  auto comps = tools_deriv::parse_components(deriv);
  auto first = tools_deriv::comp_to_string(comps[0]);
  auto f = [this, first](const Eigen::MatrixXd& uu, const Eigen::MatrixXd& pp) {
    return hfunc2_deriv_raw(uu, pp, first);
  };
  return fd_deriv(f, u, parameters, comps[1]);
}
//! @}
}
