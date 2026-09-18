// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <boost/math/quadrature/tanh_sinh.hpp>
#include <utility>

namespace vinecopulib {

namespace tools_integration {

//! @brief Integrates `f` over the unit interval.
//!
//! Suited to integrands that are singular or steep at `0` and `1`, as the
//! Kendall's tau integrals of the BB families are: tanh-sinh clusters its
//! abscissas at the endpoints. Use `integrate_zero_to_one_split()` when the
//! integrand instead peaks in the interior.
template<typename F>
inline double
integrate_zero_to_one(F&& f)
{
  // The bounds stay clear of 0 and 1: the integrands are singular there.
  const double lb = 1e-12;
  boost::math::quadrature::tanh_sinh<double> integrator;
  return integrator.integrate(f, lb, 1.0 - lb);
}

//! @brief Integrates `f` over the unit interval, splitting it at `split`.
//!
//! A narrow interior peak is invisible to a rule whose abscissas cluster at the
//! endpoints. Splitting there makes the peak an endpoint of both subintervals,
//! which is exactly where tanh-sinh puts its resolution. Needed for Tawn's tau
//! integral, whose peak is both narrow and far from `1/2` when the parameters
//! are asymmetric.
template<typename F>
inline double
integrate_zero_to_one_split(F&& f, const double split)
{
  const double lb = 1e-12;
  const double ub = 1.0 - lb;
  if (!(split > lb) || !(split < ub)) {
    return integrate_zero_to_one(std::forward<F>(f));
  }
  boost::math::quadrature::tanh_sinh<double> integrator;
  return integrator.integrate(f, lb, split) +
         integrator.integrate(f, split, ub);
}
}
}
