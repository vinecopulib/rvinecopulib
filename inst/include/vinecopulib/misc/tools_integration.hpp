// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

// replace std::sprintf by void function
// - this is necessary because `sprintf()` is now flagged as potential security
//   risk and deprecated in macOS 13.
// - in our case, `boost/odeint` makes (safe) use of the function, but we don't
//   really need it.
// - we can remove this hack if and when updates odeint (PR open).
#define sprintf _sprintf_do_nothing
namespace std
{
constexpr int _sprintf_do_nothing(char*, const char*, ...) { return 0; }
}
#include <boost/numeric/odeint.hpp>
#undef sprintf

#include <functional>

namespace vinecopulib {

namespace tools_integration {

inline double
integrate_zero_to_one(std::function<double(double)> f)
{
  boost::numeric::odeint::runge_kutta_dopri5<double> stepper;
  double lb = 1e-12;
  double ub = 1.0 - lb;
  double x = 0.0;
  auto ifunc = [f](const double /* x */, double& dxdt, const double t) {
    dxdt = f(t);
  };
  integrate_adaptive(boost::numeric::odeint::make_controlled(lb, lb, stepper),
                     ifunc,
                     x,
                     lb,
                     ub,
                     lb);
  return x;
}
}
}
