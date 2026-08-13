// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <string>
#include <vector>

namespace vinecopulib {

//! @brief A bivariate copula family identifier.
//!
//! The list below summarises each family's parameter count, parameter
//! range, available rotations, and tail-dependence behaviour. The exact
//! parameter bounds enforced at fit time are visible via
//! `Bicop::get_parameters_lower_bounds()` /
//! `Bicop::get_parameters_upper_bounds()`. The Kendall's-tau column
//! refers to the closed-form mapping implemented by
//! `Bicop::parameters_to_tau()` /  `Bicop::tau_to_parameters()`.
enum class BicopFamily
{
  //! Independence copula. 0 parameters; rotationless; no tail dependence;
  //! Kendall's tau is 0.
  indep,
  //! Gaussian copula. 1 parameter (rho in (-1, 1)); rotationless; no tail
  //! dependence; Kendall's tau is (2 / pi) * arcsin(rho).
  gaussian,
  //! Student-t copula. 2 parameters (rho in (-1, 1), df > 2); rotationless;
  //! symmetric tail dependence; Kendall's tau is (2 / pi) * arcsin(rho).
  student,
  //! Clayton copula. 1 parameter (theta > 0); rotations 0 / 90 / 180 / 270
  //! degrees; lower-tail dependence; Kendall's tau is theta / (theta + 2).
  clayton,
  //! Gumbel copula (also extreme-value). 1 parameter (theta >= 1);
  //! rotations 0 / 90 / 180 / 270 degrees; upper-tail dependence; Kendall's
  //! tau is 1 - 1 / theta.
  gumbel,
  //! Frank copula. 1 parameter (theta in R \ {0}); rotationless; no tail
  //! dependence; Kendall's tau given by the Debye-function form.
  frank,
  //! Joe copula. 1 parameter (theta >= 1); rotations 0 / 90 / 180 / 270
  //! degrees; upper-tail dependence; Kendall's tau via a series expansion.
  joe,
  //! BB1 copula (two-parameter Archimedean). 2 parameters
  //! (theta > 0, delta >= 1); rotations 0 / 90 / 180 / 270 degrees; both
  //! lower- and upper-tail dependence; Kendall's tau in closed form.
  bb1,
  //! BB6 copula (two-parameter Archimedean). 2 parameters
  //! (theta >= 1, delta >= 1); rotations 0 / 90 / 180 / 270 degrees;
  //! upper-tail dependence; Kendall's tau in closed form.
  bb6,
  //! BB7 copula (two-parameter Archimedean). 2 parameters
  //! (theta >= 1, delta > 0); rotations 0 / 90 / 180 / 270 degrees; both
  //! lower- and upper-tail dependence; Kendall's tau in closed form.
  bb7,
  //! BB8 copula (two-parameter Archimedean). 2 parameters
  //! (theta >= 1, delta in (0, 1]); rotations 0 / 90 / 180 / 270 degrees;
  //! upper-tail dependence; Kendall's tau in closed form.
  bb8,
  //! Tawn copula (extreme-value, asymmetric). 3 parameters; rotations
  //! 0 / 90 / 180 / 270 degrees; (asymmetric) upper-tail dependence;
  //! Kendall's tau via the Pickands dependence function.
  tawn,
  //! Transformation Local Likelihood (TLL) nonparametric estimator. No
  //! finite parametric form: the copula density is fit on a grid in the
  //! inverse-normal-transformed copula space. Data-driven rotation and
  //! tail behaviour; Kendall's tau is rank-based on the fitted density.
  tll
};

std::string
get_family_name(BicopFamily family);

BicopFamily
get_family_enum(std::string family);

//! Convenience definitions of sets of bivariate copula families
namespace bicop_families {

//! All implemented families
const std::vector<BicopFamily> all = {
  BicopFamily::indep,   BicopFamily::gaussian, BicopFamily::student,
  BicopFamily::clayton, BicopFamily::gumbel,   BicopFamily::frank,
  BicopFamily::joe,     BicopFamily::bb1,      BicopFamily::bb6,
  BicopFamily::bb7,     BicopFamily::bb8,      BicopFamily::tawn,
  BicopFamily::tll
};

//! All parametric families
const std::vector<BicopFamily> parametric = {
  BicopFamily::indep,   BicopFamily::gaussian, BicopFamily::student,
  BicopFamily::clayton, BicopFamily::gumbel,   BicopFamily::frank,
  BicopFamily::joe,     BicopFamily::bb1,      BicopFamily::bb6,
  BicopFamily::bb7,     BicopFamily::bb8,      BicopFamily::tawn
};

//! All nonparametric families
const std::vector<BicopFamily> nonparametric = { BicopFamily::indep,
                                                 BicopFamily::tll };

//! All one-parameter families
const std::vector<BicopFamily> one_par = {
  BicopFamily::gaussian, BicopFamily::clayton, BicopFamily::gumbel,
  BicopFamily::frank,    BicopFamily::joe,
};

//! All two-parameter families
const std::vector<BicopFamily> two_par = { BicopFamily::student,
                                           BicopFamily::bb1,
                                           BicopFamily::bb6,
                                           BicopFamily::bb7,
                                           BicopFamily::bb8 };

//! All three-parameter families
const std::vector<BicopFamily> three_par = { BicopFamily::tawn };

//! All elliptical copulas
const std::vector<BicopFamily> elliptical = { BicopFamily::gaussian,
                                              BicopFamily::student };

//! All Archimedean copulas
const std::vector<BicopFamily> archimedean = {
  BicopFamily::clayton, BicopFamily::gumbel, BicopFamily::frank,
  BicopFamily::joe,     BicopFamily::bb1,    BicopFamily::bb6,
  BicopFamily::bb7,     BicopFamily::bb8
};

//! All Extreme-value copulas
const std::vector<BicopFamily> extreme_value = { BicopFamily::tawn,
                                                 BicopFamily::gumbel };

//! All BB copulas
const std::vector<BicopFamily> bb = { BicopFamily::bb1,
                                      BicopFamily::bb6,
                                      BicopFamily::bb7,
                                      BicopFamily::bb8 };

//! @brief All copulas that don't have a rotation.
//!
//! (because they already cover positive and negative dependence)
const std::vector<BicopFamily> rotationless = { BicopFamily::indep,
                                                BicopFamily::gaussian,
                                                BicopFamily::student,
                                                BicopFamily::frank,
                                                BicopFamily::tll };

//! Families with stronger dependence in the lower tail
const std::vector<BicopFamily> lt = { BicopFamily::clayton,
                                      BicopFamily::bb1,
                                      BicopFamily::bb7,
                                      BicopFamily::tawn };

//! Families with stronger dependence in the upper tail
const std::vector<BicopFamily> ut = { BicopFamily::gumbel, BicopFamily::joe,
                                      BicopFamily::bb1,    BicopFamily::bb6,
                                      BicopFamily::bb7,    BicopFamily::bb8,
                                      BicopFamily::tawn };

//! Families for which `method = "itau"` is available in Bicop::fit()
const std::vector<BicopFamily> itau = {
  BicopFamily::indep,   BicopFamily::gaussian, BicopFamily::student,
  BicopFamily::clayton, BicopFamily::gumbel,   BicopFamily::frank,
  BicopFamily::joe
};

//! @brief Families with closed-form derivatives of the density and
//! h-functions.
//!
//! For these families, `Bicop::pdf_deriv()`, `Bicop::pdf_deriv2()`,
//! `Bicop::hfunc1_deriv()`, `Bicop::hfunc2_deriv()` (and their second-order
//! and log-density counterparts) evaluate analytical expressions; other
//! parametric families fall back to central finite differences.
const std::vector<BicopFamily> analytic_derivs = {
  BicopFamily::indep,   BicopFamily::gaussian, BicopFamily::student,
  BicopFamily::clayton, BicopFamily::gumbel,   BicopFamily::frank,
  BicopFamily::joe,     BicopFamily::bb1,      BicopFamily::bb6,
  BicopFamily::bb7,     BicopFamily::bb8,      BicopFamily::tawn
};

} // end of namespace BicopFamilies
} // end of namespace vinecopulib

#include <vinecopulib/bicop/implementation/family.ipp>
