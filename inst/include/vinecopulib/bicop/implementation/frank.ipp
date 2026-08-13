// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/misc/tools_eigen.hpp>

namespace vinecopulib {

inline FrankBicop::FrankBicop()
{
  family_ = BicopFamily::frank;
  parameters_ = Eigen::VectorXd(1);
  parameters_lower_bounds_ = Eigen::VectorXd(1);
  parameters_upper_bounds_ = Eigen::VectorXd(1);
  parameters_ << 0;
  parameters_lower_bounds_ << -35;
  parameters_upper_bounds_ << 35;
}

inline double
FrankBicop::generator(const double& u,
                      const Eigen::Ref<const Eigen::VectorXd>& parameters)
{
  double theta = parameters(0);
  return -std::log(std::expm1(-theta * u) / std::expm1(-theta));
}

inline double
FrankBicop::generator_inv(const double& u,
                          const Eigen::Ref<const Eigen::VectorXd>& parameters)
{
  double theta = parameters(0);
  return -std::log1p(std::expm1(-theta) * std::exp(-u)) / theta;
}

inline double
FrankBicop::generator_derivative(
  const double& u,
  const Eigen::Ref<const Eigen::VectorXd>& parameters)
{
  double theta = parameters(0);
  return -theta / std::expm1(theta * u);
}

inline Eigen::VectorXd
FrankBicop::pdf_raw(const Eigen::MatrixXd& u, const Eigen::MatrixXd& parameters)
{
  auto f = [](const double& u1,
              const double& u2,
              const Eigen::Ref<const Eigen::VectorXd>& par) {
    double theta = par(0);
    return (theta * std::expm1(theta) *
            std::exp(theta * u2 + theta * u1 + theta)) /
           std::pow(std::exp(theta * u2 + theta * u1) -
                      std::exp(theta * u2 + theta) -
                      std::exp(theta * u1 + theta) + std::exp(theta),
                    2.0);
  };
  return tools_eigen::binaryExpr_or_nan(u, parameters, f);
}

inline Eigen::MatrixXd
FrankBicop::tau_to_parameters(const double& tau)
{
  Eigen::VectorXd tau0 = Eigen::VectorXd::Constant(1, tau);
  auto f = [&](const Eigen::VectorXd& par) {
    return Eigen::VectorXd::Constant(1, parameters_to_tau(par));
  };
  return tools_eigen::invert_f(tau0,
                               f,
                               parameters_lower_bounds_(0) + 1e-6,
                               parameters_upper_bounds_(0) - 1e-5);
}

inline double
FrankBicop::parameters_to_tau(const Eigen::MatrixXd& parameters)
{
  double par = std::fabs(parameters(0));
  if (par < 1e-5) {
    return 0.0;
  }
  double tau = 1 - 4 / par + (4 / par) * debye1(par) / par;
  return (parameters(0) >= 0) ? tau : -tau;
}

inline Eigen::VectorXd
FrankBicop::get_start_parameters(const double tau)
{
  Eigen::VectorXd par = tau_to_parameters(tau);
  par = par.cwiseMax(parameters_lower_bounds_);
  par = par.cwiseMin(parameters_upper_bounds_);
  return par;
}

inline Eigen::MatrixXd
FrankBicop::parameters_to_taildep(const Eigen::MatrixXd&)
{
  return Eigen::MatrixXd::Zero(2, 2);
}

inline Eigen::VectorXd
FrankBicop::pdf_deriv_raw(const Eigen::MatrixXd& u,
                          const Eigen::MatrixXd& parameters,
                          const std::string& deriv)
{
  // exchangeability: route "u2"-flavored selectors through a swap
  if (tools_deriv::is_u2_only(deriv)) {
    return pdf_deriv_raw(
      tools_eigen::swap_cols(u), parameters, tools_deriv::swap_args(deriv));
  }

  if (deriv == "par1") {
    // ported from VineCopula deriv.c diffPDF (family 5 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      double theta = par(0);
      // theta -> 0 is a 0/0 form; clamp |theta| >= 1e-5 preserving sign,
      // consistent with parameters_to_tau()
      theta = (theta >= 0) ? std::max(theta, 1e-5) : std::min(theta, -1e-5);
      double t2 = std::exp(theta);
      double t3 = t2 - 1.0;
      double t4 = theta * u2;
      double t5 = theta * u1;
      double t7 = std::exp(t4 + t5 + theta);
      double t10 = std::exp(t4 + t5);
      double t12 = std::exp(t4 + theta);
      double t14 = std::exp(t5 + theta);
      double t15 = t10 - t12 - t14 + t2;
      double t16 = t15 * t15;
      double t17 = 1.0 / t16;
      double t21 = theta * t3;
      return t3 * t7 * t17 + theta * t2 * t7 * t17 +
             t21 * (u2 + u1 + 1.0) * t7 * t17 -
             2.0 * t21 * t7 / t15 / t16 *
               ((u2 + u1) * t10 - (u2 + 1.0) * t12 - (u1 + 1.0) * t14 + t2);
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  } else if (deriv == "u1") {
    // ported from VineCopula deriv.c diffPDF_u (family 5 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      double theta = par(0);
      // theta -> 0 is a 0/0 form; clamp |theta| >= 1e-5 preserving sign,
      // consistent with parameters_to_tau()
      theta = (theta >= 0) ? std::max(theta, 1e-5) : std::min(theta, -1e-5);
      double t1 = theta * theta;
      double t2 = std::exp(theta);
      double t3 = t2 - 1.0;
      double t5 = theta * u2;
      double t6 = theta * u1;
      double t8 = std::exp(t5 + t6 + theta);
      double t10 = std::exp(t5 + t6);
      double t12 = std::exp(t5 + theta);
      double t14 = std::exp(t6 + theta);
      double t15 = t10 - t12 - t14 + t2;
      double t16 = t15 * t15;
      return t1 * t3 * t8 / t16 -
             2.0 * theta * t3 * t8 / t16 / t15 * (theta * t10 - theta * t14);
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  }
  throw std::runtime_error("unexpected derivative selector: " + deriv);
}

inline Eigen::VectorXd
FrankBicop::pdf_deriv2_raw(const Eigen::MatrixXd& u,
                           const Eigen::MatrixXd& parameters,
                           const std::string& deriv)
{
  // exchangeability: route "u2"-flavored selectors through a swap
  if (tools_deriv::is_u2_only(deriv)) {
    return pdf_deriv2_raw(
      tools_eigen::swap_cols(u), parameters, tools_deriv::swap_args(deriv));
  }

  if (deriv == "par1par1") {
    // ported from VineCopula deriv2.c diff2PDF (family 5 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      double theta = par(0);
      // theta -> 0 is a 0/0 form; clamp |theta| >= 1e-5 preserving sign,
      // consistent with parameters_to_tau()
      theta = (theta >= 0) ? std::max(theta, 1e-5) : std::min(theta, -1e-5);
      double t1 = std::exp(theta);
      double t2 = theta * u2;
      double t3 = theta * u1;
      double t5 = std::exp(t2 + t3 + theta);
      double t8 = std::exp(t2 + t3);
      double t10 = std::exp(t2 + theta);
      double t12 = std::exp(t3 + theta);
      double t13 = t8 - t10 - t12 + t1;
      double t14 = t13 * t13;
      double t15 = 1.0 / t14;
      double t18 = t1 - 1.0;
      double t19 = u2 + u1 + 1.0;
      double t21 = t5 * t15;
      double t26 = 1.0 / t14 / t13;
      double t27 = u2 + u1;
      double t29 = u2 + 1.0;
      double t31 = u1 + 1.0;
      double t33 = t27 * t8 - t29 * t10 - t31 * t12 + t1;
      double t37 = theta * t1;
      double t43 = t5 * t26;
      double t44 = t43 * t33;
      double t47 = theta * t18;
      double t48 = t19 * t19;
      double t11 = t14 * t14;
      double t9 = t33 * t33;
      double t7 = t27 * t27;
      double t6 = t29 * t29;
      double t4 = t31 * t31;
      return 2.0 * t1 * t5 * t15 + 2.0 * t18 * t19 * t21 -
             4.0 * t18 * t5 * t26 * t33 + t37 * t21 +
             2.0 * t37 * t19 * t5 * t15 - 4.0 * t37 * t44 +
             t47 * t48 * t5 * t15 - 4.0 * t47 * t19 * t44 +
             6.0 * t47 * t5 / t11 * t9 -
             2.0 * t47 * t43 * (t7 * t8 - t6 * t10 - t4 * t12 + t1);
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  } else if (deriv == "par1u1") {
    // ported from VineCopula deriv2.c diff2PDF_par_u (family 5 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      double theta = par(0);
      // theta -> 0 is a 0/0 form; clamp |theta| >= 1e-5 preserving sign,
      // consistent with parameters_to_tau()
      theta = (theta >= 0) ? std::max(theta, 1e-5) : std::min(theta, -1e-5);
      double t1 = std::exp(theta);
      double t2 = t1 - 1.0;
      double t3 = theta * t2;
      double t4 = theta * u2;
      double t5 = theta * u1;
      double t7 = std::exp(t4 + t5 + theta);
      double t9 = std::exp(t4 + t5);
      double t11 = std::exp(t4 + theta);
      double t13 = std::exp(t5 + theta);
      double t14 = t9 - t11 - t13 + t1;
      double t15 = t14 * t14;
      double t16 = 1.0 / t15;
      double t17 = t7 * t16;
      double t22 = 1.0 / t15 / t14;
      double t25 = theta * t9 - theta * t13;
      double t29 = theta * theta;
      double t33 = t7 * t22;
      double t34 = t33 * t25;
      double t37 = t29 * t2;
      double t38 = u2 + u1 + 1.0;
      double t46 = u2 + u1;
      double t49 = u1 + 1.0;
      double t51 = t46 * t9 - (u2 + 1.0) * t11 - t49 * t13 + t1;
      double t56 = t15 * t15;
      return 2.0 * t3 * t17 - 2.0 * t2 * t7 * t22 * t25 + t29 * t1 * t17 -
             2.0 * theta * t1 * t34 + t37 * t38 * t7 * t16 -
             2.0 * t3 * t38 * t34 - 2.0 * t37 * t33 * t51 +
             6.0 * t3 * t7 / t56 * t51 * t25 -
             2.0 * t3 * t33 * (t9 + t46 * theta * t9 - t13 - t49 * theta * t13);
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  } else if (deriv == "u1u1") {
    // ported from VineCopula deriv2.c diff2PDF_u (family 5 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      double theta = par(0);
      // theta -> 0 is a 0/0 form; clamp |theta| >= 1e-5 preserving sign,
      // consistent with parameters_to_tau()
      theta = (theta >= 0) ? std::max(theta, 1e-5) : std::min(theta, -1e-5);
      double t1 = theta * theta;
      double t3 = std::exp(theta);
      double t4 = t3 - 1.0;
      double t6 = theta * u2;
      double t7 = theta * u1;
      double t9 = std::exp(t6 + t7 + theta);
      double t11 = std::exp(t6 + t7);
      double t13 = std::exp(t6 + theta);
      double t15 = std::exp(t7 + theta);
      double t16 = t11 - t13 - t15 + t3;
      double t17 = t16 * t16;
      double t24 = t9 / t17 / t16;
      double t27 = theta * t11 - theta * t15;
      double t31 = theta * t4;
      double t32 = t17 * t17;
      double t35 = t27 * t27;
      return t1 * theta * t4 * t9 / t17 - 4.0 * t1 * t4 * t24 * t27 +
             6.0 * t31 * t9 / t32 * t35 -
             2.0 * t31 * t24 * (t1 * t11 - t1 * t15);
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  } else if (deriv == "u1u2") {
    // ported from VineCopula deriv2.c diff2PDF_u_v (family 5 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      double theta = par(0);
      // theta -> 0 is a 0/0 form; clamp |theta| >= 1e-5 preserving sign,
      // consistent with parameters_to_tau()
      theta = (theta >= 0) ? std::max(theta, 1e-5) : std::min(theta, -1e-5);
      double t1 = theta * theta;
      double t3 = std::exp(theta);
      double t4 = t3 - 1.0;
      double t5 = t1 * theta * t4;
      double t6 = theta * u2;
      double t7 = theta * u1;
      double t9 = std::exp(t6 + t7 + theta);
      double t11 = std::exp(t6 + t7);
      double t13 = std::exp(t6 + theta);
      double t15 = std::exp(t7 + theta);
      double t16 = t11 - t13 - t15 + t3;
      double t17 = t16 * t16;
      double t21 = t1 * t4;
      double t24 = t9 / t17 / t16;
      double t25 = theta * t11;
      double t27 = t25 - theta * t13;
      double t32 = t25 - theta * t15;
      double t38 = t17 * t17;
      return t5 * t9 / t17 - 2.0 * t21 * t24 * t27 - 2.0 * t21 * t24 * t32 +
             6.0 * theta * t4 * t9 / t38 * t32 * t27 - 2.0 * t5 * t24 * t11;
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  }
  throw std::runtime_error("unexpected derivative selector: " + deriv);
}

inline Eigen::VectorXd
FrankBicop::hfunc1_deriv_raw(const Eigen::MatrixXd& u,
                             const Eigen::MatrixXd& parameters,
                             const std::string& deriv)
{
  // the VineCopula kernels differentiate h(u|v) = dC/dv; our
  // hfunc1(u1, u2) = dC/du1, so their u := our u2 and their v := our u1
  if (deriv == "par1") {
    // ported from VineCopula hfuncderiv.c diffhfunc (family 5 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      double theta = par(0);
      // theta -> 0 is a 0/0 form; clamp |theta| >= 1e-5 preserving sign,
      // consistent with parameters_to_tau()
      theta = (theta >= 0) ? std::max(theta, 1e-5) : std::min(theta, -1e-5);
      double t1 = std::exp(theta);
      double t2 = theta * u2;
      double t3 = std::exp(t2);
      double t5 = t1 * (t3 - 1.0);
      double t6 = theta * u1;
      double t8 = std::exp(t6 + t2);
      double t9 = std::exp(t6 + theta);
      double t10 = std::exp(t2 + theta);
      double t11 = t8 - t9 - t10 + t1;
      double t14 = 1.0 / t11;
      double t18 = t11 * t11;
      return -t5 * t14 - t1 * u2 * t3 * t14 +
             t5 / t18 *
               ((u1 + u2) * t8 - (u1 + 1.0) * t9 - (u2 + 1.0) * t10 + t1);
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  } else if (deriv == "u1") {
    // ported from VineCopula hfuncderiv.c diffhfunc_v (family 5 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      double theta = par(0);
      // theta -> 0 is a 0/0 form; clamp |theta| >= 1e-5 preserving sign,
      // consistent with parameters_to_tau()
      theta = (theta >= 0) ? std::max(theta, 1e-5) : std::min(theta, -1e-5);
      double t1 = std::exp(theta);
      double t2 = theta * u2;
      double t3 = std::exp(t2);
      double t6 = theta * u1;
      double t8 = std::exp(t6 + t2);
      double t10 = std::exp(t6 + theta);
      double t12 = std::exp(t2 + theta);
      double t13 = std::pow(t8 - t10 - t12 + t1, 2.0);
      return t1 * (t3 - 1.0) / t13 * (theta * t8 - theta * t10);
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  }
  throw std::runtime_error("unexpected derivative selector: " + deriv);
}

inline Eigen::VectorXd
FrankBicop::hfunc1_deriv2_raw(const Eigen::MatrixXd& u,
                              const Eigen::MatrixXd& parameters,
                              const std::string& deriv)
{
  // same argument convention as hfunc1_deriv_raw: their u := our u2,
  // their v := our u1
  if (deriv == "par1par1") {
    // ported from VineCopula hfuncderiv2.c diff2hfunc (family 5 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      double theta = par(0);
      // theta -> 0 is a 0/0 form; clamp |theta| >= 1e-5 preserving sign,
      // consistent with parameters_to_tau()
      theta = (theta >= 0) ? std::max(theta, 1e-5) : std::min(theta, -1e-5);
      double t1 = std::exp(theta);
      double t2 = theta * u2;
      double t3 = std::exp(t2);
      double t5 = t1 * (t3 - 1.0);
      double t6 = theta * u1;
      double t8 = std::exp(t6 + t2);
      double t10 = std::exp(t6 + theta);
      double t12 = std::exp(t2 + theta);
      double t13 = t8 - t10 - t12 + t1;
      double t14 = 1.0 / t13;
      double t16 = t1 * u2;
      double t18 = t3 * t14;
      double t20 = t13 * t13;
      double t21 = 1.0 / t20;
      double t23 = u1 + u2;
      double t25 = u1 + 1.0;
      double t26 = u2 + 1.0;
      double t28 = t23 * t8 - t25 * t10 - t26 * t12 + t1;
      double t32 = u2 * u2;
      double t42 = t28 * t28;
      double t44 = t23 * t23;
      double t47 = t25 * t25;
      double t49 = t26 * t26;
      return -t5 * t14 - 2.0 * t16 * t18 + 2.0 * t5 * t21 * t28 -
             t1 * t32 * t18 + 2.0 * t16 * t3 * t21 * t28 -
             2.0 * t5 / t20 / t13 * t42 +
             t5 * t21 * (t44 * t8 - t47 * t10 - t49 * t12 + t1);
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  } else if (deriv == "par1u1") {
    // ported from VineCopula hfuncderiv2.c diff2hfunc_par_v (family 5
    // branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      double theta = par(0);
      // theta -> 0 is a 0/0 form; clamp |theta| >= 1e-5 preserving sign,
      // consistent with parameters_to_tau()
      theta = (theta >= 0) ? std::max(theta, 1e-5) : std::min(theta, -1e-5);
      double t1 = std::exp(theta);
      double t2 = theta * u2;
      double t3 = std::exp(t2);
      double t5 = t1 * (t3 - 1.0);
      double t6 = theta * u1;
      double t8 = std::exp(t6 + t2);
      double t10 = std::exp(t6 + theta);
      double t12 = std::exp(t2 + theta);
      double t13 = t8 - t10 - t12 + t1;
      double t14 = t13 * t13;
      double t15 = 1.0 / t14;
      double t18 = theta * t8 - theta * t10;
      double t28 = u1 + u2;
      double t29 = u1 + 1.0;
      return t5 * t15 * t18 + t1 * u2 * t3 * t15 * t18 -
             2.0 * t5 / t14 / t13 *
               (t28 * t8 - t29 * t10 - (u2 + 1.0) * t12 + t1) * t18 +
             t5 * t15 * (t8 + t28 * theta * t8 - t10 - t29 * theta * t10);
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  } else if (deriv == "u1u1") {
    // ported from VineCopula hfuncderiv2.c diff2hfunc_v (family 5 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      double theta = par(0);
      // theta -> 0 is a 0/0 form; clamp |theta| >= 1e-5 preserving sign,
      // consistent with parameters_to_tau()
      theta = (theta >= 0) ? std::max(theta, 1e-5) : std::min(theta, -1e-5);
      double t1 = std::exp(theta);
      double t2 = theta * u2;
      double t3 = std::exp(t2);
      double t5 = t1 * (t3 - 1.0);
      double t6 = theta * u1;
      double t8 = std::exp(t6 + t2);
      double t10 = std::exp(t6 + theta);
      double t12 = std::exp(t2 + theta);
      double t13 = t8 - t10 - t12 + t1;
      double t14 = t13 * t13;
      double t20 = std::pow(theta * t8 - theta * t10, 2.0);
      double t24 = theta * theta;
      return -2.0 * t5 / t14 / t13 * t20 + t5 / t14 * (t24 * t8 - t24 * t10);
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  }
  throw std::runtime_error("unexpected derivative selector: " + deriv);
}

inline Eigen::VectorXd
FrankBicop::logpdf_deriv_raw(const Eigen::MatrixXd& u,
                             const Eigen::MatrixXd& parameters,
                             const std::string& deriv)
{
  // fused single-pass form kept for performance: the base
  // logpdf_deriv*_raw would compose this as (pdf_deriv)/pdf, recomputing the
  // shared pdf temporaries 2-3x; this leaf is on the scores/Hessian hot paths.
  if (deriv == "par1") {
    // ported from VineCopula logderiv.c difflPDF (family 5 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      double theta = par(0);
      // theta -> 0 is a 0/0 form; clamp |theta| >= 1e-5 preserving sign,
      // consistent with parameters_to_tau()
      theta = (theta >= 0) ? std::max(theta, 1e-5) : std::min(theta, -1e-5);
      double t1 = std::exp(-theta);
      double t2 = 1.0 - t1;
      double t3 = u1 + u2;
      double t5 = std::exp(-theta * t3);
      double t8 = std::exp(-theta * u1);
      double t9 = 1.0 - t8;
      double t11 = std::exp(-theta * u2);
      double t12 = 1.0 - t11;
      double t14 = 1.0 - t1 - t9 * t12;
      double t15 = t14 * t14;
      double t16 = 1.0 / t15;
      double t17 = theta * t2;
      return (t2 * t5 * t16 + theta * t1 * t5 * t16 - t17 * t3 * t5 * t16 -
              2.0 * t17 * t5 / t15 / t14 *
                (t1 - u1 * t8 * t12 - t9 * u2 * t11)) /
             theta / t2 / t5 * t15;
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  }
  return AbstractBicop::logpdf_deriv_raw(u, parameters, deriv);
}

inline Eigen::VectorXd
FrankBicop::logpdf_deriv2_raw(const Eigen::MatrixXd& u,
                              const Eigen::MatrixXd& parameters,
                              const std::string& deriv)
{
  // fused single-pass form kept for performance: the base
  // logpdf_deriv*_raw would compose this as (pdf_deriv)/pdf, recomputing the
  // shared pdf temporaries 2-3x; this leaf is on the scores/Hessian hot paths.
  if (deriv == "par1par1") {
    // ported from VineCopula logderiv.c diff2lPDF (family 5 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      double theta = par(0);
      // theta -> 0 is a 0/0 form; clamp |theta| >= 1e-5 preserving sign,
      // consistent with parameters_to_tau()
      theta = (theta >= 0) ? std::max(theta, 1e-5) : std::min(theta, -1e-5);
      double t1 = std::exp(theta);
      double t2 = theta * u2;
      double t3 = theta * u1;
      double t5 = std::exp(t2 + t3 + theta);
      double t8 = std::exp(t2 + t3);
      double t10 = std::exp(t2 + theta);
      double t12 = std::exp(t3 + theta);
      double t13 = t8 - t10 - t12 + t1;
      double t14 = t13 * t13;
      double t15 = 1.0 / t14;
      double t18 = t1 - 1.0;
      double t19 = u2 + u1 + 1.0;
      double t21 = t5 * t15;
      double t24 = t18 * t5;
      double t26 = 1.0 / t14 / t13;
      double t27 = u2 + u1;
      double t29 = u2 + 1.0;
      double t31 = u1 + 1.0;
      double t33 = t27 * t8 - t29 * t10 - t31 * t12 + t1;
      double t37 = theta * t1;
      double t38 = t37 * t21;
      double t40 = t19 * t5 * t15;
      double t43 = t5 * t26;
      double t44 = t43 * t33;
      double t47 = theta * t18;
      double t48 = t19 * t19;
      double t55 = t14 * t14;
      double t58 = t33 * t33;
      double t62 = t27 * t27;
      double t64 = t29 * t29;
      double t66 = t31 * t31;
      double t73 = 1.0 / theta;
      double t75 = 1.0 / t18;
      double t76 = 1.0 / t5;
      double t78 = t75 * t76 * t14;
      double t84 = t24 * t15 + t38 + t47 * t40 - 2.0 * t47 * t44;
      double t85 = theta * theta;
      double t89 = t84 * t73;
      double t90 = t18 * t18;
      double t93 = t76 * t14;
      double t96 = t89 * t75;
      return (2.0 * t1 * t5 * t15 + 2.0 * t18 * t19 * t21 -
              4.0 * t24 * t26 * t33 + t38 + 2.0 * t37 * t40 - 4.0 * t37 * t44 +
              t47 * t48 * t5 * t15 - 4.0 * t47 * t19 * t44 +
              6.0 * t47 * t5 / t55 * t58 -
              2.0 * t47 * t43 * (t62 * t8 - t64 * t10 - t66 * t12 + t1)) *
               t73 * t78 -
             t84 / t85 * t78 - t89 / t90 * t93 * t1 - t96 * t93 * t19 +
             2.0 * t96 * t76 * t13 * t33;
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  }
  return AbstractBicop::logpdf_deriv2_raw(u, parameters, deriv);
}

//! @brief computes the Debye function of order 1.
//! @param x the argument and upper limit of the integral. x>=0.
//! @return the Debye function. Zero if x<=0.
//! @details The Debye function is defined as \f$ \int_0^x t  / [exp(t)-1] dt
//! \f$. Code modified from implementation in
//! https://github.com/openturns/openturns.
inline double
debye1(const double& x)
{
  if (x <= 0.)
    return 0.;
  double m_1_2pi = .159154943091895335768883763373; // 1/(2pi)

  if (x >= 3.) {
    double sum = 1.64493406684822643647241516665e+00;
    static constexpr std::array<short, 14> kLim{ 0, 0, 0, 13, 10, 8, 7,
                                                 6, 5, 5, 4,  4,  4, 3 };
    const short kmax = (x < 14) ? kLim[static_cast<int>(x)] : 3;
    for (short k = 1; k <= kmax; k++) {
      const double xk = x * k;
      double ksum = 1. / xk;
      ksum += ksum / xk;
      sum -= exp(-xk) * ksum * x * x;
    }
    return sum;
  } else {
    static constexpr std::array<double, 71> koeff = {
      0.,
      1.289868133696452872944830333292e+00,
      1.646464674222763830320073930823e-01,
      3.468612396889827942903585958184e-02,
      8.154712395888678757370477017305e-03,
      1.989150255636170674291917800638e-03,
      4.921731066160965972759960954793e-04,
      1.224962701174096585170902102707e-04,
      3.056451881730374346514297527344e-05,
      7.634586529999679712923289243879e-06,
      1.907924067745592226304077366899e-06,
      4.769010054554659800072963735060e-07,
      1.192163781025189592248804158716e-07,
      2.980310965673008246931701326140e-08,
      7.450668049576914109638408036805e-09,
      1.862654864839336365743529470042e-09,
      4.656623667353010984002911951881e-10,
      1.164154417580540177848737197821e-10,
      2.910384378208396847185926449064e-11,
      7.275959094757302380474472711747e-12,
      1.818989568052777856506623677390e-12,
      4.547473691649305030453643155957e-13,
      1.136868397525517121855436593505e-13,
      2.842170965606321353966861428348e-14,
      7.105427382674227346596939068119e-15,
      1.776356842186163180619218277278e-15,
      4.440892101596083967998640188409e-16,
      1.110223024969096248744747318102e-16,
      2.775557561945046552567818981300e-17,
      6.938893904331845249488542992219e-18,
      1.734723476023986745668411013469e-18,
      4.336808689994439570027820336642e-19,
      1.084202172491329082183740080878e-19,
      2.710505431220232916297046799365e-20,
      6.776263578041593636171406200902e-21,
      1.694065894509399669649398521836e-21,
      4.235164736272389463688418879636e-22,
      1.058791184067974064762782460584e-22,
      2.646977960169798160618902050189e-23,
      6.617444900424343177893912768629e-24,
      1.654361225106068880734221123349e-24,
      4.135903062765153408791935838694e-25,
      1.033975765691286264082026643327e-25,
      2.584939414228213340076225223666e-26,
      6.462348535570530772269628236053e-27,
      1.615587133892632406631747637268e-27,
      4.038967834731580698317525293132e-28,
      1.009741958682895139216954234507e-28,
      2.524354896707237808750799932127e-29,
      6.310887241768094478219682436680e-30,
      1.577721810442023614704107565240e-30,
      3.944304526105059031370476640000e-31,
      9.860761315262647572437533499000e-32,
      2.465190328815661892443976898000e-32,
      6.162975822039154730370601500000e-33,
      1.540743955509788682510501190000e-33,
      3.851859888774471706184973900000e-34,
      9.629649721936179265360991000000e-35,
      2.407412430484044816328953000000e-35,
      6.018531076210112040809600000000e-36,
      1.504632769052528010200750000000e-36,
      3.761581922631320025497600000000e-37,
      9.403954806578300063715000000000e-38,
      2.350988701644575015901000000000e-38,
      5.877471754111437539470000000000e-39,
      1.469367938527859384580000000000e-39,
      3.673419846319648458500000000000e-40,
      9.183549615799121117000000000000e-41,
      2.295887403949780249000000000000e-41,
      5.739718509874450320000000000000e-42,
      1.434929627468612270000000000000e-42
    };

    double sum = 0, sumold = 1;
    const double x2pi = x * m_1_2pi;
    for (short k = 1; (k < 70) && (sum != sumold); k++) {
      sumold = sum;
      sum += (2. + koeff[k]) * std::pow(x2pi, 2. * k) / (2 * k + 1.);
      k++;
      sum -= (2. + koeff[k]) * std::pow(x2pi, 2. * k) / (2 * k + 1.);
    }
    return x * (sum + 1. - x / 4.);
  }
}
}
