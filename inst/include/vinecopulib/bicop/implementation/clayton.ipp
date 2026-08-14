// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <algorithm>
#include <stdexcept>

#include <vinecopulib/misc/tools_eigen.hpp>

namespace vinecopulib {
inline ClaytonBicop::ClaytonBicop()
{
  family_ = BicopFamily::clayton;
  parameters_ = Eigen::VectorXd(1);
  parameters_lower_bounds_ = Eigen::VectorXd(1);
  parameters_upper_bounds_ = Eigen::VectorXd(1);
  parameters_ << 1e-10;
  parameters_lower_bounds_ << 1e-10;
  parameters_upper_bounds_ << 28;
}

inline double
ClaytonBicop::generator(const double& u,
                        const Eigen::Ref<const Eigen::VectorXd>& parameters)
{
  double theta = parameters(0);
  return (std::pow(u, -theta) - 1) / theta;
}

inline double
ClaytonBicop::generator_inv(const double& u,
                            const Eigen::Ref<const Eigen::VectorXd>& parameters)
{
  double theta = parameters(0);
  return std::pow(1 + theta * u, -1 / theta);
}

inline double
ClaytonBicop::generator_derivative(
  const double& u,
  const Eigen::Ref<const Eigen::VectorXd>& parameters)
{
  return (-1) * std::pow(u, -1 - parameters(0));
}

inline Eigen::VectorXd
ClaytonBicop::pdf_raw(const Eigen::MatrixXd& u,
                      const Eigen::MatrixXd& parameters)
{
  auto f = [](const double& u1,
              const double& u2,
              const Eigen::Ref<const Eigen::VectorXd>& par) {
    double theta = par(0);
    // avoid numerical issues when copula is too close to independence
    if (theta < 1e-10) {
      return 1.0;
    }
    double temp = std::log1p(theta) - (1.0 + theta) * std::log(u1 * u2);
    temp = temp - (2.0 + 1.0 / (theta)) *
                    std::log(std::pow(u1, -theta) + std::pow(u2, -theta) - 1.0);
    return std::exp(temp);
  };
  return tools_eigen::binaryExpr_or_nan(u, parameters, f);
}

inline Eigen::VectorXd
ClaytonBicop::pdf_deriv_raw(const Eigen::MatrixXd& u,
                            const Eigen::MatrixXd& parameters,
                            const std::string& deriv)
{
  // the copula is exchangeable: route "u2"-flavored selectors through a
  // column/argument swap
  if (tools_deriv::is_u2_only(deriv)) {
    return pdf_deriv_raw(
      tools_eigen::swap_cols(u), parameters, tools_deriv::swap_args(deriv));
  }
  if (deriv == "par1") {
    // ported from VineCopula deriv.c diffPDF (family 3 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      // clamp theta: the formula has log / 1/theta^2 0*inf forms as theta -> 0
      double theta = std::max(par(0), 1e-5);
      double t1 = u1 * u2;
      double t2 = -theta - 1.0;
      double t3 = std::pow(t1, 1.0 * t2);
      double t4 = std::pow(u1, -1.0 * theta);
      double t5 = std::pow(u2, -1.0 * theta);
      double t6 = t4 + t5 - 1.0;
      double t7 = -2.0 - 1 / theta;
      double t8 = std::pow(t6, 1.0 * t7);
      double t9 = -t2 * t3;
      double t10 = std::log(t1);
      double t11 = theta * theta;
      double t12 = std::log(t6);
      double t13 = std::log(u1);
      double t14 = std::log(u2);
      return t3 * t8 - t9 * t10 * t8 +
             t9 * t8 * (1 / t11 * t12 + t7 * (-t4 * t13 - t5 * t14) / t6);
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  } else if (deriv == "u1") {
    // ported from VineCopula deriv.c diffPDF_u (family 3 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      // clamp theta: the formula has log / 1/theta^2 0*inf forms as theta -> 0
      double theta = std::max(par(0), 1e-5);
      double t1 = 1.0 + theta;
      double t3 = std::pow(u1 * u2, -1.0 * t1);
      double t4 = t1 * t3;
      double t5 = 1 / u1;
      double t7 = std::pow(u1, -1.0 * theta);
      double t8 = std::pow(u2, -1.0 * theta);
      double t9 = t7 + t8 - 1.0;
      double t11 = -2.0 - 1 / theta;
      double t12 = std::pow(t9, 1.0 * t11);
      return -t4 * t1 * t5 * t12 - t4 * t12 * t11 * t7 * theta * t5 / t9;
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  }
  throw std::runtime_error("unexpected derivative selector: " + deriv);
}

inline Eigen::VectorXd
ClaytonBicop::pdf_deriv2_raw(const Eigen::MatrixXd& u,
                             const Eigen::MatrixXd& parameters,
                             const std::string& deriv)
{
  // the copula is exchangeable: route "u2"-flavored selectors through a
  // column/argument swap
  if (tools_deriv::is_u2_only(deriv)) {
    return pdf_deriv2_raw(
      tools_eigen::swap_cols(u), parameters, tools_deriv::swap_args(deriv));
  }
  if (deriv == "par1par1") {
    // ported from VineCopula deriv2.c diff2PDF (family 3 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      // clamp theta: the formula has log / 1/theta^2 0*inf forms as theta -> 0
      double theta = std::max(par(0), 1e-5);
      double t1 = u1 * u2;
      double t2 = -theta - 1.0;
      double t3 = std::pow(t1, 1.0 * t2);
      double t4 = std::log(t1);
      double t6 = std::pow(u1, -1.0 * theta);
      double t7 = std::pow(u2, -1.0 * theta);
      double t8 = t6 + t7 - 1.0;
      double t10 = -2.0 - 1 / theta;
      double t11 = std::pow(t8, 1.0 * t10);
      double t15 = theta * theta;
      double t16 = 1 / t15;
      double t17 = std::log(t8);
      double t19 = std::log(u1);
      double t21 = std::log(u2);
      double t24 = -t6 * t19 - t7 * t21;
      double t26 = 1 / t8;
      double t27 = t16 * t17 + t10 * t24 * t26;
      double t30 = -t2 * t3;
      double t32 = t4 * t4;
      double t14 = t27 * t27;
      double t13 = t19 * t19;
      double t12 = t21 * t21;
      double t9 = t24 * t24;
      double t5 = t8 * t8;
      return -2.0 * t3 * t4 * t11 + 2.0 * t3 * t11 * t27 + t30 * t32 * t11 -
             2.0 * t30 * t4 * t11 * t27 + t30 * t11 * t14 +
             t30 * t11 *
               (-2.0 / t15 / theta * t17 + 2.0 * t16 * t24 * t26 +
                t10 * (t6 * t13 + t7 * t12) * t26 - t10 * t9 / t5);
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  } else if (deriv == "par1u1") {
    // ported from VineCopula deriv2.c diff2PDF_par_u (family 3 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      // clamp theta: the formula has log / 1/theta^2 0*inf forms as theta -> 0
      double theta = std::max(par(0), 1e-5);
      double t1 = u1 * u2;
      double t2 = -theta - 1.0;
      double t3 = std::pow(t1, 1.0 * t2);
      double t5 = 1 / u1;
      double t6 = std::pow(u1, -1.0 * theta);
      double t7 = std::pow(u2, -1.0 * theta);
      double t8 = t6 + t7 - 1.0;
      double t9 = 1 / theta;
      double t10 = -2.0 - t9;
      double t11 = std::pow(t8, 1.0 * t10);
      double t12 = t5 * t11;
      double t16 = t6 * theta;
      double t17 = 1 / t8;
      double t18 = t5 * t17;
      double t21 = -t3 * t2;
      double t22 = t21 * t2;
      double t23 = std::log(t1);
      double t35 = theta * theta;
      double t37 = std::log(t8);
      double t39 = std::log(u1);
      double t41 = std::log(u2);
      double t44 = t10 * (-t6 * t39 - t7 * t41);
      double t46 = 1 / t35 * t37 + t44 * t17;
      double t62 = t8 * t8;
      return t3 * t2 * t12 - t3 * t11 * t10 * t16 * t18 - t22 * t5 * t23 * t11 -
             t21 * t12 + t21 * t23 * t11 * t10 * t6 * theta * t5 * t17 +
             t22 * t12 * t46 - t21 * t11 * t10 * t16 * t18 * t46 +
             t21 * t11 *
               (-t9 * t6 * t18 + t10 * (t16 * t5 * t39 - t6 * t5) * t17 +
                t44 / t62 * t16 * t5);
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  } else if (deriv == "u1u1") {
    // ported from VineCopula deriv2.c diff2PDF_u (family 3 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      // clamp theta: the formula has log / 1/theta^2 0*inf forms as theta -> 0
      double theta = std::max(par(0), 1e-5);
      double t1 = 1.0 + theta;
      double t3 = std::pow(u1 * u2, -1.0 * t1);
      double t4 = t1 * t3;
      double t5 = t1 * t1;
      double t6 = u1 * u1;
      double t7 = 1 / t6;
      double t9 = std::pow(u1, -1.0 * theta);
      double t10 = std::pow(u2, -1.0 * theta);
      double t11 = t9 + t10 - 1.0;
      double t13 = -2.0 - 1 / theta;
      double t14 = std::pow(t11, 1.0 * t13);
      double t17 = -t1 * t7;
      double t21 = t14 * t13;
      double t22 = t9 * theta;
      double t23 = 1 / t11;
      double t28 = t13 * t13;
      double t31 = t9 * t9;
      double t32 = theta * theta;
      double t34 = t11 * t11;
      double t37 = t31 * t32 * t7 / t34;
      double t39 = t4 * t21;
      double t41 = t7 * t23;
      return t4 * t5 * t7 * t14 - t4 * t17 * t14 -
             2.0 * t4 * t17 * t21 * t22 * t23 + t4 * t14 * t28 * t37 +
             t39 * t9 * t32 * t41 + t39 * t22 * t41 - t39 * t37;
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  } else if (deriv == "u1u2") {
    // ported from VineCopula deriv2.c diff2PDF_u_v (family 3 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      // clamp theta: the formula has log / 1/theta^2 0*inf forms as theta -> 0
      double theta = std::max(par(0), 1e-5);
      double t1 = 1.0 + theta;
      double t3 = std::pow(u1 * u2, -1.0 * t1);
      double t4 = t1 * t3;
      double t5 = t1 * t1;
      double t7 = 1 / u2;
      double t8 = 1 / u1;
      double t10 = std::pow(u1, -1.0 * theta);
      double t11 = std::pow(u2, -1.0 * theta);
      double t12 = t10 + t11 - 1.0;
      double t14 = -2.0 - 1 / theta;
      double t15 = std::pow(t12, 1.0 * t14);
      double t23 = 1 / t12;
      double t35 = t14 * t14;
      double t39 = theta * theta;
      double t41 = t12 * t12;
      double t42 = 1 / t41;
      return t4 * t5 * t7 * t8 * t15 +
             t4 * t1 * t8 * t15 * t14 * t11 * theta * t7 * t23 +
             t4 * t1 * t7 * t15 * t14 * t10 * theta * t8 * t23 +
             t4 * t15 * t35 * t11 * t39 * t7 * t42 * t10 * t8 -
             t4 * t15 * t14 * t10 * t39 * t8 * t42 * t11 * t7;
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  }
  throw std::runtime_error("unexpected derivative selector: " + deriv);
}

inline Eigen::VectorXd
ClaytonBicop::hfunc1_deriv_raw(const Eigen::MatrixXd& u,
                               const Eigen::MatrixXd& parameters,
                               const std::string& deriv)
{
  // VineCopula's kernels differentiate h(u|v) = dC/dv; our hfunc1 conditions
  // on the first argument, so their (u, v) = our (u2, u1)
  if (deriv == "par1") {
    // ported from VineCopula hfuncderiv.c diffhfunc (family 3 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      // clamp theta: the formula has log / 1/theta^2 0*inf forms as theta -> 0
      double theta = std::max(par(0), 1e-5);
      double t1 = std::pow(u1, -1.0 * theta - 1.0);
      double t2 = std::log(u1);
      double t3 = std::pow(u2, -1.0 * theta);
      double t4 = std::pow(u1, -1.0 * theta);
      double t5 = t3 + t4 - 1.0;
      double t6 = -1.0 - 1 / theta;
      double t7 = std::pow(t5, 1.0 * t6);
      double t8 = theta * theta;
      double t9 = std::log(t5);
      double t10 = std::log(u2);
      return -t1 * t2 * t7 +
             t1 * t7 * (1 / t8 * t9 + t6 * (-t3 * t10 - t4 * t2) / t5);
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  } else if (deriv == "u1") {
    // ported from VineCopula hfuncderiv.c diffhfunc_v (family 3 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      // clamp theta: the formula has log / 1/theta^2 0*inf forms as theta -> 0
      double theta = std::max(par(0), 1e-5);
      double t1 = -theta - 1.0;
      double t2 = std::pow(u1, 1.0 * t1);
      double t4 = 1 / u1;
      double t5 = std::pow(u2, -1.0 * theta);
      double t6 = std::pow(u1, -1.0 * theta);
      double t7 = t5 + t6 - 1.0;
      double t9 = -1.0 - 1 / theta;
      double t10 = std::pow(t7, 1.0 * t9);
      return t10 * t4 * t1 * t2 - 1 / t7 * t4 * theta * t6 * t9 * t10 * t2;
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  }
  throw std::runtime_error("unexpected derivative selector: " + deriv);
}

inline Eigen::VectorXd
ClaytonBicop::hfunc1_deriv2_raw(const Eigen::MatrixXd& u,
                                const Eigen::MatrixXd& parameters,
                                const std::string& deriv)
{
  // VineCopula's kernels differentiate h(u|v) = dC/dv; our hfunc1 conditions
  // on the first argument, so their (u, v) = our (u2, u1)
  if (deriv == "par1par1") {
    // ported from VineCopula hfuncderiv2.c diff2hfunc (family 3 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      // clamp theta: the formula has log / 1/theta^2 0*inf forms as theta -> 0
      double theta = std::max(par(0), 1e-5);
      double t2 = std::pow(u1, -1.0 * theta - 1.0);
      double t3 = std::log(u1);
      double t4 = t3 * t3;
      double t6 = std::pow(u2, -1.0 * theta);
      double t7 = std::pow(u1, -1.0 * theta);
      double t8 = t6 + t7 - 1.0;
      double t10 = -1.0 - 1 / theta;
      double t11 = std::pow(t8, 1.0 * t10);
      double t14 = theta * theta;
      double t15 = 1 / t14;
      double t16 = std::log(t8);
      double t18 = std::log(u2);
      double t21 = -t6 * t18 - t7 * t3;
      double t23 = 1 / t8;
      double t25 = t15 * t16 + t10 * t21 * t23;
      double t29 = t2 * t11;
      double t32 = t25 * t25;
      double t12 = t18 * t18;
      double t9 = t21 * t21;
      double t5 = t8 * t8;
      return t2 * t4 * t11 - 2.0 * t2 * t3 * t11 * t25 + t29 * t32 +
             t29 * (-2.0 / t14 / theta * t16 + 2.0 * t15 * t21 * t23 +
                    t10 * (t6 * t12 + t7 * t4) * t23 - t10 * t9 / t5);
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  } else if (deriv == "par1u1") {
    // ported from VineCopula hfuncderiv2.c diff2hfunc_par_v (family 3 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      // clamp theta: the formula has log / 1/theta^2 0*inf forms as theta -> 0
      double theta = std::max(par(0), 1e-5);
      double t1 = -theta - 1.0;
      double t2 = std::pow(u1, 1.0 * t1);
      double t3 = t2 * t1;
      double t4 = 1 / u1;
      double t5 = std::log(u1);
      double t6 = t4 * t5;
      double t7 = std::pow(u2, -1.0 * theta);
      double t8 = std::pow(u1, -1.0 * theta);
      double t9 = t7 + t8 - 1.0;
      double t10 = 1 / theta;
      double t11 = -1.0 - t10;
      double t12 = std::pow(t9, 1.0 * t11);
      double t20 = t8 * theta;
      double t21 = 1 / t9;
      double t22 = t4 * t21;
      double t26 = theta * theta;
      double t28 = std::log(t9);
      double t30 = std::log(u2);
      double t34 = t11 * (-t7 * t30 - t8 * t5);
      double t36 = 1 / t26 * t28 + t34 * t21;
      double t39 = t2 * t12;
      double t53 = t9 * t9;
      return -t3 * t6 * t12 - t2 * t4 * t12 + t2 * t5 * t12 * t11 * t20 * t22 +
             t3 * t4 * t12 * t36 - t39 * t11 * t8 * theta * t4 * t21 * t36 +
             t39 * (-t10 * t8 * t22 + t11 * (t20 * t6 - t8 * t4) * t21 +
                    t34 / t53 * t20 * t4);
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  } else if (deriv == "u1u1") {
    // ported from VineCopula hfuncderiv2.c diff2hfunc_v (family 3 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      // clamp theta: the formula has log / 1/theta^2 0*inf forms as theta -> 0
      double theta = std::max(par(0), 1e-5);
      double t1 = -theta - 1.0;
      double t2 = std::pow(u1, 1.0 * t1);
      double t3 = t1 * t1;
      double t5 = u1 * u1;
      double t6 = 1 / t5;
      double t7 = std::pow(u2, -1.0 * theta);
      double t8 = std::pow(u1, -1.0 * theta);
      double t9 = t7 + t8 - 1.0;
      double t11 = -1.0 - 1 / theta;
      double t12 = std::pow(t9, 1.0 * t11);
      double t13 = t6 * t12;
      double t16 = t2 * t1 * t13;
      double t18 = 1 / t9;
      double t23 = t2 * t12;
      double t24 = t11 * t11;
      double t26 = t8 * t8;
      double t27 = theta * theta;
      double t29 = t9 * t9;
      double t32 = t26 * t27 * t6 / t29;
      double t34 = t23 * t11;
      double t36 = t6 * t18;
      return t2 * t3 * t13 - t16 - 2.0 * t16 * t11 * t8 * theta * t18 +
             t23 * t24 * t32 + t34 * t8 * t27 * t36 + t34 * t8 * theta * t36 -
             t34 * t32;
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  }
  throw std::runtime_error("unexpected derivative selector: " + deriv);
}

inline Eigen::VectorXd
ClaytonBicop::logpdf_deriv_raw(const Eigen::MatrixXd& u,
                               const Eigen::MatrixXd& parameters,
                               const std::string& deriv)
{
  // fused single-pass form kept for performance: the base
  // logpdf_deriv*_raw would compose this as (pdf_deriv)/pdf, recomputing the
  // shared pdf temporaries 2-3x; this leaf is on the scores/Hessian hot paths.
  if (deriv == "par1") {
    // ported from VineCopula logderiv.c difflPDF (family 3 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      // clamp theta: the formula has log / 1/theta^2 0*inf forms as theta -> 0
      double theta = std::max(par(0), 1e-5);
      double t4 = std::log(u1 * u2);
      double t5 = std::pow(u1, -1.0 * theta);
      double t6 = std::pow(u2, -1.0 * theta);
      double t7 = t5 + t6 - 1.0;
      double t8 = std::log(t7);
      double t9 = theta * theta;
      double t14 = std::log(u1);
      double t16 = std::log(u2);
      return 1 / (1.0 + theta) - t4 + t8 / t9 +
             (1 / theta + 2.0) * (t5 * t14 + t6 * t16) / t7;
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  }
  return AbstractBicop::logpdf_deriv_raw(u, parameters, deriv);
}

inline Eigen::VectorXd
ClaytonBicop::logpdf_deriv2_raw(const Eigen::MatrixXd& u,
                                const Eigen::MatrixXd& parameters,
                                const std::string& deriv)
{
  // fused single-pass form kept for performance: the base
  // logpdf_deriv*_raw would compose this as (pdf_deriv)/pdf, recomputing the
  // shared pdf temporaries 2-3x; this leaf is on the scores/Hessian hot paths.
  if (deriv == "par1par1") {
    // ported from VineCopula logderiv.c diff2lPDF (family 3 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      // clamp theta: the formula has log / 1/theta^2 0*inf forms as theta -> 0
      double theta = std::max(par(0), 1e-5);
      double t1 = u1 * u2;
      double t2 = -theta - 1.0;
      double t3 = std::pow(t1, 1.0 * t2);
      double t4 = std::log(t1);
      double t6 = std::pow(u1, -1.0 * theta);
      double t7 = std::pow(u2, -1.0 * theta);
      double t8 = t6 + t7 - 1.0;
      double t10 = -2.0 - 1 / theta;
      double t11 = std::pow(t8, 1.0 * t10);
      double t14 = t3 * t11;
      double t15 = theta * theta;
      double t16 = 1 / t15;
      double t17 = std::log(t8);
      double t19 = std::log(u1);
      double t21 = std::log(u2);
      double t23 = -t6 * t19 - t7 * t21;
      double t25 = 1 / t8;
      double t27 = t16 * t17 + t10 * t23 * t25;
      double t30 = -t2 * t3;
      double t31 = t4 * t4;
      double t34 = t4 * t11;
      double t38 = t27 * t27;
      double t48 = t19 * t19;
      double t50 = t21 * t21;
      double t55 = t23 * t23;
      double t57 = t8 * t8;
      double t64 = -1 / t2;
      double t68 = 1 / t3 / t11;
      double t73 = t14 - t30 * t34 + t30 * t11 * t27;
      double t74 = t2 * t2;
      double t78 = t73 * t64;
      return (-2.0 * t3 * t4 * t11 + 2.0 * t14 * t27 + t30 * t31 * t11 -
              2.0 * t30 * t34 * t27 + t30 * t11 * t38 +
              t30 * t11 *
                (-2.0 / t15 / theta * t17 + 2.0 * t16 * t23 * t25 +
                 t10 * (t6 * t48 + t7 * t50) * t25 - t10 * t55 / t57)) *
               t64 * t68 -
             t73 / t74 * t68 + t78 * t68 * t4 - t78 * t68 * t27;
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  }
  return AbstractBicop::logpdf_deriv2_raw(u, parameters, deriv);
}

inline Eigen::VectorXd
ClaytonBicop::hinv1_raw(const Eigen::MatrixXd& u,
                        const Eigen::MatrixXd& parameters)
{
  // theta is bounded above by 28 < 75, so the closed form is always used
  const Eigen::Index n = u.rows();
  Eigen::ArrayXd theta =
    tools_eigen::parameter_as_vector(parameters, 0, n).array();
  Eigen::ArrayXd u0 = u.col(0).array();
  Eigen::ArrayXd u1 = u.col(1).array();
  Eigen::ArrayXd hinv = u0.pow(theta + 1.0);
  hinv = u1 * hinv;
  hinv = hinv.pow(-theta / (theta + 1.0));
  Eigen::ArrayXd x = u0.pow(-theta);
  hinv = hinv - x + 1.0;
  hinv = hinv.pow(-1.0 / theta);
  return hinv.matrix();
}

inline Eigen::MatrixXd
ClaytonBicop::tau_to_parameters(const double& tau)
{
  Eigen::VectorXd parameters(1);
  parameters(0) = 2 * std::fabs(tau) / (1 - std::fabs(tau));
  return parameters.cwiseMax(parameters_lower_bounds_)
    .cwiseMin(parameters_upper_bounds_);
}

inline double
ClaytonBicop::parameters_to_tau(const Eigen::MatrixXd& parameters)
{
  return parameters(0) / (2 + std::fabs(parameters(0)));
}

inline Eigen::MatrixXd
ClaytonBicop::parameters_to_taildep(const Eigen::MatrixXd& parameters)
{
  Eigen::MatrixXd taildep = Eigen::MatrixXd::Zero(2, 2);
  taildep(0, 0) = std::pow(2.0, -1.0 / parameters(0)); // lower tail dependence
  return taildep;
}

inline Eigen::VectorXd
ClaytonBicop::get_start_parameters(const double tau)
{
  Eigen::VectorXd par = tau_to_parameters(tau);
  par = par.cwiseMax(parameters_lower_bounds_);
  par = par.cwiseMin(parameters_upper_bounds_);
  return par;
}
}
