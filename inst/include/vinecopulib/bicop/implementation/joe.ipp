// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <boost/math/special_functions/digamma.hpp>
#include <cmath>
#include <stdexcept>
#include <vinecopulib/misc/tools_eigen.hpp>

namespace vinecopulib {
inline JoeBicop::JoeBicop()
{
  family_ = BicopFamily::joe;
  parameters_ = Eigen::VectorXd(1);
  parameters_lower_bounds_ = Eigen::VectorXd(1);
  parameters_upper_bounds_ = Eigen::VectorXd(1);
  parameters_ << 1;
  parameters_lower_bounds_ << 1;
  parameters_upper_bounds_ << 30;
}

inline double
JoeBicop::generator(const double& u,
                    const Eigen::Ref<const Eigen::VectorXd>& parameters)
{
  return (-1) * std::log1p(-std::pow(1 - u, parameters(0)));
}

inline double
JoeBicop::generator_inv(const double& u,
                        const Eigen::Ref<const Eigen::VectorXd>& parameters)
{
  return 1 - std::pow(-std::expm1(-u), 1 / parameters(0));
}

inline double
JoeBicop::generator_derivative(
  const double& u,
  const Eigen::Ref<const Eigen::VectorXd>& parameters)
{
  double theta = parameters(0);
  return (-theta) * std::pow(1 - u, theta - 1) / (1 - std::pow(1 - u, theta));
}

inline Eigen::VectorXd
JoeBicop::pdf_raw(const Eigen::MatrixXd& u, const Eigen::MatrixXd& parameters)
{
  auto f = [](const double& u1,
              const double& u2,
              const Eigen::Ref<const Eigen::VectorXd>& par) {
    double theta = par(0);
    double t1 = std::pow(1 - u1, theta);
    double t2 = std::pow(1 - u2, theta);
    return std::pow(t1 + t2 - t1 * t2, 1 / theta - 2) *
           std::pow(1 - u1, theta - 1) * std::pow(1 - u2, theta - 1) *
           (theta - 1 + t1 + t2 - t1 * t2);
  };
  return tools_eigen::binaryExpr_or_nan(u, parameters, f);
}

// inverse h-function
inline Eigen::VectorXd
JoeBicop::hinv1_raw(const Eigen::MatrixXd& u, const Eigen::MatrixXd& parameters)
{
  auto qcondjoe_func =
    [](const double& u1,
       const double& u2,
       const Eigen::Ref<const Eigen::VectorXd>& par) -> double {
    return qcondjoe(u2, u1, par(0));
  };
  return tools_eigen::binaryExpr_or_nan(u, parameters, qcondjoe_func);
}

inline Eigen::VectorXd
JoeBicop::pdf_deriv_raw(const Eigen::MatrixXd& u,
                        const Eigen::MatrixXd& parameters,
                        const std::string& deriv)
{
  if (tools_deriv::is_u2_only(deriv)) {
    // exchangeability: c(u1, u2) = c(u2, u1)
    return pdf_deriv_raw(
      tools_eigen::swap_cols(u), parameters, tools_deriv::swap_args(deriv));
  }
  if (deriv == "par1") {
    // ported from VineCopula deriv.c diffPDF (family 6 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      double theta = par(0);
      double t1 = 1.0 - u1;
      double t2 = std::pow(t1, 1.0 * theta);
      double t3 = 1.0 - u2;
      double t4 = std::pow(t3, 1.0 * theta);
      double t5 = t2 * t4;
      double t6 = t2 + t4 - t5;
      double t8 = 1 / theta - 2.0;
      double t9 = std::pow(t6, 1.0 * t8);
      double t10 = theta * theta;
      double t11 = std::log(t6);
      double t12 = std::log(t1);
      double t13 = t2 * t12;
      double t14 = std::log(t3);
      double t15 = t4 * t14;
      double t16 = t13 * t4;
      double t19 = t5 * t14;
      double t21 = theta - 1.0;
      double t27 = std::pow(t1, 1.0 * t21);
      double t28 = std::pow(t3, 1.0 * t21);
      double t30 = theta - 1.0 + t2 + t4 - t5;
      double t33 = t9 * t27;
      return t9 * (-1 / t10 * t11 + t8 * (t13 + t15 - t16 - t19) / t6) * t27 *
               t28 * t30 +
             t33 * t12 * t28 * t30 + t33 * t28 * t14 * t30 +
             t33 * t28 * (1.0 + t13 + t15 - t16 - t19);
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  } else if (deriv == "u1") {
    // ported from VineCopula deriv.c diffPDF_u (family 6 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      double theta = par(0);
      double t1 = 1.0 - u1;
      double t2 = std::pow(t1, 1.0 * theta);
      double t3 = 1.0 - u2;
      double t4 = std::pow(t3, 1.0 * theta);
      double t5 = t2 * t4;
      double t6 = t2 + t4 - t5;
      double t8 = 1 / theta - 2.0;
      double t9 = std::pow(t6, 1.0 * t8);
      double t11 = t2 * theta;
      double t12 = 1 / t1;
      double t16 = -t11 * t12 + t11 * t12 * t4;
      double t19 = theta - 1.0;
      double t20 = std::pow(t1, 1.0 * t19);
      double t22 = std::pow(t3, 1.0 * t19);
      double t23 = theta - 1.0 + t2 + t4 - t5;
      double t27 = t9 * t20;
      return t9 * t8 * t16 / t6 * t20 * t22 * t23 -
             t27 * t19 * t12 * t22 * t23 + t27 * t22 * t16;
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  }
  throw std::runtime_error("unexpected derivative selector: " + deriv);
}

inline Eigen::VectorXd
JoeBicop::pdf_deriv2_raw(const Eigen::MatrixXd& u,
                         const Eigen::MatrixXd& parameters,
                         const std::string& deriv)
{
  if (tools_deriv::is_u2_only(deriv)) {
    // exchangeability: c(u1, u2) = c(u2, u1)
    return pdf_deriv2_raw(
      tools_eigen::swap_cols(u), parameters, tools_deriv::swap_args(deriv));
  }
  if (deriv == "par1par1") {
    // ported from VineCopula deriv2.c diff2PDF (family 6 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      double theta = par(0);
      double t1 = 1.0 - u1;
      double t2 = std::pow(t1, 1.0 * theta);
      double t3 = 1.0 - u2;
      double t4 = std::pow(t3, 1.0 * theta);
      double t5 = t2 * t4;
      double t6 = t2 + t4 - t5;
      double t8 = 1 / theta - 2.0;
      double t9 = std::pow(t6, 1.0 * t8);
      double t10 = theta * theta;
      double t11 = 1 / t10;
      double t12 = std::log(t6);
      double t14 = std::log(t1);
      double t15 = t2 * t14;
      double t16 = std::log(t3);
      double t17 = t4 * t16;
      double t18 = t15 * t4;
      double t19 = t5 * t16;
      double t20 = t15 + t17 - t18 - t19;
      double t22 = 1 / t6;
      double t24 = -t11 * t12 + t8 * t20 * t22;
      double t25 = t24 * t24;
      double t27 = theta - 1.0;
      double t28 = std::pow(t1, 1.0 * t27);
      double t29 = std::pow(t3, 1.0 * t27);
      double t30 = t28 * t29;
      double t31 = theta - 1.0 + t2 + t4 - t5;
      double t32 = t30 * t31;
      double t41 = t14 * t14;
      double t42 = t2 * t41;
      double t43 = t16 * t16;
      double t49 = t42 + t4 * t43 - t42 * t4 - 2.0 * t15 * t17 - t5 * t43;
      double t51 = t20 * t20;
      double t53 = t6 * t6;
      double t60 = t9 * t24;
      double t61 = t60 * t28;
      double t62 = t14 * t29;
      double t66 = t29 * t16;
      double t67 = t66 * t31;
      double t70 = 1.0 + t15 + t17 - t18 - t19;
      double t74 = t9 * t28;
      return t9 * t25 * t32 +
             t9 *
               (2.0 / t10 / theta * t12 - 2.0 * t11 * t20 * t22 +
                t8 * t49 * t22 - t8 * t51 / t53) *
               t32 +
             2.0 * t61 * t62 * t31 + 2.0 * t61 * t67 + 2.0 * t60 * t30 * t70 +
             t74 * t41 * t29 * t31 + 2.0 * t74 * t14 * t67 +
             2.0 * t74 * t62 * t70 + t74 * t29 * t43 * t31 +
             2.0 * t74 * t66 * t70 + t74 * t29 * t49;
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  } else if (deriv == "par1u1") {
    // ported from VineCopula deriv2.c diff2PDF_par_u (family 6 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      double theta = par(0);
      double t1 = 1.0 - u1;
      double t2 = std::pow(t1, 1.0 * theta);
      double t3 = 1.0 - u2;
      double t4 = std::pow(t3, 1.0 * theta);
      double t5 = t2 * t4;
      double t6 = t2 + t4 - t5;
      double t8 = 1 / theta - 2.0;
      double t9 = std::pow(t6, 1.0 * t8);
      double t10 = t9 * t8;
      double t11 = t2 * theta;
      double t12 = 1 / t1;
      double t14 = t12 * t4;
      double t16 = -t11 * t12 + t11 * t14;
      double t17 = 1 / t6;
      double t19 = t10 * t16 * t17;
      double t20 = theta * theta;
      double t21 = 1 / t20;
      double t22 = std::log(t6);
      double t24 = std::log(t1);
      double t25 = t2 * t24;
      double t26 = std::log(t3);
      double t27 = t4 * t26;
      double t28 = t25 * t4;
      double t29 = t5 * t26;
      double t31 = t8 * (t25 + t27 - t28 - t29);
      double t33 = -t21 * t22 + t31 * t17;
      double t34 = theta - 1.0;
      double t35 = std::pow(t1, 1.0 * t34);
      double t37 = std::pow(t3, 1.0 * t34);
      double t38 = theta - 1.0 + t2 + t4 - t5;
      double t39 = t37 * t38;
      double t44 = t12 * t24;
      double t46 = t2 * t12;
      double t52 =
        -t11 * t44 - t46 + t11 * t44 * t4 + t46 * t4 + t11 * t14 * t26;
      double t55 = t6 * t6;
      double t61 = t35 * t37;
      double t64 = t9 * t33;
      double t71 = t9 * t35;
      double t80 = t71 * t34;
      double t81 = t12 * t37;
      double t87 = t26 * t38;
      double t94 = 1.0 + t25 + t27 - t28 - t29;
      return t19 * t33 * t35 * t39 +
             t9 * (-t21 * t16 * t17 + t8 * t52 * t17 - t31 / t55 * t16) * t61 *
               t38 -
             t64 * t35 * t34 * t12 * t39 + t64 * t61 * t16 +
             t19 * t35 * t24 * t39 - t80 * t44 * t39 - t71 * t81 * t38 +
             t71 * t24 * t37 * t16 + t19 * t61 * t87 - t80 * t81 * t87 +
             t71 * t37 * t26 * t16 + t10 * t16 * t17 * t35 * t37 * t94 -
             t80 * t81 * t94 + t71 * t37 * t52;
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  } else if (deriv == "u1u1") {
    // ported from VineCopula deriv2.c diff2PDF_u (family 6 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      double theta = par(0);
      double t1 = 1.0 - u1;
      double t2 = std::pow(t1, 1.0 * theta);
      double t3 = 1.0 - u2;
      double t4 = std::pow(t3, 1.0 * theta);
      double t5 = t2 * t4;
      double t6 = t2 + t4 - t5;
      double t8 = 1 / theta - 2.0;
      double t9 = std::pow(t6, 1.0 * t8);
      double t10 = t8 * t8;
      double t12 = t2 * theta;
      double t13 = 1 / t1;
      double t17 = -t12 * t13 + t12 * t13 * t4;
      double t18 = t17 * t17;
      double t20 = t6 * t6;
      double t22 = theta - 1.0;
      double t23 = std::pow(t1, 1.0 * t22);
      double t25 = std::pow(t3, 1.0 * t22);
      double t26 = theta - 1.0 + t2 + t4 - t5;
      double t27 = t25 * t26;
      double t28 = 1 / t20 * t23 * t27;
      double t30 = t9 * t8;
      double t31 = theta * theta;
      double t32 = t2 * t31;
      double t33 = t1 * t1;
      double t34 = 1 / t33;
      double t37 = t34 * t4;
      double t40 = t32 * t34 - t12 * t34 - t32 * t37 + t12 * t37;
      double t42 = 1 / t6;
      double t43 = t42 * t23;
      double t46 = t30 * t18;
      double t16 = t13 * t25;
      double t15 = t9 * t23;
      double t14 = t22 * t22;
      double t11 = t34 * t25 * t26;
      double t7 = t15 * t22;
      return t9 * t10 * t18 * t28 + t30 * t40 * t43 * t27 - t46 * t28 -
             2.0 * t30 * t17 * t42 * t23 * t22 * t16 * t26 +
             2.0 * t46 * t43 * t25 + t15 * t14 * t11 - t7 * t11 -
             2.0 * t7 * t16 * t17 + t15 * t25 * t40;
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  } else if (deriv == "u1u2") {
    // ported from VineCopula deriv2.c diff2PDF_u_v (family 6 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      double theta = par(0);
      double t1 = 1.0 - u1;
      double t2 = std::pow(t1, 1.0 * theta);
      double t3 = 1.0 - u2;
      double t4 = std::pow(t3, 1.0 * theta);
      double t5 = t2 * t4;
      double t6 = t2 + t4 - t5;
      double t8 = 1 / theta - 2.0;
      double t9 = std::pow(t6, 1.0 * t8);
      double t10 = t8 * t8;
      double t13 = 1 / t3;
      double t17 = -t4 * theta * t13 + t5 * theta * t13;
      double t18 = t6 * t6;
      double t19 = 1 / t18;
      double t22 = t2 * theta;
      double t23 = 1 / t1;
      double t27 = -t22 * t23 + t22 * t23 * t4;
      double t28 = theta - 1.0;
      double t29 = std::pow(t1, 1.0 * t28);
      double t31 = std::pow(t3, 1.0 * t28);
      double t32 = theta - 1.0 + t2 + t4 - t5;
      double t36 = t9 * t8;
      double t37 = theta * theta;
      double t41 = t4 * t13;
      double t42 = 1 / t6;
      double t45 = t29 * t31;
      double t55 = t28 * t13;
      double t71 = t23 * t31;
      double t72 = t9 * t29;
      double t7 = t28 * t28;
      return t9 * t10 * t17 * t19 * t27 * t29 * t31 * t32 -
             t36 * t2 * t37 * t23 * t41 * t42 * t45 * t32 -
             t36 * t27 * t19 * t45 * t32 * t17 -
             t36 * t27 * t42 * t45 * t55 * t32 +
             2.0 * t36 * t27 * t42 * t29 * t31 * t17 -
             t36 * t17 * t42 * t29 * t28 * t71 * t32 +
             t72 * t7 * t71 * t13 * t32 - t72 * t28 * t71 * t17 -
             t72 * t31 * t55 * t27 - t72 * t31 * t2 * t37 * t23 * t41;
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  }
  throw std::runtime_error("unexpected derivative selector: " + deriv);
}

// The VineCopula h-function kernels differentiate h(u|v) = dC/dv,
// conditioning on their second argument v, while hfunc1(u1, u2) = dC/du1
// conditions on the first; hence their u := our u2 and their v := our u1
// below, and their "_v" derivative corresponds to our "u1".
inline Eigen::VectorXd
JoeBicop::hfunc1_deriv_raw(const Eigen::MatrixXd& u,
                           const Eigen::MatrixXd& parameters,
                           const std::string& deriv)
{
  if (deriv == "par1") {
    // ported from VineCopula hfuncderiv.c diffhfunc (family 6 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      double theta = par(0);
      double t1 = 1.0 - u2;
      double t2 = std::pow(t1, 1.0 * theta);
      double t3 = 1.0 - u1;
      double t4 = std::pow(t3, 1.0 * theta);
      double t5 = t2 * t4;
      double t6 = t2 + t4 - t5;
      double t8 = 1 / theta - 1.0;
      double t9 = std::pow(t6, 1.0 * t8);
      double t10 = theta * theta;
      double t12 = std::log(t6);
      double t14 = std::log(t1);
      double t15 = t2 * t14;
      double t16 = std::log(t3);
      double t27 = std::pow(t3, 1.0 * theta - 1.0);
      double t7 = 1.0 - t2;
      double t11 = t9 * t27;
      return t9 *
               (-1.0 / t10 * t12 +
                t8 * (t15 + t4 * t16 - t15 * t4 - t5 * t16) / t6) *
               t27 * t7 +
             t11 * t16 * t7 - t11 * t15;
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  } else if (deriv == "u1") {
    // ported from VineCopula hfuncderiv.c diffhfunc_v (family 6 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      double theta = par(0);
      double t2 = std::pow(1.0 - u2, 1.0 * theta);
      double t3 = 1.0 - u1;
      double t4 = std::pow(t3, 1.0 * theta);
      double t5 = t2 * t4;
      double t6 = t2 + t4 - t5;
      double t8 = 1 / theta - 1.0;
      double t9 = std::pow(t6, 1.0 * t8);
      double t12 = 1 / t3;
      double t19 = theta - 1.0;
      double t20 = std::pow(t3, 1.0 * t19);
      double t22 = 1.0 - t2;
      return t9 * t8 * (-t4 * theta * t12 + t5 * theta * t12) / t6 * t20 * t22 -
             t9 * t20 * t19 * t12 * t22;
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  }
  throw std::runtime_error("unexpected derivative selector: " + deriv);
}

inline Eigen::VectorXd
JoeBicop::hfunc1_deriv2_raw(const Eigen::MatrixXd& u,
                            const Eigen::MatrixXd& parameters,
                            const std::string& deriv)
{
  if (deriv == "par1par1") {
    // ported from VineCopula hfuncderiv2.c diff2hfunc (family 6 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      double theta = par(0);
      double t1 = 1.0 - u2;
      double t2 = std::pow(t1, 1.0 * theta);
      double t3 = 1.0 - u1;
      double t4 = std::pow(t3, 1.0 * theta);
      double t5 = t2 * t4;
      double t6 = t2 + t4 - t5;
      double t8 = 1 / theta - 1.0;
      double t9 = std::pow(t6, 1.0 * t8);
      double t10 = theta * theta;
      double t11 = 1 / t10;
      double t12 = std::log(t6);
      double t14 = std::log(t1);
      double t15 = t2 * t14;
      double t16 = std::log(t3);
      double t18 = t4 * t16;
      double t20 = t15 + t18 - t15 * t4 - t5 * t16;
      double t21 = 1 / t6;
      double t23 = -t11 * t12 + t8 * t20 * t21;
      double t25 = t23 * t23;
      double t28 = std::pow(t3, 1.0 * theta - 1.0);
      double t29 = 1.0 - t2;
      double t30 = t28 * t29;
      double t39 = t14 * t14;
      double t40 = t2 * t39;
      double t42 = t16 * t16;
      double t50 = t20 * t20;
      double t56 = t6 * t6;
      double t13 = t9 * t23;
      double t7 = t9 * t28;
      return t9 * t25 * t30 +
             t9 *
               (2.0 / t10 / theta * t12 - 2.0 * t11 * t20 * t21 +
                t8 * (t40 + t4 * t42 - t40 * t4 - 2.0 * t15 * t18 - t5 * t42) *
                  t21 -
                t8 * t50 / t56) *
               t30 +
             2.0 * t13 * t28 * t16 * t29 - 2.0 * t13 * t28 * t2 * t14 +
             t7 * t42 * t29 - 2.0 * t7 * t16 * t2 * t14 - t7 * t40;
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  } else if (deriv == "par1u1") {
    // ported from VineCopula hfuncderiv2.c diff2hfunc_par_v (family 6
    // branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      double theta = par(0);
      double t1 = 1.0 - u2;
      double t2 = std::pow(t1, 1.0 * theta);
      double t3 = 1.0 - u1;
      double t4 = std::pow(t3, 1.0 * theta);
      double t5 = t2 * t4;
      double t6 = t2 + t4 - t5;
      double t8 = 1 / theta - 1.0;
      double t9 = std::pow(t6, 1.0 * t8);
      double t11 = t4 * theta;
      double t12 = 1 / t3;
      double t13 = t11 * t12;
      double t14 = theta * t12;
      double t16 = -t13 + t5 * t14;
      double t17 = t9 * t8 * t16;
      double t18 = 1 / t6;
      double t19 = theta * theta;
      double t20 = 1 / t19;
      double t21 = std::log(t6);
      double t23 = std::log(t1);
      double t24 = t2 * t23;
      double t25 = std::log(t3);
      double t30 = t8 * (t24 + t4 * t25 - t24 * t4 - t5 * t25);
      double t32 = -t20 * t21 + t30 * t18;
      double t34 = theta - 1.0;
      double t35 = std::pow(t3, 1.0 * t34);
      double t36 = 1.0 - t2;
      double t37 = t35 * t36;
      double t42 = t12 * t25;
      double t57 = t6 * t6;
      double t64 = t18 * t35;
      double t67 = t9 * t35;
      double t69 = t67 * t34;
      return t17 * t18 * t32 * t37 +
             t9 *
               (-t20 * t16 * t18 +
                t8 *
                  (-t11 * t42 - t4 * t12 + t24 * t13 + t5 * t14 * t25 +
                   t5 * t12) *
                  t18 -
                t30 / t57 * t16) *
               t37 -
             t9 * t32 * t35 * t34 * t12 * t36 + t17 * t64 * t25 * t36 -
             t69 * t42 * t36 - t67 * t12 * t36 - t17 * t64 * t24 +
             t69 * t12 * t2 * t23;
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  } else if (deriv == "u1u1") {
    // ported from VineCopula hfuncderiv2.c diff2hfunc_v (family 6 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      double theta = par(0);
      double t2 = std::pow(1.0 - u2, 1.0 * theta);
      double t3 = 1.0 - u1;
      double t4 = std::pow(t3, 1.0 * theta);
      double t5 = t2 * t4;
      double t6 = t2 + t4 - t5;
      double t8 = 1 / theta - 1.0;
      double t9 = std::pow(t6, 1.0 * t8);
      double t10 = t8 * t8;
      double t12 = t4 * theta;
      double t13 = 1 / t3;
      double t14 = -t12 * t13 + t5 * theta * t13;
      double t18 = t14 * t14;
      double t20 = t6 * t6;
      double t22 = theta - 1.0;
      double t23 = std::pow(t3, 1.0 * t22);
      double t24 = 1.0 - t2;
      double t26 = 1 / t20 * t23 * t24;
      double t27 = t9 * t8;
      double t29 = theta * theta;
      double t31 = t3 * t3;
      double t32 = 1 / t31;
      double t41 = 1 / t6;
      double t51 = t9 * t23;
      double t55 = t22 * t22;
      return t9 * t10 * t18 * t26 +
             t27 *
               (t4 * t29 * t32 - t12 * t32 - t5 * t29 * t32 +
                t5 * theta * t32) *
               t41 * t23 * t24 -
             t27 * t18 * t26 - 2.0 * t27 * t14 * t41 * t23 * t22 * t13 * t24 +
             t51 * t55 * t32 * t24 - t51 * t22 * t32 * t24;
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  }
  throw std::runtime_error("unexpected derivative selector: " + deriv);
}

inline Eigen::VectorXd
JoeBicop::logpdf_deriv_raw(const Eigen::MatrixXd& u,
                           const Eigen::MatrixXd& parameters,
                           const std::string& deriv)
{
  // fused single-pass form kept for performance: the base
  // logpdf_deriv*_raw would compose this as (pdf_deriv)/pdf, recomputing the
  // shared pdf temporaries 2-3x; this leaf is on the scores/Hessian hot paths.
  if (deriv != "par1") {
    return AbstractBicop::logpdf_deriv_raw(u, parameters, deriv);
  }
  // ported from VineCopula logderiv.c difflPDF (family 6 branch)
  auto f = [](const double& u1,
              const double& u2,
              const Eigen::Ref<const Eigen::VectorXd>& par) {
    double theta = par(0);
    double t1 = 1.0 - u1;
    double t2 = std::pow(t1, 1.0 * theta);
    double t3 = 1.0 - u2;
    double t4 = std::pow(t3, 1.0 * theta);
    double t5 = t2 * t4;
    double t6 = t2 + t4 - t5;
    double t8 = 1 / theta - 2.0;
    double t9 = std::pow(t6, 1.0 * t8);
    double t10 = theta * theta;
    double t12 = std::log(t6);
    double t14 = std::log(t1);
    double t15 = t2 * t14;
    double t16 = std::log(t3);
    double t17 = t4 * t16;
    double t18 = t15 * t4;
    double t19 = t5 * t16;
    double t26 = theta - 1.0;
    double t27 = std::pow(t1, 1.0 * t26);
    double t28 = std::pow(t3, 1.0 * t26);
    double t30 = theta - 1.0 + t2 + t4 - t5;
    double t33 = t9 * t27;
    return (t9 * (-1 / t10 * t12 + t8 * (t15 + t17 - t18 - t19) / t6) * t27 *
              t28 * t30 +
            t33 * t14 * t28 * t30 + t33 * t28 * t16 * t30 +
            t33 * t28 * (1.0 + t15 + t17 - t18 - t19)) /
           t9 / t27 / t28 / t30;
  };
  return tools_eigen::binaryExpr_or_nan(u, parameters, f);
}

inline Eigen::VectorXd
JoeBicop::logpdf_deriv2_raw(const Eigen::MatrixXd& u,
                            const Eigen::MatrixXd& parameters,
                            const std::string& deriv)
{
  // fused single-pass form kept for performance: the base
  // logpdf_deriv*_raw would compose this as (pdf_deriv)/pdf, recomputing the
  // shared pdf temporaries 2-3x; this leaf is on the scores/Hessian hot paths.
  if (deriv != "par1par1") {
    return AbstractBicop::logpdf_deriv2_raw(u, parameters, deriv);
  }
  // ported from VineCopula logderiv.c diff2lPDF (family 6 branch)
  auto f = [](const double& u1,
              const double& u2,
              const Eigen::Ref<const Eigen::VectorXd>& par) {
    double theta = par(0);
    double t1 = 1.0 - u1;
    double t2 = std::pow(t1, 1.0 * theta);
    double t3 = 1.0 - u2;
    double t4 = std::pow(t3, 1.0 * theta);
    double t5 = t2 * t4;
    double t6 = t2 + t4 - t5;
    double t8 = 1 / theta - 2.0;
    double t9 = std::pow(t6, 1.0 * t8);
    double t10 = theta * theta;
    double t11 = 1 / t10;
    double t12 = std::log(t6);
    double t14 = std::log(t1);
    double t15 = t2 * t14;
    double t16 = std::log(t3);
    double t17 = t4 * t16;
    double t18 = t15 * t4;
    double t19 = t5 * t16;
    double t20 = t15 + t17 - t18 - t19;
    double t22 = 1 / t6;
    double t24 = -t11 * t12 + t8 * t20 * t22;
    double t25 = t24 * t24;
    double t27 = theta - 1.0;
    double t28 = std::pow(t1, 1.0 * t27);
    double t29 = std::pow(t3, 1.0 * t27);
    double t30 = t28 * t29;
    double t31 = theta - 1.0 + t2 + t4 - t5;
    double t32 = t30 * t31;
    double t41 = t14 * t14;
    double t42 = t2 * t41;
    double t43 = t16 * t16;
    double t49 = t42 + t4 * t43 - t42 * t4 - 2.0 * t15 * t17 - t5 * t43;
    double t52 = t20 * t20;
    double t54 = t6 * t6;
    double t60 = t9 * t24;
    double t61 = t60 * t28;
    double t62 = t14 * t29;
    double t63 = t62 * t31;
    double t66 = t29 * t16;
    double t67 = t66 * t31;
    double t70 = 1.0 + t15 + t17 - t18 - t19;
    double t74 = t9 * t28;
    double t92 = t9 * t25 * t32 +
                 t9 *
                   (2.0 / t10 / theta * t12 - 2.0 * t11 * t20 * t22 +
                    t8 * t49 * t22 - t8 * t52 / t54) *
                   t32 +
                 2.0 * t61 * t63 + 2.0 * t61 * t67 + 2.0 * t60 * t30 * t70 +
                 t74 * t41 * t29 * t31 + 2.0 * t74 * t14 * t67 +
                 2.0 * t74 * t62 * t70 + t74 * t29 * t43 * t31 +
                 2.0 * t74 * t66 * t70 + t74 * t29 * t49;
    double t93 = 1 / t9;
    double t95 = 1 / t28;
    double t96 = 1 / t29;
    double t98 = 1 / t31;
    double t108 =
      (t60 * t32 + t74 * t63 + t74 * t67 + t74 * t29 * t70) * t93 * t95;
    double t109 = t96 * t98;
    double t116 = t31 * t31;
    return t92 * t93 * t95 * t96 * t98 - t108 * t109 * t24 - t108 * t109 * t14 -
           t108 * t109 * t16 - t108 * t96 / t116 * t70;
  };
  return tools_eigen::binaryExpr_or_nan(u, parameters, f);
}

// link between Kendall's tau and the par_bicop parameter
inline Eigen::MatrixXd
JoeBicop::tau_to_parameters(const double& tau)
{
  Eigen::VectorXd tau0 = Eigen::VectorXd::Constant(1, std::fabs(tau));
  auto f = [&](const Eigen::VectorXd& v) {
    return Eigen::VectorXd::Constant(1, std::fabs(parameters_to_tau(v)));
  };
  return tools_eigen::invert_f(tau0,
                               f,
                               parameters_lower_bounds_(0) + 1e-6,
                               parameters_upper_bounds_(0) - 1e-6);
}

inline double
JoeBicop::parameters_to_tau(const Eigen::MatrixXd& parameters)
{
  double par = parameters(0);
  // The closed form 1 + 2 [psi(2) - psi(2/par + 1)] / (2 - par) has a removable
  // singularity at par = 2, where the bracket and the denominator vanish
  // together: the quotient evaluates to NaN there and loses most of its
  // significant digits nearby. The limit is 1 - psi'(2) = 2 - pi^2 / 6.
  constexpr double eps = 1e-5;
  if (std::fabs(par - 2.0) < eps) {
    // tau(2 + d) = 1 - psi'(2) + d [psi'(2) / 2 + psi''(2) / 4] + O(d^2)
    constexpr double tri2 = 0.6449340668482264;  // psi'(2) = pi^2 / 6 - 1
    constexpr double tet2 = -0.4041138063191886; // psi''(2)
    const double d = par - 2.0;
    return 1.0 - tri2 + d * (tri2 / 2.0 + tet2 / 4.0);
  }
  double tau = 2 / par + 1;
  tau = boost::math::digamma(2.0) - boost::math::digamma(tau);
  return 1 + 2 * tau / (2 - par);
}

inline Eigen::MatrixXd
JoeBicop::parameters_to_taildep(const Eigen::MatrixXd& par)
{
  Eigen::MatrixXd taildep = Eigen::MatrixXd::Zero(2, 2);
  taildep(1, 1) = 2 - std::pow(2.0, 1.0 / par(0)); // upper tail dependence
  return taildep;
}

inline Eigen::VectorXd
JoeBicop::get_start_parameters(const double tau)
{
  Eigen::VectorXd par = tau_to_parameters(tau);
  par = par.cwiseMax(parameters_lower_bounds_);
  par = par.cwiseMin(parameters_upper_bounds_);
  return par;
}
}

// This is copy&paste from the VineCopula package
inline double
qcondjoe(const double& q, const double& u, const double& de)
{
  double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t13, t15, t16, t19, t23,
    t28, t31;
  double c21, pdf;
  int iter;
  double diff, v, de1, dtem, de1inv, tem;

  t1 = 1.0 - u;
  t2 = std::pow(t1, 1.0 * (de));
  t7 = 1. / (de);
  t10 = t2 * (de);
  t11 = 1. / t1;
  t19 = (de) * (de);
  de1 = de - 1; // may need better modification for large delta
  dtem = -de1 / (1. + de1);
  de1inv = -1. / de1;

  // v = 0.5 * (q+u); // starting guess

  // Use a better starting point based on reflected B4 copula
  // A good starting point is crucial when delta is large because
  //    C_{2|1} will be steep
  // C_{R,2|1}(v|u)=1-C_{2|1}(1-v|1-u),
  // C_{R,2|1}^{-1}(q|u)=1-C_{2|1}^{-1}(1-q|1-u)
  tem = std::pow(1. - q, dtem) - 1.;
  tem = tem * std::pow(1. - u, -de1) + 1.;
  v = std::pow(tem, de1inv);
  v = 1. - v;
  diff = 1;
  iter = 0;
  while (fabs(diff) > 1.e-6 && iter < 20) {
    t3 = 1. - v;
    t4 = std::pow(t3, de);
    t5 = t2 * t4;
    t6 = t2 + t4 - t5;
    t8 = std::pow(t6, t7);
    t9 = t7 * t8;
    t13 = t11 * t4;
    t15 = -t10 * t11 + t10 * t13;
    t16 = 1. / t6;
    t23 = 1. / t3;
    t28 = t6 * t6;
    t31 = (-t4 * (de)*t23 + t5 * (de)*t23) / t28 * t15;
    c21 = -t9 * t15 * t16;
    pdf = -t8 / t19 * t31 + t8 * (de)*t2 * t13 * t23 * t16 + t9 * t31;
    iter++;
    if ((std::isnan)(pdf) || (std::isnan)(c21)) {
      diff /= -2.;
    } // added for de>=30
    else
      diff = (c21 - q) / pdf;
    v -= diff;
    int iter2 = 0;
    while ((v <= 0 || v >= 1 || fabs(diff) > 0.25) & (iter2 < 20)) {
      ++iter2;
      diff /= 2.;
      v += diff;
    }
  }

  // make sure that boundaries are respected
  if (v <= 0) {
    v = 1e-10;
  } else if (v >= 1) {
    v = 1 - 1e-10;
  }

  return v;
}
