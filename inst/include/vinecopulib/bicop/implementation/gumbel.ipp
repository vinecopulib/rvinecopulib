// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <cmath>
#include <vinecopulib/misc/tools_eigen.hpp>

namespace vinecopulib {
inline GumbelBicop::GumbelBicop()
{
  family_ = BicopFamily::gumbel;
  parameters_ = Eigen::VectorXd(1);
  parameters_lower_bounds_ = Eigen::VectorXd(1);
  parameters_upper_bounds_ = Eigen::VectorXd(1);
  parameters_ << 1;
  parameters_lower_bounds_ << 1;
  parameters_upper_bounds_ << 50;
}

inline double
GumbelBicop::generator(const double& u,
                       const Eigen::Ref<const Eigen::VectorXd>& parameters)
{
  return std::pow(std::log(1 / u), parameters(0));
}

inline double
GumbelBicop::generator_inv(const double& u,
                           const Eigen::Ref<const Eigen::VectorXd>& parameters)
{
  return std::exp(-std::pow(u, 1 / parameters(0)));
}

inline double
GumbelBicop::generator_derivative(
  const double& u,
  const Eigen::Ref<const Eigen::VectorXd>& parameters)
{
  double theta = parameters(0);
  return std::pow(std::log(1 / u), theta - 1) * (-theta / u);
}

inline Eigen::VectorXd
GumbelBicop::pdf_raw(const Eigen::MatrixXd& u,
                     const Eigen::MatrixXd& parameters)
{
  auto f = [](const double& u1,
              const double& u2,
              const Eigen::Ref<const Eigen::VectorXd>& par) {
    double theta = par(0);
    double thetha1 = 1.0 / theta;
    double t1 = std::pow(-std::log(u1), theta) + std::pow(-std::log(u2), theta);
    double temp = -std::pow(t1, thetha1) + (2 * thetha1 - 2.0) * std::log(t1) +
                  (theta - 1.0) * std::log(std::log(u1) * std::log(u2)) -
                  std::log(u1 * u2) +
                  std::log1p((theta - 1.0) * std::pow(t1, -thetha1));
    return std::exp(temp);
  };
  return tools_eigen::binaryExpr_or_nan(u, parameters, f);
}

inline Eigen::VectorXd
GumbelBicop::hinv1_raw(const Eigen::MatrixXd& u,
                       const Eigen::MatrixXd& parameters)
{
  auto qcondgum_func =
    [](const double& u1,
       const double& u2,
       const Eigen::Ref<const Eigen::VectorXd>& par) -> double {
    return qcondgum(u2, u1, par(0));
  };
  return tools_eigen::binaryExpr_or_nan(u, parameters, qcondgum_func);
}

inline Eigen::VectorXd
GumbelBicop::pdf_deriv_raw(const Eigen::MatrixXd& u,
                           const Eigen::MatrixXd& parameters,
                           const std::string& deriv)
{
  if (tools_deriv::is_u2_only(deriv)) {
    // exchangeability: c(u1, u2) = c(u2, u1)
    return pdf_deriv_raw(
      tools_eigen::swap_cols(u), parameters, tools_deriv::swap_args(deriv));
  }
  if (deriv == "par1") {
    // ported from VineCopula deriv.c diffPDF (family 4 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      double theta = par(0);
      double t3 = std::log(u1);
      double t4 = std::pow(-t3, 1.0 * theta);
      double t5 = std::log(u2);
      double t6 = std::pow(-t5, 1.0 * theta);
      double t7 = t4 + t6;
      double t8 = 1 / theta;
      double t9 = std::pow(t7, 1.0 * t8);
      double t10 = theta * theta;
      double t12 = std::log(t7);
      double t13 = 1 / t10 * t12;
      double t14 = std::log(-t3);
      double t16 = std::log(-t5);
      double t18 = t4 * t14 + t6 * t16;
      double t20 = 1 / t7;
      double t22 = -t13 + t8 * t18 * t20;
      double t24 = std::exp(-t9);
      double t26 = t24 / u1;
      double t28 = 1 / u2;
      double t29 = -1.0 + t8;
      double t30 = std::pow(t7, 2.0 * t29);
      double t32 = t3 * t5;
      double t33 = theta - 1.0;
      double t34 = std::pow(t32, 1.0 * t33);
      double t35 = std::pow(t7, -1.0 * t8);
      double t36 = t33 * t35;
      double t17 = 1.0 + t36;
      double t15 = t34 * t17;
      double t11 = t26 * t28;
      double t2 = t30 * t34;
      double t1 = std::log(t32);
      return -t9 * t22 * t26 * t28 * t30 * t15 +
             t11 * t30 * (-2.0 * t13 + 2.0 * t29 * t18 * t20) * t15 +
             t11 * t2 * t1 * t17 + t11 * t2 * (t35 - t36 * t22);
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  } else if (deriv == "u1") {
    // ported from VineCopula deriv.c diffPDF_u (family 4 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      double theta = par(0);
      double t3 = std::log(u1);
      double t4 = std::pow(-t3, 1.0 * theta);
      double t5 = std::log(u2);
      double t6 = std::pow(-t5, 1.0 * theta);
      double t7 = t4 + t6;
      double t8 = 1 / theta;
      double t9 = std::pow(t7, 1.0 * t8);
      double t11 = u1 * u1;
      double t12 = 1 / t11;
      double t13 = 1 / t3;
      double t15 = 1 / t7;
      double t18 = std::exp(-t9);
      double t19 = 1 / u2;
      double t21 = -1.0 + t8;
      double t22 = std::pow(t7, 2.0 * t21);
      double t24 = theta - 1.0;
      double t25 = std::pow(t3 * t5, 1.0 * t24);
      double t27 = std::pow(t7, -1.0 * t8);
      double t28 = t24 * t27;
      double t29 = 1.0 + t28;
      double t30 = t22 * t25 * t29;
      double t33 = t18 * t12;
      double t36 = t19 * t22;
      return -t9 * t4 * t12 * t13 * t15 * t18 * t19 * t30 - t33 * t19 * t30 +
             2.0 * t33 * t36 * t21 * t4 * theta * t13 * t15 * t25 * t29 +
             t33 * t36 * t25 * t24 * t13 * t29 -
             t33 * t36 * t25 * t28 * t4 * t13 * t15;
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  }
  throw std::runtime_error("unexpected derivative selector: " + deriv);
}

inline Eigen::VectorXd
GumbelBicop::pdf_deriv2_raw(const Eigen::MatrixXd& u,
                            const Eigen::MatrixXd& parameters,
                            const std::string& deriv)
{
  if (tools_deriv::is_u2_only(deriv)) {
    // exchangeability: c(u1, u2) = c(u2, u1)
    return pdf_deriv2_raw(
      tools_eigen::swap_cols(u), parameters, tools_deriv::swap_args(deriv));
  }
  if (deriv == "par1par1") {
    // ported from VineCopula deriv2.c diff2PDF (family 4 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      double theta = par(0);
      double t3 = std::log(u1);
      double t4 = std::pow(-t3, 1.0 * theta);
      double t5 = std::log(u2);
      double t6 = std::pow(-t5, 1.0 * theta);
      double t7 = t4 + t6;
      double t8 = 1 / theta;
      double t9 = std::pow(t7, 1.0 * t8);
      double t10 = theta * theta;
      double t11 = 1 / t10;
      double t12 = std::log(t7);
      double t13 = t11 * t12;
      double t14 = std::log(-t3);
      double t16 = std::log(-t5);
      double t18 = t4 * t14 + t6 * t16;
      double t20 = 1 / t7;
      double t22 = -t13 + t8 * t18 * t20;
      double t23 = t22 * t22;
      double t25 = std::exp(-t9);
      double t27 = t25 / u1;
      double t29 = 1 / u2;
      double t30 = -1.0 + t8;
      double t31 = std::pow(t7, 2.0 * t30);
      double t32 = t29 * t31;
      double t33 = t3 * t5;
      double t34 = theta - 1.0;
      double t35 = std::pow(t33, 1.0 * t34);
      double t36 = std::pow(t7, -1.0 * t8);
      double t37 = t34 * t36;
      double t38 = 1.0 + t37;
      double t39 = t35 * t38;
      double t40 = t32 * t39;
      double t44 = 1 / t10 / theta * t12;
      double t47 = t11 * t18 * t20;
      double t49 = t14 * t14;
      double t51 = t16 * t16;
      double t53 = t4 * t49 + t6 * t51;
      double t56 = t18 * t18;
      double t58 = t7 * t7;
      double t59 = 1 / t58;
      double t61 = 2.0 * t44 - 2.0 * t47 + t8 * t53 * t20 - t8 * t56 * t59;
      double t65 = t9 * t9;
      double t70 = t9 * t22 * t27;
      double t74 = -2.0 * t13 + 2.0 * t30 * t18 * t20;
      double t75 = t74 * t35;
      double t80 = std::log(t33);
      double t87 = t36 - t37 * t22;
      double t88 = t35 * t87;
      double t17 = t27 * t29;
      double t15 = t74 * t74;
      double t2 = t31 * t35;
      double t1 = t80 * t80;
      return -t9 * t23 * t27 * t40 - t9 * t61 * t27 * t40 +
             t65 * t23 * t27 * t40 - 2.0 * t70 * t32 * t75 * t38 -
             2.0 * t70 * t32 * t35 * t80 * t38 - 2.0 * t70 * t32 * t88 +
             t17 * t31 * t15 * t39 +
             t17 * t31 *
               (4.0 * t44 - 4.0 * t47 + 2.0 * t30 * t53 * t20 -
                2.0 * t30 * t56 * t59) *
               t39 +
             2.0 * t27 * t32 * t75 * t80 * t38 + 2.0 * t17 * t31 * t74 * t88 +
             t17 * t2 * t1 * t38 + 2.0 * t17 * t2 * t80 * t87 +
             t17 * t2 * (-2.0 * t36 * t22 + t37 * t23 - t37 * t61);
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  } else if (deriv == "par1u1") {
    // ported from VineCopula deriv2.c diff2PDF_par_u (family 4 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      double theta = par(0);
      double t3 = std::log(u1);
      double t4 = std::pow(-t3, 1.0 * theta);
      double t5 = std::log(u2);
      double t6 = std::pow(-t5, 1.0 * theta);
      double t7 = t4 + t6;
      double t8 = 1 / theta;
      double t9 = std::pow(t7, 1.0 * t8);
      double t11 = u1 * u1;
      double t12 = 1 / t11;
      double t13 = 1 / t3;
      double t15 = 1 / t7;
      double t17 = t9 * t4 * t12 * t13 * t15;
      double t18 = theta * theta;
      double t20 = std::log(t7);
      double t21 = 1 / t18 * t20;
      double t22 = std::log(-t3);
      double t24 = std::log(-t5);
      double t26 = t4 * t22 + t6 * t24;
      double t29 = -t21 + t8 * t26 * t15;
      double t30 = std::exp(-t9);
      double t32 = 1 / u2;
      double t34 = -1.0 + t8;
      double t35 = std::pow(t7, 2.0 * t34);
      double t36 = t3 * t5;
      double t37 = theta - 1.0;
      double t38 = std::pow(t36, 1.0 * t37);
      double t39 = t35 * t38;
      double t40 = std::pow(t7, -1.0 * t8);
      double t41 = t37 * t40;
      double t42 = t41 + 1.0;
      double t43 = t39 * t42;
      double t47 = 1 / u1;
      double t48 = t47 * t13;
      double t49 = t48 * t15;
      double t50 = t8 * t4 * t49;
      double t51 = t4 * theta;
      double t55 = t4 * t47 * t13;
      double t56 = t51 * t48 * t22 + t55;
      double t59 = t7 * t7;
      double t60 = 1 / t59;
      double t63 = -t50 + t8 * t56 * t15 - t26 * t60 * t55;
      double t65 = t30 * t47;
      double t67 = t32 * t35;
      double t68 = t38 * t42;
      double t69 = t67 * t68;
      double t71 = t9 * t9;
      double t80 = t9 * t29;
      double t81 = t30 * t12;
      double t87 = t80 * t30 * t12 * t32 * t35;
      double t94 = t81 * t32;
      double t97 = t37 * t13 * t42;
      double t100 = t38 * t37;
      double t103 = t4 * t13 * t15;
      double t104 = t100 * t40 * t103;
      double t106 = t30 * t32;
      double t107 = t106 * t35;
      double t109 = 2.0 * t34 * t26;
      double t111 = -2.0 * t21 + t109 * t15;
      double t112 = t111 * t38;
      double t113 = t112 * t42;
      double t121 = 2.0 * t94 * t35 * t34 * t4;
      double t123 = theta * t13 * t15;
      double t126 = t65 * t32;
      double t137 = t81 * t67;
      double t140 =
        -t17 * t29 * t30 * t32 * t43 - t9 * t63 * t65 * t69 +
        t71 * t29 * t4 * t12 * t13 * t15 * t30 * t32 * t43 + t80 * t81 * t69 -
        2.0 * t87 * t34 * t4 * theta * t13 * t15 * t68 - t80 * t94 * t39 * t97 +
        t87 * t104 - t17 * t107 * t113 - t94 * t35 * t111 * t68 +
        t121 * t123 * t113 +
        t126 * t35 *
          (-2.0 * t50 + 2.0 * t34 * t56 * t15 - t109 * t60 * t51 * t48) * t68 +
        t137 * t112 * t97;
      double t144 = std::log(t36);
      double t146 = t38 * t144 * t42;
      double t10 = t40 - t41 * t29;
      double t2 = t39 * t10;
      double t1 =
        -t81 * t67 * t111 * t104 - t17 * t107 * t146 - t94 * t39 * t144 * t42 +
        t121 * t123 * t146 + t137 * t100 * t13 * t144 * t42 +
        t94 * t39 * t13 * t42 - t81 * t67 * t38 * t144 * t37 * t40 * t103 -
        t17 * t106 * t2 - t94 * t2 +
        2.0 * t81 * t67 * t34 * t51 * t13 * t15 * t38 * t10 +
        t137 * t100 * t13 * t10 +
        t126 * t39 * (-t40 * t4 * t49 + t41 * t4 * t48 * t15 * t29 - t41 * t63);
      return t140 + t1;
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  } else if (deriv == "u1u1") {
    // ported from VineCopula deriv2.c diff2PDF_u (family 4 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      double theta = par(0);
      double t3 = std::log(u1);
      double t4 = std::pow(-t3, 1.0 * theta);
      double t5 = std::log(u2);
      double t6 = std::pow(-t5, 1.0 * theta);
      double t7 = t4 + t6;
      double t8 = 1 / theta;
      double t9 = std::pow(t7, 1.0 * t8);
      double t10 = std::exp(-t9);
      double t11 = u1 * u1;
      double t13 = 1 / t11 / u1;
      double t14 = t10 * t13;
      double t15 = 1 / u2;
      double t16 = t14 * t15;
      double t17 = -1.0 + t8;
      double t18 = std::pow(t7, 2.0 * t17);
      double t20 = theta - 1.0;
      double t21 = std::pow(t3 * t5, 1.0 * t20);
      double t23 = std::pow(t7, -1.0 * t8);
      double t24 = t20 * t23;
      double t25 = t24 + 1.0;
      double t26 = t18 * t21 * t25;
      double t29 = t9 * t4;
      double t30 = 1 / t3;
      double t32 = 1 / t7;
      double t35 = t10 * t15;
      double t36 = t35 * t26;
      double t39 = t3 * t3;
      double t40 = 1 / t39;
      double t41 = t13 * t40;
      double t43 = t29 * t41 * t32;
      double t45 = t15 * t18;
      double t46 = t14 * t45;
      double t47 = t21 * t20;
      double t52 = t4 * t4;
      double t53 = t9 * t52;
      double t54 = t7 * t7;
      double t55 = 1 / t54;
      double t56 = t41 * t55;
      double t57 = t53 * t56;
      double t66 = t35 * t18;
      double t68 = t21 * t25 * theta;
      double t71 = t9 * t9;
      double t79 = 2.0 * t45 * t17;
      double t83 = t47 * t25;
      double t87 = t47 * t23;
      double t91 = t14 * t79;
      double t92 = t4 * theta;
      double t95 = t32 * t21 * t25;
      double t100 = t14 * t45 * t21;
      double t106 =
        2.0 * t16 * t26 + 3.0 * t29 * t13 * t30 * t32 * t36 + t43 * t36 -
        3.0 * t46 * t47 * t30 * t25 - t57 * t36 -
        t29 * theta * t13 * t40 * t32 * t10 * t15 * t26 + t57 * t66 * t68 +
        t71 * t52 * t56 * t36 - 2.0 * t53 * t13 * t40 * t55 * t10 * t79 * t68 -
        2.0 * t43 * t66 * t83 + 2.0 * t57 * t66 * t87 -
        3.0 * t91 * t92 * t30 * t95 + 3.0 * t100 * t24 * t4 * t30 * t32;
      double t113 = theta * theta;
      double t118 = t52 * t113 * t40 * t55 * t21 * t25;
      double t125 = 2.0 * t18 * t17;
      double t128 = theta * t40;
      double t129 = t128 * t32;
      double t19 = t128 * t55;
      double t12 = t40 * t25;
      double t2 = t20 * t20;
      double t1 = -t91 * t92 * t40 * t95 + 4.0 * t14 * t45 * t17 * t17 * t118 +
                  t91 * t4 * t113 * t40 * t95 - t91 * t118 +
                  2.0 * t16 * t125 * t4 * t129 * t83 -
                  2.0 * t16 * t125 * t52 * t19 * t87 - t46 * t47 * t12 +
                  t46 * t21 * t2 * t12 -
                  2.0 * t100 * t2 * t40 * t23 * t4 * t32 +
                  t100 * t24 * t4 * t40 * t32 + t100 * t24 * t52 * t40 * t55 -
                  t100 * t24 * t4 * t129 + t100 * t24 * t52 * t19;
      return t106 + t1;
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  } else if (deriv == "u1u2") {
    // ported from VineCopula deriv2.c diff2PDF_u_v (family 4 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      double theta = par(0);
      double t3 = std::log(u1);
      double t4 = std::pow(-t3, 1.0 * theta);
      double t5 = std::log(u2);
      double t6 = std::pow(-t5, 1.0 * theta);
      double t7 = t4 + t6;
      double t8 = 1 / theta;
      double t9 = std::pow(t7, 1.0 * t8);
      double t10 = std::exp(-t9);
      double t11 = u1 * u1;
      double t12 = 1 / t11;
      double t13 = t10 * t12;
      double t14 = u2 * u2;
      double t15 = 1 / t14;
      double t16 = t13 * t15;
      double t17 = -1.0 + t8;
      double t18 = std::pow(t7, 2.0 * t17);
      double t20 = theta - 1.0;
      double t21 = std::pow(t3 * t5, 1.0 * t20);
      double t22 = t18 * t21;
      double t23 = std::pow(t7, -1.0 * t8);
      double t24 = t20 * t23;
      double t25 = 1.0 + t24;
      double t26 = t22 * t25;
      double t28 = t9 * t4;
      double t29 = 1 / t3;
      double t30 = t12 * t29;
      double t31 = 1 / t7;
      double t34 = t10 * t15;
      double t37 = t9 * t6;
      double t38 = t37 * t15;
      double t39 = 1 / t5;
      double t40 = t7 * t7;
      double t41 = 1 / t40;
      double t48 = t28 * t12;
      double t49 = t29 * t41;
      double t51 = t48 * t49 * t10;
      double t52 = t15 * t18;
      double t53 = t52 * t21;
      double t55 = theta * t39;
      double t59 = t9 * t9;
      double t64 = t15 * t39;
      double t70 = 2.0 * t18 * t17;
      double t71 = t70 * t6;
      double t72 = t21 * t25;
      double t84 = t6 * t39;
      double t85 = t24 * t84;
      double t94 = 2.0 * t13 * t52 * t17;
      double t98 = t31 * t21 * t25;
      double t101 = t13 * t52;
      double t102 = t21 * t20;
      double t104 = t102 * t39 * t25;
      double t106 = t13 * t53;
      double t107 = t84 * t31;
      double t110 = t4 * theta;
      double t114 =
        t16 * t26 + t28 * t30 * t31 * t34 * t26 -
        t38 * t39 * t41 * t4 * t30 * t10 * t26 + t51 * t53 * t25 * t6 * t55 +
        t59 * t4 * t12 * t49 * t6 * t64 * t10 * t26 -
        2.0 * t48 * t49 * t34 * t71 * t55 * t72 -
        t48 * t29 * t31 * t10 * t53 * t20 * t39 * t25 + 2.0 * t51 * t53 * t85 +
        t37 * t64 * t31 * t13 * t26 - t94 * t6 * theta * t39 * t98 -
        t101 * t104 + t106 * t24 * t107 - t94 * t110 * t29 * t98;
      double t118 = theta * theta;
      double t128 = t4 * t29;
      double t129 = t16 * t70 * t4;
      double t32 = theta * t29 * t31 * t104;
      double t27 = t20 * t20;
      double t19 = t128 * t31;
      double t2 = t16 * t22 * t20;
      double t1 =
        4.0 * t16 * t18 * t17 * t17 * t6 * t118 * t39 * t41 * t128 * t72 -
        t129 * t118 * t29 * t41 * t72 * t84 + t129 * t32 -
        2.0 * t16 * t70 * t110 * t49 * t21 * t85 - t101 * t102 * t29 * t25 -
        t38 * t39 * t31 * t10 * t12 * t18 * t21 * t20 * t29 * t25 +
        t16 * t71 * t32 + t101 * t21 * t27 * t39 * t29 * t25 -
        t106 * t27 * t29 * t23 * t107 + t106 * t24 * t19 -
        t106 * t27 * t39 * t23 * t19 + t2 * t23 * t6 * t39 * t41 * t4 * t29 +
        t2 * t23 * t4 * t29 * t41 * t6 * t55;
      return t114 + t1;
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  }
  throw std::runtime_error("unexpected derivative selector: " + deriv);
}

inline Eigen::VectorXd
GumbelBicop::hfunc1_deriv_raw(const Eigen::MatrixXd& u,
                              const Eigen::MatrixXd& parameters,
                              const std::string& deriv)
{
  // VineCopula's h-function kernels differentiate h(u|v) = dC(u, v)/dv;
  // our hfunc1(u1, u2) = dC(u1, u2)/du1, so their u is our u2 and their v
  // is our u1 in the lambdas below.
  if (deriv == "par1") {
    // ported from VineCopula hfuncderiv.c diffhfunc (family 4 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      double theta = par(0);
      double t1 = std::log(u1);
      double t2 = std::pow(-t1, 1.0 * theta);
      double t3 = std::log(u2);
      double t4 = std::pow(-t3, 1.0 * theta);
      double t5 = t2 + t4;
      double t6 = 1 / theta;
      double t7 = std::pow(t5, 1.0 * t6);
      double t8 = theta * theta;
      double t9 = std::log(t5);
      double t10 = 1 / t8 * t9;
      double t11 = std::log(-t1);
      double t14 = std::log(-t3);
      double t16 = t2 * t11 + t4 * t14;
      double t18 = 1 / t5;
      double t22 = std::exp(-t7);
      double t24 = t6 - 1.0;
      double t25 = std::pow(t5, 1.0 * t24);
      double t27 = 1 / u1;
      double t28 = 1 / t1;
      double t32 = t22 * t25;
      return t7 * (-t10 + t6 * t16 * t18) * t22 * t25 * t2 * t27 * t28 -
             t32 * (-t10 + t24 * t16 * t18) * t2 * t27 * t28 -
             t32 * t2 * t11 * t27 * t28;
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  } else if (deriv == "u1") {
    // ported from VineCopula hfuncderiv.c diffhfunc_v (family 4 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      double theta = par(0);
      double t3 = std::log(u2);
      double t4 = std::pow(-t3, 1.0 * theta);
      double t5 = std::log(u1);
      double t6 = std::pow(-t5, 1.0 * theta);
      double t7 = t4 + t6;
      double t8 = 1 / theta;
      double t9 = std::pow(t7, 1.0 * t8);
      double t10 = t6 * t6;
      double t12 = u1 * u1;
      double t13 = 1 / t12;
      double t15 = t5 * t5;
      double t16 = 1 / t15;
      double t18 = t16 / t7;
      double t19 = std::exp(-t9);
      double t20 = t8 - 1.0;
      double t21 = std::pow(t7, 1.0 * t20);
      double t22 = t19 * t21;
      double t27 = theta * t13;
      double t33 = t6 * t13;
      return t9 * t10 * t13 * t18 * t22 - t22 * t20 * t10 * t27 * t18 -
             t22 * t6 * t27 * t16 + t22 * t33 / t5 + t22 * t33 * t16;
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  }
  throw std::runtime_error("unexpected derivative selector: " + deriv);
}

inline Eigen::VectorXd
GumbelBicop::hfunc1_deriv2_raw(const Eigen::MatrixXd& u,
                               const Eigen::MatrixXd& parameters,
                               const std::string& deriv)
{
  // VineCopula's h-function kernels differentiate h(u|v) = dC(u, v)/dv;
  // our hfunc1(u1, u2) = dC(u1, u2)/du1, so their u is our u2 and their v
  // is our u1 in the lambdas below.
  if (deriv == "par1par1") {
    // ported from VineCopula hfuncderiv2.c diff2hfunc (family 4 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      double theta = par(0);
      double t1 = std::log(u1);
      double t2 = std::pow(-t1, 1.0 * theta);
      double t3 = std::log(u2);
      double t4 = std::pow(-t3, 1.0 * theta);
      double t5 = t2 + t4;
      double t6 = 1 / theta;
      double t7 = std::pow(t5, 1.0 * t6);
      double t8 = theta * theta;
      double t9 = 1 / t8;
      double t10 = std::log(t5);
      double t11 = t9 * t10;
      double t12 = std::log(-t1);
      double t13 = t2 * t12;
      double t14 = std::log(-t3);
      double t16 = t13 + t4 * t14;
      double t18 = 1 / t5;
      double t20 = -t11 + t6 * t16 * t18;
      double t21 = t20 * t20;
      double t23 = std::exp(-t7);
      double t25 = t6 - 1.0;
      double t26 = std::pow(t5, 1.0 * t25);
      double t28 = 1 / u1;
      double t29 = 1 / t1;
      double t30 = t28 * t29;
      double t31 = t26 * t2 * t30;
      double t36 = 2.0 / t8 / theta * t10;
      double t39 = 2.0 * t9 * t16 * t18;
      double t40 = t12 * t12;
      double t42 = t14 * t14;
      double t44 = t2 * t40 + t4 * t42;
      double t47 = t16 * t16;
      double t49 = t5 * t5;
      double t50 = 1 / t49;
      double t56 = t7 * t7;
      double t61 = t23 * t26;
      double t62 = t7 * t20 * t61;
      double t65 = -t11 + t25 * t16 * t18;
      double t70 = t13 * t30;
      double t73 = t65 * t65;
      double t15 = t2 * t28 * t29;
      return t7 * t21 * t23 * t31 +
             t7 * (t36 - t39 + t6 * t44 * t18 - t6 * t47 * t50) * t23 * t31 -
             t56 * t21 * t23 * t31 + 2.0 * t62 * t65 * t2 * t30 +
             2.0 * t62 * t70 - t61 * t73 * t15 -
             t61 * (t36 - t39 + t25 * t44 * t18 - t25 * t47 * t50) * t15 -
             2.0 * t61 * t65 * t70 - t61 * t2 * t40 * t28 * t29;
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  } else if (deriv == "par1u1") {
    // ported from VineCopula hfuncderiv2.c diff2hfunc_par_v (family 4 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      double theta = par(0);
      double t1 = std::log(u1);
      double t2 = std::pow(-t1, 1.0 * theta);
      double t3 = std::log(u2);
      double t4 = std::pow(-t3, 1.0 * theta);
      double t5 = t2 + t4;
      double t6 = 1 / theta;
      double t7 = std::pow(t5, 1.0 * t6);
      double t8 = t2 * t2;
      double t10 = u1 * u1;
      double t11 = 1 / t10;
      double t12 = t1 * t1;
      double t13 = 1 / t12;
      double t14 = t11 * t13;
      double t15 = t7 * t8 * t14;
      double t16 = 1 / t5;
      double t17 = theta * theta;
      double t19 = std::log(t5);
      double t20 = 1 / t17 * t19;
      double t21 = std::log(-t1);
      double t23 = std::log(-t3);
      double t25 = t2 * t21 + t4 * t23;
      double t28 = -t20 + t6 * t25 * t16;
      double t30 = std::exp(-t7);
      double t31 = t6 - 1.0;
      double t32 = std::pow(t5, 1.0 * t31);
      double t33 = t30 * t32;
      double t37 = 1 / u1;
      double t38 = 1 / t1;
      double t39 = t37 * t38;
      double t41 = t6 * t2 * t39 * t16;
      double t42 = t2 * theta;
      double t46 = t2 * t37 * t38;
      double t47 = t42 * t39 * t21 + t46;
      double t50 = t5 * t5;
      double t51 = 1 / t50;
      double t57 = t32 * t2;
      double t60 = t7 * t7;
      double t64 = t13 * t16;
      double t67 = t7 * t28;
      double t75 = t42 * t14;
      double t77 = t67 * t30;
      double t83 = t16 * t30;
      double t84 = t31 * t25;
      double t86 = -t20 + t84 * t16;
      double t91 = t33 * t31 * t8;
      double t92 = theta * t11;
      double t104 = t33 * t86;
      double t106 = t2 * t11;
      double t109 = t106 * t13;
      double t117 = t33 * t2;
      double t122 = t21 * t11;
      return t15 * t16 * t28 * t33 +
             t7 * (-t41 + t6 * t47 * t16 - t25 * t51 * t46) * t30 * t57 * t39 -
             t60 * t28 * t8 * t11 * t64 * t33 +
             t67 * t33 * t31 * t8 * theta * t14 * t16 + t67 * t33 * t75 -
             t77 * t57 * t11 * t38 - t77 * t57 * t14 + t15 * t83 * t32 * t86 -
             t91 * t92 * t64 * t86 -
             t33 * (-t41 + t31 * t47 * t16 - t84 * t51 * t42 * t39) * t46 -
             t104 * t75 + t104 * t106 * t38 + t104 * t109 +
             t15 * t83 * t32 * t21 - t91 * t92 * t64 * t21 -
             t117 * t92 * t13 * t21 - t33 * t109 + t117 * t122 * t38 +
             t117 * t122 * t13;
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  } else if (deriv == "u1u1") {
    // ported from VineCopula hfuncderiv2.c diff2hfunc_v (family 4 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      double theta = par(0);
      double t1 = std::log(u1);
      double t2 = std::pow(-t1, 1.0 * theta);
      double t3 = std::log(u2);
      double t4 = std::pow(-t3, 1.0 * theta);
      double t5 = t2 + t4;
      double t6 = 1 / theta;
      double t7 = std::pow(t5, 1.0 * t6);
      double t8 = std::exp(-t7);
      double t9 = t6 - 1.0;
      double t10 = std::pow(t5, 1.0 * t9);
      double t11 = t8 * t10;
      double t12 = u1 * u1;
      double t14 = 1 / t12 / u1;
      double t15 = t2 * t14;
      double t20 = t1 * t1;
      double t21 = 1 / t20;
      double t26 = 1 / t20 / t1;
      double t30 = t7 * t7;
      double t31 = t2 * t2;
      double t32 = t31 * t2;
      double t34 = t5 * t5;
      double t36 = 1 / t34;
      double t37 = t26 * t36;
      double t38 = t37 * t11;
      double t40 = t11 * t2;
      double t41 = theta * t14;
      double t48 = t7 * t31;
      double t49 = t48 * t14;
      double t50 = 1 / t5;
      double t51 = t21 * t50;
      double t55 = t26 * t50;
      double t59 = t7 * t32;
      double t62 = t14 * t26;
      double t65 = t10 * theta;
      double t69 = t59 * t62;
      double t70 = t36 * t8;
      double t79 = t11 * t9 * t31;
      double t86 = t9 * t9;
      double t18 = theta * theta;
      double t16 = t18 * t14;
      double t13 = t16 * t37;
      return -2.0 * t11 * t15 / t1 - 3.0 * t11 * t15 * t21 -
             2.0 * t11 * t15 * t26 - t30 * t32 * t14 * t38 +
             3.0 * t40 * t41 * t21 + 3.0 * t40 * t41 * t26 -
             3.0 * t49 * t51 * t11 - 3.0 * t49 * t55 * t11 + t59 * t14 * t38 +
             3.0 * t48 * t62 * t50 * t8 * t65 - t69 * t70 * t65 +
             2.0 * t69 * t70 * t10 * t9 * theta + 3.0 * t79 * t41 * t51 +
             3.0 * t79 * t41 * t55 - t11 * t86 * t32 * t13 -
             3.0 * t79 * t16 * t55 + t11 * t9 * t32 * t13 - t40 * t16 * t26;
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  }
  throw std::runtime_error("unexpected derivative selector: " + deriv);
}

inline Eigen::VectorXd
GumbelBicop::logpdf_deriv_raw(const Eigen::MatrixXd& u,
                              const Eigen::MatrixXd& parameters,
                              const std::string& deriv)
{
  // fused single-pass form kept for performance: the base
  // logpdf_deriv*_raw would compose this as (pdf_deriv)/pdf, recomputing the
  // shared pdf temporaries 2-3x; this leaf is on the scores/Hessian hot paths.
  if (deriv == "par1") {
    // ported from VineCopula logderiv.c difflPDF (family 4 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      double theta = par(0);
      double t1 = std::log(u1);
      double t2 = std::pow(-t1, 1.0 * theta);
      double t3 = std::log(u2);
      double t4 = std::pow(-t3, 1.0 * theta);
      double t5 = t2 + t4;
      double t6 = 1 / theta;
      double t7 = std::pow(t5, 1.0 * t6);
      double t8 = theta * theta;
      double t10 = std::log(t5);
      double t11 = 1 / t8 * t10;
      double t12 = std::log(-t1);
      double t14 = std::log(-t3);
      double t16 = t2 * t12 + t4 * t14;
      double t18 = 1 / t5;
      double t20 = -t11 + t6 * t16 * t18;
      double t22 = std::exp(-t7);
      double t23 = -1.0 + t6;
      double t24 = std::pow(t5, 2.0 * t23);
      double t25 = t22 * t24;
      double t27 = t1 * t3;
      double t28 = theta - 1.0;
      double t29 = std::pow(t27, 1.0 * t28);
      double t30 = std::pow(t5, -1.0 * t6);
      double t31 = t28 * t30;
      double t32 = 1.0 + t31;
      double t34 = 1 / u1;
      double t35 = 1 / u2;
      double t36 = t34 * t35;
      double t37 = t29 * t32 * t36;
      double t45 = t25 * t29;
      double t46 = std::log(t27);
      return (-t7 * t20 * t25 * t37 +
              t25 * (-2.0 * t11 + 2.0 * t23 * t16 * t18) * t37 +
              t45 * t46 * t32 * t36 + t45 * (t30 - t31 * t20) * t34 * t35) /
             t22 / t24 / t29 / t32 * u1 * u2;
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  }
  return AbstractBicop::logpdf_deriv_raw(u, parameters, deriv);
}

inline Eigen::VectorXd
GumbelBicop::logpdf_deriv2_raw(const Eigen::MatrixXd& u,
                               const Eigen::MatrixXd& parameters,
                               const std::string& deriv)
{
  // fused single-pass form kept for performance: the base
  // logpdf_deriv*_raw would compose this as (pdf_deriv)/pdf, recomputing the
  // shared pdf temporaries 2-3x; this leaf is on the scores/Hessian hot paths.
  if (deriv == "par1par1") {
    // ported from VineCopula logderiv.c diff2lPDF (family 4 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      double theta = par(0);
      double t3 = std::log(u1);
      double t4 = std::pow(-t3, 1.0 * theta);
      double t5 = std::log(u2);
      double t6 = std::pow(-t5, 1.0 * theta);
      double t7 = t4 + t6;
      double t8 = 1 / theta;
      double t9 = std::pow(t7, 1.0 * t8);
      double t10 = theta * theta;
      double t11 = 1 / t10;
      double t12 = std::log(t7);
      double t13 = t11 * t12;
      double t14 = std::log(-t3);
      double t16 = std::log(-t5);
      double t18 = t4 * t14 + t6 * t16;
      double t20 = 1 / t7;
      double t22 = -t13 + t8 * t18 * t20;
      double t23 = t22 * t22;
      double t25 = std::exp(-t9);
      double t27 = t25 / u1;
      double t29 = 1 / u2;
      double t30 = -1.0 + t8;
      double t31 = std::pow(t7, 2.0 * t30);
      double t32 = t29 * t31;
      double t33 = t3 * t5;
      double t34 = theta - 1.0;
      double t35 = std::pow(t33, 1.0 * t34);
      double t36 = std::pow(t7, -1.0 * t8);
      double t37 = t34 * t36;
      double t38 = 1.0 + t37;
      double t39 = t35 * t38;
      double t40 = t32 * t39;
      double t44 = 1 / t10 / theta * t12;
      double t47 = t11 * t18 * t20;
      double t49 = t14 * t14;
      double t51 = t16 * t16;
      double t53 = t4 * t49 + t6 * t51;
      double t56 = t18 * t18;
      double t58 = t7 * t7;
      double t59 = 1 / t58;
      double t61 = 2.0 * t44 - 2.0 * t47 + t8 * t53 * t20 - t8 * t56 * t59;
      double t65 = t9 * t9;
      double t70 = t9 * t22 * t27;
      double t74 = -2.0 * t13 + 2.0 * t30 * t18 * t20;
      double t75 = t74 * t35;
      double t80 = std::log(t33);
      double t87 = t36 - t37 * t22;
      double t88 = t35 * t87;
      double t92 = t27 * t29;
      double t93 = t74 * t74;
      double t108 = t80 * t38;
      double t112 = t31 * t74;
      double t116 = t31 * t35;
      double t117 = t80 * t80;
      double t132 = -t9 * t23 * t27 * t40 - t9 * t61 * t27 * t40 +
                    t65 * t23 * t27 * t40 - 2.0 * t70 * t32 * t75 * t38 -
                    2.0 * t70 * t32 * t35 * t80 * t38 - 2.0 * t70 * t32 * t88 +
                    t92 * t31 * t93 * t39 +
                    t92 * t31 *
                      (4.0 * t44 - 4.0 * t47 + 2.0 * t30 * t53 * t20 -
                       2.0 * t30 * t56 * t59) *
                      t39 +
                    2.0 * t27 * t32 * t75 * t108 + 2.0 * t92 * t112 * t88 +
                    t92 * t116 * t117 * t38 + 2.0 * t92 * t116 * t80 * t87 +
                    t92 * t116 * (-2.0 * t36 * t22 + t37 * t23 - t37 * t61);
      double t133 = 1 / t25;
      double t136 = 1 / t31;
      double t138 = 1 / t35;
      double t139 = 1 / t38;
      double t153 =
        (-t70 * t40 + t92 * t112 * t39 + t92 * t116 * t108 + t92 * t116 * t87) *
        t133 * u1 * u2;
      double t154 = t136 * t138;
      double t165 = t38 * t38;
      return t132 * t133 * u1 * u2 * t136 * t138 * t139 +
             t153 * t154 * t139 * t9 * t22 - t153 * t154 * t139 * t74 -
             t153 * t154 * t139 * t80 - t153 * t154 / t165 * t87;
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  }
  return AbstractBicop::logpdf_deriv2_raw(u, parameters, deriv);
}

inline Eigen::MatrixXd
GumbelBicop::tau_to_parameters(const double& tau)
{
  auto par = Eigen::VectorXd::Constant(1, 1.0 / (1 - std::fabs(tau)));
  return par.cwiseMax(parameters_lower_bounds_)
    .cwiseMin(parameters_upper_bounds_);
}

inline double
GumbelBicop::parameters_to_tau(const Eigen::MatrixXd& parameters)
{
  return (parameters(0) - 1) / parameters(0);
}

inline Eigen::MatrixXd
GumbelBicop::parameters_to_taildep(const Eigen::MatrixXd& parameters)
{
  Eigen::MatrixXd taildep = Eigen::MatrixXd::Zero(2, 2);
  taildep(1, 1) = 2 - std::pow(2.0, 1.0 / parameters(0)); // upper tail dep.
  return taildep;
}

inline Eigen::VectorXd
GumbelBicop::get_start_parameters(const double tau)
{
  Eigen::VectorXd par = tau_to_parameters(tau);
  par = par.cwiseMax(parameters_lower_bounds_);
  par = par.cwiseMin(parameters_upper_bounds_);
  return par;
}
}

// This is copy&paste from the VineCopula package
inline double
qcondgum(const double& q, const double& u, const double& de)
{
  double a, p, z1, z2, con, de1, dif;
  double mxdif;
  int iter;

  p = 1 - q;
  z1 = -std::log(u);
  con = std::log(1. - p) - z1 + (1. - de) * std::log(z1);
  de1 = de - 1.;
  a = std::pow(2. * std::pow(z1, de), 1. / (de));
  mxdif = 1;
  iter = 0;
  dif = .1; // needed in case first step leads to NaN
  while ((mxdif > 1.e-6) && (iter < 20)) {
    double g = a + de1 * std::log(a) + con;
    double gp = 1. + de1 / a;
    if ((std::isnan)(g) || (std::isnan)(gp) || (std::isnan)(g / gp)) {
      // added for de>50
      dif /= -2.;
    } else {
      dif = g / gp;
    }
    a -= dif;
    iter++;
    int it = 0;
    while ((a <= z1) && (it < 20)) {
      dif /= 2.;
      a += dif;
      ++it;
    }
    mxdif = fabs(dif);
  }
  z2 = std::pow(std::pow(a, de) - std::pow(z1, de), 1. / (de));
  return (std::exp(-z2));
}
