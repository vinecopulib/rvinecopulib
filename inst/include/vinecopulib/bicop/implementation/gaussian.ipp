// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <cmath>
#include <vinecopulib/misc/tools_eigen.hpp>
#include <vinecopulib/misc/tools_stats.hpp>

namespace vinecopulib {
inline GaussianBicop::GaussianBicop()
{
  family_ = BicopFamily::gaussian;
  parameters_ = Eigen::VectorXd(1);
  parameters_lower_bounds_ = Eigen::VectorXd(1);
  parameters_upper_bounds_ = Eigen::VectorXd(1);
  parameters_ << 0;
  parameters_lower_bounds_ << -1;
  parameters_upper_bounds_ << 1;
}

inline Eigen::VectorXd
GaussianBicop::pdf_raw(const Eigen::MatrixXd& u,
                       const Eigen::MatrixXd& parameters)
{
  const Eigen::Index n = u.rows();
  if (parameters.rows() == 1) {
    // broadcast: scalar inverse-Cholesky form (bit-identical to the historical
    // state-based path)
    double rho = parameters(0, 0);
    Eigen::Matrix2d L;
    L(0, 0) = 1;
    L(1, 1) = 1 / sqrt(1.0 - pow(rho, 2.0));
    L(0, 1) = -rho * L(1, 1);
    L(1, 0) = 0;

    Eigen::VectorXd f = Eigen::VectorXd::Ones(n);
    Eigen::MatrixXd tmp = tools_stats::qnorm(u);
    f = f.cwiseQuotient(tools_stats::dnorm(tmp).rowwise().prod());
    tmp = tmp * L;
    f = f.cwiseProduct(tools_stats::dnorm(tmp).rowwise().prod());
    return f / sqrt(1.0 - pow(rho, 2.0));
  }

  Eigen::ArrayXd rho =
    tools_eigen::parameter_as_vector(parameters, 0, n).array();
  Eigen::ArrayXd s = (1.0 - rho.square()).sqrt(); // sqrt(1 - rho^2)

  Eigen::MatrixXd z = tools_stats::qnorm(u);
  Eigen::MatrixXd zL(n, 2);
  zL.col(0) = z.col(0);
  zL.col(1) = (z.col(1).array() - rho * z.col(0).array()) / s;

  Eigen::VectorXd f = tools_stats::dnorm(zL).rowwise().prod();
  f = f.cwiseQuotient(tools_stats::dnorm(z).rowwise().prod());
  f = (f.array() / s).matrix();
  return f;
}

inline Eigen::VectorXd
GaussianBicop::cdf(const Eigen::MatrixXd& u, const Eigen::MatrixXd& parameters)
{
  Eigen::MatrixXd z = tools_stats::qnorm(u);
  if (parameters.rows() == 1) {
    return tools_stats::pbvnorm(z, parameters(0, 0));
  }
  const Eigen::Index n = u.rows();
  Eigen::VectorXd p(n);
  for (Eigen::Index i = 0; i < n; ++i) {
    p(i) = tools_stats::pbvnorm(z.row(i), parameters(i, 0))(0);
  }
  return p;
}

inline Eigen::VectorXd
GaussianBicop::hfunc1_raw(const Eigen::MatrixXd& u,
                          const Eigen::MatrixXd& parameters)
{
  Eigen::MatrixXd z = tools_stats::qnorm(u);
  if (parameters.rows() == 1) {
    // broadcast: scalar form (computes sqrt(1 - rho^2) once; bit-identical to
    // the historical state-based path)
    double rho = parameters(0, 0);
    Eigen::VectorXd h = (z.col(1) - rho * z.col(0)) / sqrt(1.0 - pow(rho, 2.0));
    return tools_stats::pnorm(h);
  }
  const Eigen::Index n = u.rows();
  Eigen::ArrayXd rho =
    tools_eigen::parameter_as_vector(parameters, 0, n).array();
  Eigen::VectorXd h =
    ((z.col(1).array() - rho * z.col(0).array()) / (1.0 - rho.square()).sqrt())
      .matrix();
  return tools_stats::pnorm(h);
}

inline Eigen::VectorXd
GaussianBicop::hinv1_raw(const Eigen::MatrixXd& u,
                         const Eigen::MatrixXd& parameters)
{
  Eigen::MatrixXd z = tools_stats::qnorm(u);
  if (parameters.rows() == 1) {
    double rho = parameters(0, 0);
    Eigen::VectorXd hinv =
      z.col(1) * sqrt(1.0 - pow(rho, 2.0)) + rho * z.col(0);
    return tools_stats::pnorm(hinv);
  }
  const Eigen::Index n = u.rows();
  Eigen::ArrayXd rho =
    tools_eigen::parameter_as_vector(parameters, 0, n).array();
  Eigen::VectorXd hinv =
    (z.col(1).array() * (1.0 - rho.square()).sqrt() + rho * z.col(0).array())
      .matrix();
  return tools_stats::pnorm(hinv);
}

inline Eigen::VectorXd
GaussianBicop::pdf_deriv_raw(const Eigen::MatrixXd& u,
                             const Eigen::MatrixXd& parameters,
                             const std::string& deriv)
{
  // the copula is exchangeable: route "u2" through a swap of the arguments
  if (tools_deriv::is_u2_only(deriv)) {
    return pdf_deriv_raw(
      tools_eigen::swap_cols(u), parameters, tools_deriv::swap_args(deriv));
  }

  if (deriv == "par1") {
    // ported from VineCopula deriv.c diffPDF (family 1 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      const boost::math::normal_distribution<> norm_dist;
      const double theta = par(0);
      const double t1 = boost::math::quantile(norm_dist, u1);
      const double t2 = boost::math::quantile(norm_dist, u2);
      const double t4 = t1 * t1;
      const double t5 = t2 * t2;
      const double t3 = t4 + t5;
      const double t7 = theta * theta;
      const double t8 = 1.0 - t7;
      const double t9 = 1.0 / t8 / 2.0;
      const double t15 = t7 * t3 - 2.0 * theta * t1 * t2;
      const double t22 = std::exp(-t15 * t9);
      const double t24 = std::sqrt(t8);
      return (-2.0 * (theta * t3 - t1 * t2) * t9 - t15 / (t8 * t8) * theta) *
               t22 / t24 +
             t22 / t24 / t8 * theta;
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  } else if (deriv == "u1") {
    // ported from VineCopula deriv.c diffPDF_u (family 1 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      const boost::math::normal_distribution<> norm_dist;
      const double theta = par(0);
      const double t1 = boost::math::quantile(norm_dist, u1);
      const double t2 = boost::math::quantile(norm_dist, u2);
      const double t3 = theta * theta;
      const double t4 = 1.0 - t3;
      const double t5 = std::sqrt(t4);
      const double t6 = t1 * t1;
      const double t7 = t2 * t2;
      const double t8 = t3 * (t6 + t7) - (2.0 * theta * t1 * t2);
      const double t9 = std::exp(-t8 / t4 / 2.0);
      const double t10 = t9 / t5;
      const double t11 = std::sqrt(2.0 * constant::pi);
      const double t12 = std::exp(-t6 / 2.0);
      const double t13 = t11 / t12;
      return -t10 * (theta * t13 / t4) * (theta * t1 - t2);
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  }
  throw std::runtime_error("unexpected derivative selector: " + deriv);
}

inline Eigen::VectorXd
GaussianBicop::pdf_deriv2_raw(const Eigen::MatrixXd& u,
                              const Eigen::MatrixXd& parameters,
                              const std::string& deriv)
{
  // the copula is exchangeable: route "u2" through a swap of the arguments
  if (tools_deriv::is_u2_only(deriv)) {
    return pdf_deriv2_raw(
      tools_eigen::swap_cols(u), parameters, tools_deriv::swap_args(deriv));
  }

  if (deriv == "par1par1") {
    // ported from VineCopula deriv2.c diff2PDF (family 1 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      const boost::math::normal_distribution<> norm_dist;
      const double theta = par(0);
      const double t6 = boost::math::quantile(norm_dist, u1);
      const double t7 = boost::math::quantile(norm_dist, u2);
      const double t1 = t6 * t7;
      const double t2 = theta * theta;
      const double t3 = 1.0 - t2;
      const double t4 = 4.0 * t3 * t3;
      const double t5 = 1.0 / t4;
      const double t12 = t6 * t6;
      const double t13 = t7 * t7;
      const double t14 = 2.0 * theta * t6 * t7 - t12 - t13;
      const double t21 = t14 * t5;
      const double t26 = 1.0 / t3 / 2.0;
      const double t29 = std::exp(t12 / 2.0 + t13 / 2.0 + t14 * t26);
      const double t31 = std::sqrt(t3);
      const double t32 = 1.0 / t31;
      const double t38 = 2.0 * t1 * t26 + 4.0 * t21 * theta;
      const double t39 = t38 * t38;
      const double t44 = 1.0 / t31 / t3;
      const double t48 = t3 * t3;
      return (16.0 * t1 * t5 * theta + 16.0 * t14 / t4 / t3 * t2 + 4.0 * t21) *
               t29 * t32 +
             t39 * t29 * t32 + 2.0 * t38 * t29 * t44 * theta +
             3.0 * t29 / t31 / t48 * t2 + t29 * t44;
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  } else if (deriv == "par1u1") {
    // ported from VineCopula deriv2.c diff2PDF_par_u (family 1 branch); the
    // preamble computes c (LL) and diffc (diffPDF) inline
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      const boost::math::normal_distribution<> norm_dist;
      const double theta = par(0);
      const double t1 = boost::math::quantile(norm_dist, u1);
      const double t2 = boost::math::quantile(norm_dist, u2);
      const double s3 = t1 * t1 + t2 * t2;
      const double s7 = theta * theta;
      const double s8 = 1.0 - s7;
      const double s9 = 1.0 / s8 / 2.0;
      const double s15 = s7 * s3 - 2.0 * theta * t1 * t2;
      const double c = std::exp(-s15 * s9) / std::sqrt(s8);
      const double diffc =
        (-2.0 * (theta * s3 - t1 * t2) * s9 - s15 / (s8 * s8) * theta) * c +
        c / s8 * theta;
      const double t3 = 1.0 / boost::math::pdf(norm_dist, t1);
      const double t4 = theta * (theta * t1 - t2);
      const double t5 = 1.0 - theta * theta;
      const double t6 = -t4 / t5;
      const double t7 = -2.0 * theta * t1 + t2 + t2 * theta * theta;
      const double t8 = std::pow(-1.0 + theta * theta, 2.0);
      const double t9 = t7 / t8;
      return diffc * t6 * t3 + c * t3 * t9;
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  } else if (deriv == "u1u1") {
    // ported from VineCopula deriv2.c diff2PDF_u (family 1 branch); the
    // preamble computes c (LL) and diffc (diffPDF_u_mod) inline
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      const boost::math::normal_distribution<> norm_dist;
      const double theta = par(0);
      const double t1 = boost::math::quantile(norm_dist, u1);
      const double t2 = boost::math::quantile(norm_dist, u2);
      const double t3 = boost::math::pdf(norm_dist, t1);
      const double t5 = 1.0 - theta * theta;
      const double s8 =
        theta * theta * (t1 * t1 + t2 * t2) - 2.0 * theta * t1 * t2;
      const double c = std::exp(-s8 / t5 / 2.0) / std::sqrt(t5);
      const double diffc = -c * (theta / t3 / t5) * (theta * t1 - t2);
      const double t6 = theta * t1 - t2;
      const double t7 = theta + t6 * t1;
      const double t8 = t7 / t3 / t3;
      return -theta / t5 * (diffc * t6 / t3 + c * t8);
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  } else if (deriv == "u1u2") {
    // ported from VineCopula deriv2.c diff2PDF_u_v (family 1 branch); the
    // preamble computes c (LL) and diffc (diffPDF_v_mod) inline
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      const boost::math::normal_distribution<> norm_dist;
      const double theta = par(0);
      const double t1 = boost::math::quantile(norm_dist, u1);
      const double t2 = boost::math::quantile(norm_dist, u2);
      const double t3 = boost::math::pdf(norm_dist, t1);
      const double t4 = boost::math::pdf(norm_dist, t2);
      const double t5 = 1.0 - theta * theta;
      const double s8 =
        theta * theta * (t1 * t1 + t2 * t2) - 2.0 * theta * t1 * t2;
      const double c = std::exp(-s8 / t5 / 2.0) / std::sqrt(t5);
      const double diffc = -c * (theta / t4 / t5) * (theta * t2 - t1);
      const double t7 = theta * t1 - t2;
      return -theta / t3 / t5 * (diffc * t7 - c / t4);
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  }
  throw std::runtime_error("unexpected derivative selector: " + deriv);
}

inline Eigen::VectorXd
GaussianBicop::hfunc1_deriv_raw(const Eigen::MatrixXd& u,
                                const Eigen::MatrixXd& parameters,
                                const std::string& deriv)
{
  // VineCopula's h-function kernels differentiate h(u|v), which conditions
  // on their second argument v; our hfunc1(u1, u2) conditions on the first,
  // so their u is our u2 and their v is our u1.
  if (deriv == "par1") {
    // ported from VineCopula hfuncderiv.c diffhfunc (family 1 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      const boost::math::normal_distribution<> norm_dist;
      const double theta = par(0);
      const double t1 = boost::math::quantile(norm_dist, u2);
      const double t2 = boost::math::quantile(norm_dist, u1);
      const double t3 = t1 - theta * t2;
      const double t4 = 1.0 - theta * theta;
      const double t5 = std::sqrt(t4);
      const double t6 = t3 / t5;
      const double t7 = boost::math::pdf(norm_dist, t6);
      const double t8 = -1.0 * t2 * t5 + 1.0 * t3 * theta / t5;
      const double t9 = t8 / t4;
      return t7 * t9;
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  } else if (deriv == "u1") {
    // ported from VineCopula hfuncderiv.c diffhfunc_v (family 1 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      const boost::math::normal_distribution<> norm_dist;
      const double theta = par(0);
      const double t1 = boost::math::quantile(norm_dist, u2);
      const double t2 = boost::math::quantile(norm_dist, u1);
      const double t3 = t1 - theta * t2;
      const double t4 = 1.0 - theta * theta;
      const double t5 = std::sqrt(t4);
      const double t6 = t3 / t5;
      const double t7 = boost::math::pdf(norm_dist, t6);
      const double t8 = std::sqrt(2.0 * constant::pi);
      const double t9 = t2 * t2;
      const double t10 = std::exp(-t9 / 2.0);
      return t7 * t8 * (-theta) / t5 / t10;
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  }
  throw std::runtime_error("unexpected derivative selector: " + deriv);
}

inline Eigen::VectorXd
GaussianBicop::hfunc1_deriv2_raw(const Eigen::MatrixXd& u,
                                 const Eigen::MatrixXd& parameters,
                                 const std::string& deriv)
{
  // see hfunc1_deriv_raw for the argument convention (their u/v = our u2/u1)
  if (deriv == "par1par1") {
    // ported from VineCopula hfuncderiv2.c diff2hfunc (family 1 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      const boost::math::normal_distribution<> norm_dist;
      const double theta = par(0);
      const double t1 = boost::math::quantile(norm_dist, u2);
      const double t2 = boost::math::quantile(norm_dist, u1);
      const double t6 = t1 - theta * t2;
      const double t7 = 1.0 - theta * theta;
      const double t8 = std::sqrt(t7);
      const double t9 = t6 / t8;
      const double t10 = boost::math::pdf(norm_dist, t9);
      const double t11 = t6 * t2 / t7;
      const double t12 = t6 * t6 * theta / (t7 * t7);
      const double t13 = t10 * (t11 - t12);
      const double t14 = -t2 * t8 + (t1 - theta * t2) * theta / t8;
      const double t15 = t14 / t7;
      const double t16 = 2.0 * t1 * theta * theta - 3.0 * theta * t2 + t1;
      const double t18 = std::pow(t7, 2.5);
      const double t20 = t16 / t18;
      return t13 * t15 + t10 * t20;
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  } else if (deriv == "par1u1") {
    // ported from VineCopula hfuncderiv2.c diff2hfunc_par_v (family 1 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      const boost::math::normal_distribution<> norm_dist;
      const double theta = par(0);
      const double t1 = boost::math::quantile(norm_dist, u2);
      const double t2 = boost::math::quantile(norm_dist, u1);
      const double t4 = boost::math::pdf(norm_dist, t2);
      const double t5 = 1.0 / t4;
      const double t6 = t1 - theta * t2;
      const double t7 = 1.0 - theta * theta;
      const double t8 = std::sqrt(t7);
      const double t9 = t6 / t8;
      const double t10 = boost::math::pdf(norm_dist, t9);
      const double t11 = t6 * t2 / t7;
      const double t12 = t6 * t6 * theta / (t7 * t7);
      const double t13 = -theta / t8;
      const double t14 = -1.0 / t8;
      const double t15 = theta * theta / std::pow(t7, 1.5);
      return t5 * (t10 * (t11 - t12) * t13 + t10 * (t14 - t15));
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  } else if (deriv == "u1u1") {
    // ported from VineCopula hfuncderiv2.c diff2hfunc_v (family 1 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      const boost::math::normal_distribution<> norm_dist;
      const double theta = par(0);
      const double t1 = boost::math::quantile(norm_dist, u2);
      const double t2 = boost::math::quantile(norm_dist, u1);
      const double t4 = boost::math::pdf(norm_dist, t2);
      const double t5 = 1.0 - theta * theta;
      const double t6 = std::sqrt(t5);
      const double t7 = (t1 - theta * t2) / t6;
      const double t8 = boost::math::pdf(norm_dist, t7);
      const double t9 = t4 * t4;
      return -theta * t8 / t6 / t9 * (t7 / t6 * theta + t2);
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  }
  throw std::runtime_error("unexpected derivative selector: " + deriv);
}

inline Eigen::VectorXd
GaussianBicop::logpdf_deriv_raw(const Eigen::MatrixXd& u,
                                const Eigen::MatrixXd& parameters,
                                const std::string& deriv)
{
  if (deriv == "par1") {
    // ported from VineCopula logderiv.c difflPDF (family 1 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      const boost::math::normal_distribution<> norm_dist;
      const double theta = par(0);
      const double t1 = boost::math::quantile(norm_dist, u1);
      const double t2 = boost::math::quantile(norm_dist, u2);
      const double t3 = theta * theta;
      const double t4 = 1.0 - t3;
      const double t5 = t1 * t1;
      const double t6 = t2 * t2;
      const double t7 = t4 * t4;
      return (theta * t4 - theta * (t5 + t6) + (1.0 + t3) * t1 * t2) / t7;
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  }
  return AbstractBicop::logpdf_deriv_raw(u, parameters, deriv);
}

inline Eigen::VectorXd
GaussianBicop::logpdf_deriv2_raw(const Eigen::MatrixXd& u,
                                 const Eigen::MatrixXd& parameters,
                                 const std::string& deriv)
{
  if (deriv == "par1par1") {
    // ported from VineCopula logderiv.c diff2lPDF (family 1 branch)
    auto f = [](const double& u1,
                const double& u2,
                const Eigen::Ref<const Eigen::VectorXd>& par) {
      const boost::math::normal_distribution<> norm_dist;
      const double theta = par(0);
      const double t5 = boost::math::quantile(norm_dist, u1);
      const double t6 = boost::math::quantile(norm_dist, u2);
      const double t1 = theta * theta;
      const double t3 = t5 * t5;
      const double t4 = t6 * t6;
      const double t9 = 1.0 - t1;
      const double t10 = t9 * t9;
      return (1.0 - 3.0 * t1 - t3 - t4 + 2.0 * theta * t5 * t6) / t10 +
             4.0 * (theta * t9 - theta * (t3 + t4) + (1.0 + t1) * t5 * t6) /
               t10 / t9 * theta;
    };
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  }
  return AbstractBicop::logpdf_deriv2_raw(u, parameters, deriv);
}

inline Eigen::VectorXd
GaussianBicop::get_start_parameters(const double tau)
{
  return tau_to_parameters(tau);
}

inline Eigen::MatrixXd
GaussianBicop::tau_to_parameters(const double& tau)
{
  Eigen::VectorXd parameters = this->parameters_;
  parameters(0) = std::sin(tau * constant::pi / 2);
  return parameters;
}

inline Eigen::MatrixXd
GaussianBicop::parameters_to_taildep(const Eigen::MatrixXd&)
{
  return Eigen::MatrixXd::Zero(2, 2);
}
}
