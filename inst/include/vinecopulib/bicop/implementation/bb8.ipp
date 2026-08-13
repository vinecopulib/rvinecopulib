// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/misc/tools_eigen.hpp>
#include <vinecopulib/misc/tools_integration.hpp>

namespace vinecopulib {
inline Bb8Bicop::Bb8Bicop()
{
  family_ = BicopFamily::bb8;
  parameters_ = Eigen::VectorXd(2);
  parameters_lower_bounds_ = Eigen::VectorXd(2);
  parameters_upper_bounds_ = Eigen::VectorXd(2);
  parameters_ << 1, 1;
  parameters_lower_bounds_ << 1, 1e-4;
  parameters_upper_bounds_ << 8, 1;
}

inline double
Bb8Bicop::generator(const double& u,
                    const Eigen::Ref<const Eigen::VectorXd>& parameters)
{
  double theta = parameters(0);
  double delta = parameters(1);
  double res = (1 - std::pow(1 - delta * u, theta));
  return -std::log(res / (1 - std::pow(1 - delta, theta)));
}

inline double
Bb8Bicop::generator_inv(const double& u,
                        const Eigen::Ref<const Eigen::VectorXd>& parameters)
{
  double theta = parameters(0);
  double delta = parameters(1);
  double res = std::exp(-u) * (std::pow(1 - delta, theta) - 1);
  return (1 - std::pow(1 + res, 1 / theta)) / delta;
}

inline double
Bb8Bicop::generator_derivative(
  const double& u,
  const Eigen::Ref<const Eigen::VectorXd>& parameters)
{
  double theta = parameters(0);
  double delta = parameters(1);
  double res = delta * theta * std::pow(1 - delta * u, theta - 1);
  return -res / (1 - std::pow(1 - delta * u, theta));
}

inline Eigen::VectorXd
Bb8Bicop::pdf_raw(const Eigen::MatrixXd& u, const Eigen::MatrixXd& parameters)
{
  auto f = [](const double& u1,
              const double& u2,
              const Eigen::Ref<const Eigen::VectorXd>& par) {
    double theta = par(0);
    double delta = par(1);
    double t10 = 1.0 - delta;
    double t16 = 1.0 / theta;
    double t38 = 2.0 * theta;
    double t39 = std::pow(t10, t38);
    double t59 = std::pow(t10, 3.0 * theta);
    double t2 = 1.0 - delta * u1;
    double t3 = std::pow(t2, theta);
    double t11 = std::pow(t10, theta);
    double t12 = 1.0 - t11;
    double t33 = theta * t3;
    double t49 = std::pow(t2, t38);
    double t69 = t12 * t12;
    double t6 = 1.0 - delta * u2;
    double t7 = std::pow(t6, theta);
    double t25 = t3 * t7;
    double t26 = t11 - t7 - t3 + t25;
    double t29 = std::pow(-t26 / t12, t16);
    double t44 = std::pow(t6, t38);
    double t45 = t3 * t44;
    double t50 = t49 * t7;
    double t54 = t49 * t44;
    double t62 = -2.0 * t25 * t11 + t25 - t33 * t7 + 3.0 * t33 * t7 * t11 -
                 3.0 * t33 * t7 * t39 + t25 * t39 + 2.0 * t45 * t11 -
                 t45 * t39 + 2.0 * t50 * t11 - t50 * t39 - 2.0 * t54 * t11 +
                 t54 * t39 + t54 - t50 - t45 + t33 * t7 * t59;
    double t67 = t26 * t26;

    return -delta * t29 * t62 / t6 / t2 / t67 / t69;
  };
  return tools_eigen::binaryExpr_or_nan(u, parameters, f);
}

inline double
Bb8Bicop::parameters_to_tau(const Eigen::MatrixXd& parameters)
{
  double theta = parameters(0);
  double delta = parameters(1);
  auto f = [theta, delta](const double t) {
    double tmp = std::pow(1 - t * delta, theta);
    double res = std::log((tmp - 1) / (std::pow(1 - delta, theta) - 1));
    return res * (1 - t * delta - std::pow(1 - t * delta, 1 - theta));
  };
  return 1 - 4 / (delta * theta) * tools_integration::integrate_zero_to_one(f);
}

inline Eigen::MatrixXd
Bb8Bicop::parameters_to_taildep(const Eigen::MatrixXd& par)
{
  double theta = par(0);
  double delta = par(1);
  Eigen::MatrixXd taildep = Eigen::MatrixXd::Zero(2, 2);
  // BB8 only has upper tail dependence in the limiting case delta = 1
  // (where it reduces to the Joe copula), and none otherwise.
  if (delta >= 1.0) {
    taildep(1, 1) = 2 - std::pow(2.0, 1.0 / theta);
  }
  return taildep;
}

inline Eigen::MatrixXd
Bb8Bicop::tau_to_parameters(const double& tau)
{
  return no_tau_to_parameters(tau);
}
}
