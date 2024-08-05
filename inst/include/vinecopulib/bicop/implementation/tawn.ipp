// Copyright © 2016-2023 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/misc/tools_eigen.hpp>
#include <vinecopulib/misc/tools_integration.hpp>

namespace vinecopulib {
inline TawnBicop::TawnBicop()
{
  family_ = BicopFamily::tawn;
  parameters_ = Eigen::VectorXd(3);
  parameters_lower_bounds_ = Eigen::VectorXd(3);
  parameters_upper_bounds_ = Eigen::VectorXd(3);
  parameters_ << 0.5, 0.5, 6;
  parameters_lower_bounds_ << 1e-10, 1e-10, 1 + 1e-10;
  parameters_upper_bounds_ << 1 - 1e-10, 1 - 1e-10, 60;
}

inline double
TawnBicop::pickands(const double& t)
{
  double psi1 = double(this->parameters_(0));
  double psi2 = double(this->parameters_(1));
  double theta = double(this->parameters_(2));

  double temp = std::pow(psi2 * t, theta) + 
                  std::pow(psi1 * (1 - t), theta);
  return (1 - psi1) * (1 - t) + (1 - psi2) * t + 
            std::pow(temp, 1 / theta);
}

inline double 
TawnBicop::pickands_derivative(const double& t)
{
  double psi1 = double(this->parameters_(0));
  double psi2 = double(this->parameters_(1));
  double theta = double(this->parameters_(2));

  double temp = std::pow(psi2 * t, theta) + 
                  std::pow(psi1 * (1 - t), theta);
  double temp2 = psi2 * std::pow(psi2 * t, theta - 1) - 
                  psi1 * std::pow(psi1 * (1 - t), theta - 1);
  return psi1 - psi2 + std::pow(temp, 1 / theta - 1) * temp2;
}

inline double
TawnBicop::pickands_derivative2(const double& t)
{
  double psi1 = double(this->parameters_(0));
  double psi2 = double(this->parameters_(1));
  double theta = double(this->parameters_(2));

  double temp = std::pow(psi2 * t, theta) + 
                  std::pow(psi1 * (1 - t), theta);
  double temp2 = psi2 * std::pow(psi2 * t, theta - 1) - 
                  psi1 * std::pow(psi1 * (1 - t), theta - 1);
  double temp3 = std::pow(psi2, 2) * std::pow(psi2 * t, theta - 2) + 
                  std::pow(psi1, 2) * std::pow(psi1 * (1 - t), theta - 2);
  return (1 - theta) * std::pow(temp, 1 / theta - 2) * std::pow(temp2, 2) + 
          std::pow(temp, 1 / theta - 1) * (theta - 1) * temp3;
}

inline double
TawnBicop::parameters_to_tau(const Eigen::MatrixXd& par)
{
  double psi1 = par(0);
  double psi2 = par(1);
  double theta = par(2);

  auto f = [psi1, psi2, theta](const double t) {
    double t2 = std::pow(psi2 * t, theta) + 
                  std::pow(psi1 * (1 - t), theta);
    double p = (1 - psi1) * (1 - t) + (1 - psi2) * t + 
                std::pow(t2, 1 / theta);
    double t3 = psi2 * std::pow(psi2 * t, theta - 1) - 
                  psi1 * std::pow(psi1 * (1 - t), theta - 1);
    double t4 = std::pow(psi2, 2) * std::pow(psi2 * t, theta - 2) + 
                  std::pow(psi1, 2) * std::pow(psi1 * (1 - t), theta - 2);
    double p2 = (1 - theta) * std::pow(t2, 1 / theta - 2) * std::pow(t3, 2) + 
                  std::pow(t2, 1 / theta - 1) * (theta - 1) * t4;

    return t * (1 - t) * p2 / p;
  };
  return tools_integration::integrate_zero_to_one(f);
}

inline Eigen::MatrixXd
TawnBicop::tau_to_parameters(const double& tau)
{
  return no_tau_to_parameters(tau);
}
}