// Copyright © 2016-2023 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <boost/math/distributions.hpp>
#include <memory>
#include <set>
#include <unsupported/Eigen/SpecialFunctions>
#include <vinecopulib/misc/tools_eigen.hpp>

namespace vinecopulib {

namespace tools_stats {

//! @brief Density function of the Standard normal distribution.
//!
//! @param x Evaluation points.
//!
//! @return An \f$ n \times d \f$ matrix of evaluated densities.
inline Eigen::MatrixXd
dnorm(const Eigen::MatrixXd& x)
{
  static constexpr double inv_sqrt_2pi = 0.39894228040143270286;
  return inv_sqrt_2pi * (-0.5 * x.array().square()).exp();
}

//! @brief Distribution function of the Standard normal distribution.
//!
//! @param x Evaluation points.
//!
//! @return An \f$ n \times d \f$ matrix of evaluated probabilities.
inline Eigen::MatrixXd
pnorm(const Eigen::MatrixXd& x)
{
  static const double sqrt2 = std::sqrt(2.0);
  return 0.5 * (1.0 + (x.array() / sqrt2).erf());
}

//! @brief Quantile function of the Standard normal distribution.
//!
//! @param x Evaluation points.
//!
//! @return An \f$ n \times d \f$ matrix of evaluated quantiles.
inline Eigen::MatrixXd
qnorm(const Eigen::MatrixXd& x)
{
  return x.array().ndtri();
}

//! @brief Density function of the Student t distribution.
//!
//! @param x Evaluation points.
//! @param nu Degrees of freedom parameter.
//!
//! @return An \f$ n \times d \f$ matrix of evaluated densities.
inline Eigen::MatrixXd
dt(const Eigen::MatrixXd& x, double nu)
{
  boost::math::students_t dist(nu);
  auto f = [&dist](double y) { return boost::math::pdf(dist, y); };
  return tools_eigen::unaryExpr_or_nan(x, f);
}

//! @brief Distribution function of the Student t distribution.
//!
//! @param x Evaluation points.
//! @param nu Degrees of freedom parameter.
//!
//! @return An \f$ n \times d \f$ matrix of evaluated probabilities.
inline Eigen::MatrixXd
pt(const Eigen::MatrixXd& x, double nu)
{
  boost::math::students_t dist(nu);
  auto f = [&dist](double y) { return boost::math::cdf(dist, y); };
  return tools_eigen::unaryExpr_or_nan(x, f);
}

//! @brief Quantile function of the Student t distribution.
//!
//! @param x Evaluation points.
//! @param nu Degrees of freedom parameter.
//!
//! @return An \f$ n \times d \f$ matrix of evaluated quantiles.
inline Eigen::MatrixXd
qt(const Eigen::MatrixXd& x, double nu)
{
  boost::math::students_t dist(nu);
  auto f = [&dist](double y) { return boost::math::quantile(dist, y); };
  return tools_eigen::unaryExpr_or_nan(x, f);
}

Eigen::MatrixXd
simulate_uniform(const size_t& n,
                 const size_t& d,
                 bool qrng = false,
                 std::vector<int> seeds = std::vector<int>());

Eigen::MatrixXd
simulate_normal(const size_t& n,
                const size_t& d,
                bool qrng = false,
                std::vector<int> seeds = std::vector<int>());

Eigen::VectorXd
to_pseudo_obs_1d(Eigen::VectorXd x,
                 const std::string& ties_method = "average",
                 const Eigen::VectorXd& weights = Eigen::VectorXd(),
                 std::vector<int> seeds = std::vector<int>());

Eigen::MatrixXd
to_pseudo_obs(Eigen::MatrixXd x,
              const std::string& ties_method = "average",
              const Eigen::VectorXd& weights = Eigen::VectorXd(),
              std::vector<int> seeds = std::vector<int>());

// Covers the unit hypercube with boxes and assigns each sample to a box.
// Used internally for recovering the latent sample of a discrete copula.
class BoxCovering
{
public:
  explicit BoxCovering(const Eigen::MatrixXd& u, uint16_t K = 40);
  std::vector<size_t> get_box_indices(const Eigen::VectorXd& lower,
                                      const Eigen::VectorXd& upper) const;
  void swap_sample(size_t i, const Eigen::VectorXd& new_sample);

private:
  struct Box
  {
  public:
    Box(const std::vector<double>& lower, const std::vector<double>& upper);
    std::vector<double> lower_;
    std::vector<double> upper_;
    std::set<size_t> indices_;
  };

  Eigen::MatrixXd u_;
  size_t n_;
  uint16_t K_;
  std::vector<std::vector<std::unique_ptr<Box>>> boxes_;
};

Eigen::MatrixXd
find_latent_sample(const Eigen::MatrixXd& u, double b, size_t niter = 3);

double
pairwise_mcor(const Eigen::MatrixXd& x,
              const Eigen::VectorXd& weights = Eigen::VectorXd());

Eigen::MatrixXd
dependence_matrix(const Eigen::MatrixXd& x, const std::string& measure);

Eigen::MatrixXd
ghalton(const size_t& n,
        const size_t& d,
        const std::vector<int>& seeds = std::vector<int>());

Eigen::MatrixXd
sobol(const size_t& n,
      const size_t& d,
      const std::vector<int>& seeds = std::vector<int>());

Eigen::VectorXd
pbvt(const Eigen::MatrixXd& z, int nu, double rho);

Eigen::VectorXd
pbvnorm(const Eigen::MatrixXd& z, double rho);
}
}

#include <vinecopulib/misc/implementation/tools_stats.ipp>
