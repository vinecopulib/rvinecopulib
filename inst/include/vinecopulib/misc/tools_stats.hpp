// Copyright © 2016-2023 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <boost/random.hpp>
#include <random>
#include <boost/math/distributions.hpp>
#include <vinecopulib/misc/tools_eigen.hpp>

namespace vinecopulib {

namespace tools_stats {

class SeedState {
public:
  SeedState(std::vector<int> seeds = {}) {
    if (seeds.size() == 0) {
      // no seeds provided, seed randomly
      std::random_device rd{};
      seeds = std::vector<int>(5);
      std::generate(
        seeds.begin(), seeds.end(), [&]() { return static_cast<int>(rd()); });
    }
    boost::random::seed_seq seq(seeds.begin(), seeds.end());
    engine_ = boost::mt19937(seq);
  }

  std::vector<int> next() {
    std::vector<int> seeds(5);
    std::generate(
      seeds.begin(), seeds.end(), [&]() { return static_cast<int>(engine_()); });
    return seeds;
  }

private:
  boost::mt19937 engine_;
};


//! @brief Density function of the Standard normal distribution.
//!
//! @param x Evaluation points.
//!
//! @return An \f$ n \times d \f$ matrix of evaluated densities.
inline Eigen::MatrixXd
dnorm(const Eigen::MatrixXd& x)
{
  boost::math::normal dist;
  auto f = [&dist](double y) { return boost::math::pdf(dist, y); };
  return tools_eigen::unaryExpr_or_nan(x, f);
}

//! @brief Distribution function of the Standard normal distribution.
//!
//! @param x Evaluation points.
//!
//! @return An \f$ n \times d \f$ matrix of evaluated probabilities.
inline Eigen::MatrixXd
pnorm(const Eigen::MatrixXd& x)
{
  boost::math::normal dist;
  auto f = [&dist](double y) { return boost::math::cdf(dist, y); };
  return tools_eigen::unaryExpr_or_nan(x, f);
}

//! @brief Quantile function of the Standard normal distribution.
//!
//! @param x Evaluation points.
//!
//! @return An \f$ n \times d \f$ matrix of evaluated quantiles.
inline Eigen::MatrixXd
qnorm(const Eigen::MatrixXd& x)
{
  boost::math::normal dist;
  auto f = [&dist](double y) { return boost::math::quantile(dist, y); };
  return tools_eigen::unaryExpr_or_nan(x, f);
}

inline Eigen::MatrixXd
safe_qnorm(const Eigen::MatrixXd& x)
{
  boost::math::normal dist;
  auto f = [&dist](double y) {
    if (y <= 0) {
      return -std::numeric_limits<double>::infinity();
    } else if (y >= 1) {
      return std::numeric_limits<double>::infinity();
    } else {
      return boost::math::quantile(dist, y);
    }
  };

  return tools_eigen::unaryExpr_or_nan(x, f);
}

inline Eigen::MatrixXd
safe_pnorm(const Eigen::MatrixXd& x)
{
  boost::math::normal dist;
  auto f = [&dist](double y) {
    if (y >= std::numeric_limits<double>::max()) {
      return 1.0;
    } else if (y <= -std::numeric_limits<double>::max()) {
      return 0.0;
    } else {
      return boost::math::cdf(dist, y);
    }
  };

  return tools_eigen::unaryExpr_or_nan(x, f);
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
                std::vector<int> seeds = std::vector<int>());

Eigen::VectorXd
to_pseudo_obs_1d(Eigen::VectorXd x, const std::string& ties_method = "average");

Eigen::MatrixXd
to_pseudo_obs(Eigen::MatrixXd x, const std::string& ties_method = "average");

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
