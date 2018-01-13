// Copyright © 2018 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <boost/math/distributions.hpp>
#include <random>
#include <vinecopulib/misc/tools_eigen.hpp>

namespace vinecopulib {

namespace tools_stats {


//! Density function of the Standard normal distribution
//!
//! @param x evaluation points.
//!
//! @return An \f$ n \times d \f$ matrix of evaluated densities.
inline Eigen::MatrixXd dnorm(const Eigen::MatrixXd &x)
{
    boost::math::normal dist;
    auto f = [&dist](double y) { return boost::math::pdf(dist, y); };
    return tools_eigen::unaryExpr_or_nan(x, f);
}

//! Distribution function of the Standard normal distribution
//!
//! @param x evaluation points.
//!
//! @return An \f$ n \times d \f$ matrix of evaluated probabilities.
inline Eigen::MatrixXd pnorm(const Eigen::MatrixXd &x)
{
    boost::math::normal dist;
    auto f = [&dist](double y) { return boost::math::cdf(dist, y); };
    return tools_eigen::unaryExpr_or_nan(x, f);
}

//! Quantile function of the Standard normal distribution
//!
//! @param x evaluation points.
//!
//! @return An \f$ n \times d \f$ matrix of evaluated quantiles.
inline Eigen::MatrixXd qnorm(const Eigen::MatrixXd &x)
{
    boost::math::normal dist;
    auto f = [&dist](double y) { return boost::math::quantile(dist, y); };
    return tools_eigen::unaryExpr_or_nan(x, f);
}

//! Density function of the Student t distribution
//!
//! @param x evaluation points.
//! @param nu degrees of freedom parameter.
//!
//! @return An \f$ n \times d \f$ matrix of evaluated densities.
inline Eigen::MatrixXd dt(const Eigen::MatrixXd &x, double nu)
{
    boost::math::students_t dist(nu);
    auto f = [&dist](double y) { return boost::math::pdf(dist, y); };
    return tools_eigen::unaryExpr_or_nan(x, f);
}

//! Distribution function of the Student t distribution
//!
//! @param x evaluation points.
//! @param nu degrees of freedom parameter.
//!
//! @return An \f$ n \times d \f$ matrix of evaluated probabilities.
inline Eigen::MatrixXd pt(const Eigen::MatrixXd &x, double nu)
{
    boost::math::students_t dist(nu);
    auto f = [&dist](double y) { return boost::math::cdf(dist, y); };
    return tools_eigen::unaryExpr_or_nan(x, f);
}

//! Quantile function of the Student t distribution
//!
//! @param x evaluation points.
//! @param nu degrees of freedom parameter.
//!
//! @return An \f$ n \times d \f$ matrix of evaluated quantiles.
inline Eigen::MatrixXd qt(const Eigen::MatrixXd &x, double nu)
{
    boost::math::students_t dist(nu);
    auto f = [&dist](double y) { return boost::math::quantile(dist, y); };
    return tools_eigen::unaryExpr_or_nan(x, f);
}

Eigen::MatrixXd simulate_uniform(size_t n, size_t d);

Eigen::VectorXd to_pseudo_obs_1d(Eigen::VectorXd x,
                                 std::string ties_method = "average");

Eigen::MatrixXd to_pseudo_obs(Eigen::MatrixXd x,
                              std::string ties_method = "average");

double pairwise_tau(const Eigen::Matrix<double, Eigen::Dynamic, 2>& x);

double pairwise_cor(const Eigen::Matrix<double, Eigen::Dynamic, 2>& x);

double pairwise_rho(const Eigen::Matrix<double, Eigen::Dynamic, 2>& x);

double pairwise_hoeffd(const Eigen::Matrix<double, Eigen::Dynamic, 2>& x);

Eigen::MatrixXd dependence_matrix(const Eigen::MatrixXd &x,
                                  const std::string &measure);

Eigen::MatrixXd ghalton(size_t n, size_t d);

Eigen::VectorXd pbvt(const Eigen::Matrix<double, Eigen::Dynamic, 2> &z,
                     int nu, double rho);

Eigen::VectorXd pbvnorm(const Eigen::Matrix<double, Eigen::Dynamic, 2> &z,
                        double rho);
}

}

#include <vinecopulib/misc/implementation/tools_stats.ipp>
