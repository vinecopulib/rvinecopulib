// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <boost/math/distributions.hpp>
#include <random>

#include <vinecopulib/misc/tools_eigen.hpp>

namespace vinecopulib {

namespace tools_stats
{
    //! @name Standard normal distribution
    //! 
    //! @param x evaluation points.
    //! @{

    //! evaluates the density function.
    inline Eigen::MatrixXd dnorm(const Eigen::MatrixXd& x)
    {
        boost::math::normal dist;
        return x.unaryExpr([&dist](double y) {return boost::math::pdf(dist, y);});
    };

    //! evaluates the distribution function.
    inline Eigen::MatrixXd pnorm(const Eigen::MatrixXd& x)
    {
        boost::math::normal dist;
        return x.unaryExpr([&dist](double y) {return boost::math::cdf(dist, y);});
    };
    
    //! evaluates the quantile function.
    inline Eigen::MatrixXd qnorm(const Eigen::MatrixXd& x)
    {
        boost::math::normal dist;
        return x.unaryExpr([&dist](double y) {return boost::math::quantile(dist, y);});
    };
    //! @}
    
    //! @name Student t distribution
    //! The mean is assumed to be zero.
    //! @param x evaluation points.
    //! @param nu degrees of freedom parameter.
    //! @{
    
    //! evaluates the density function.
    inline Eigen::MatrixXd dt(const Eigen::MatrixXd& x, double nu)
    {
        boost::math::students_t dist(nu);
        return x.unaryExpr([&dist](double y) {return boost::math::pdf(dist, y);});
    };

    //! evaluates the distribution function.
    inline Eigen::MatrixXd pt(const Eigen::MatrixXd& x, double nu)
    {
        boost::math::students_t dist(nu);
        return x.unaryExpr([&dist](double y) {return boost::math::cdf(dist, y);});
    };
    
    //! evaluates the  quantile function.
    inline Eigen::MatrixXd qt(const Eigen::MatrixXd& x, double nu)
    {
        boost::math::students_t dist(nu);
        return x.unaryExpr([&dist](double y) {return boost::math::quantile(dist, y);});
    };
    //! @}

    Eigen::MatrixXd simulate_uniform(size_t n, size_t d);
    Eigen::VectorXd to_pseudo_obs_1d(Eigen::VectorXd x,
                                     std::string ties_method = "average");
    Eigen::MatrixXd to_pseudo_obs(Eigen::MatrixXd x,
                                  std::string ties_method = "average");

    double pairwise_tau(Eigen::Matrix<double, Eigen::Dynamic, 2>& u);
    double pairwise_cor(const Eigen::Matrix<double, Eigen::Dynamic, 2>& z);
    double pairwise_rho(Eigen::Matrix<double, Eigen::Dynamic, 2> z);
    double pairwise_hoeffd(Eigen::Matrix<double, Eigen::Dynamic, 2> x);
    Eigen::MatrixXd dependence_matrix(const Eigen::MatrixXd& x, 
                                      const std::string& measure);
    
    Eigen::MatrixXd ghalton(size_t n, size_t d);
}

}
