// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#include "bicop/gaussian.hpp"
#include "misc/tools_stats.hpp"
#ifndef M_PI
#define M_PI       3.14159265358979323846
#endif

namespace vinecopulib
{
    GaussianBicop::GaussianBicop()
    {
        family_ = BicopFamily::gaussian;
        parameters_ = Eigen::VectorXd(1);
        parameters_lower_bounds_ = Eigen::VectorXd(1);
        parameters_upper_bounds_ = Eigen::VectorXd(1);
        parameters_ << 0;
        parameters_lower_bounds_ << -1;
        parameters_upper_bounds_ << 1;
    }

    Eigen::VectorXd GaussianBicop::pdf(
        const Eigen::Matrix<double, Eigen::Dynamic, 2>& u
    )
    {
        // Inverse Cholesky of the correlation matrix
        double rho = double(this->parameters_(0));
        Eigen::Matrix2d L;
        L(0,0) = 1;
        L(1,1) = 1/sqrt(1.0-pow(rho,2.0));
        L(0,1) = -rho*L(1,1);

        // Compute copula density
        Eigen::VectorXd f = Eigen::VectorXd::Ones(u.rows());
        Eigen::Matrix<double, Eigen::Dynamic, 2> tmp = tools_stats::qnorm(u);
        f = f.cwiseQuotient(tools_stats::dnorm(tmp).rowwise().prod());
        tmp = tmp*L;
        f = f.cwiseProduct(tools_stats::dnorm(tmp).rowwise().prod());
        return f / sqrt(1.0-pow(rho,2.0));
    }

    Eigen::VectorXd GaussianBicop::hfunc1(
        const Eigen::Matrix<double, Eigen::Dynamic, 2>& u
    )
    {
        double rho = double(this->parameters_(0));
        Eigen::VectorXd h = Eigen::VectorXd::Zero(u.rows());
        Eigen::Matrix<double, Eigen::Dynamic, 2> tmp = tools_stats::qnorm(u);
        h = (tmp.col(1) - rho * tmp.col(0)) / sqrt(1.0 - pow(rho, 2.0));
        return tools_stats::pnorm(h);
    }

    Eigen::VectorXd GaussianBicop::hinv1(
        const Eigen::Matrix<double, Eigen::Dynamic, 2>& u
    )
    {
        double rho = double(this->parameters_(0));
        Eigen::VectorXd hinv = Eigen::VectorXd::Zero(u.rows());
        Eigen::Matrix<double, Eigen::Dynamic, 2> tmp = tools_stats::qnorm(u);
        hinv = tmp.col(1) * sqrt(1.0 - pow(rho, 2.0)) + rho * tmp.col(0);
        return tools_stats::pnorm(hinv);
    }

    Eigen::VectorXd GaussianBicop::get_start_parameters(const double tau)
    {
        return tau_to_parameters(tau);
    }

    Eigen::MatrixXd GaussianBicop::tau_to_parameters(const double& tau)
    {
        Eigen::VectorXd parameters = this->parameters_;
        parameters(0) = sin(tau * M_PI / 2);
        return parameters;
    }
}
