// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/bicop/gaussian.hpp>
#include <vinecopulib/misc/tools_stats.hpp>
#include <vinecopulib/misc/tools_c.h>
#include <boost/math/constants/constants.hpp>

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
        L(1,0) = 0;

        // Compute copula density
        Eigen::VectorXd f = Eigen::VectorXd::Ones(u.rows());
        Eigen::Matrix<double, Eigen::Dynamic, 2> tmp = tools_stats::qnorm(u);
        f = f.cwiseQuotient(tools_stats::dnorm(tmp).rowwise().prod());
        tmp = tmp*L;
        f = f.cwiseProduct(tools_stats::dnorm(tmp).rowwise().prod());
        return f / sqrt(1.0-pow(rho,2.0));
    }

    Eigen::VectorXd GaussianBicop::cdf(
            const Eigen::Matrix<double, Eigen::Dynamic, 2>& u
    )
    {
        ptrdiff_t n = u.rows();
        Eigen::VectorXd p = Eigen::VectorXd::Ones(u.rows());
        double rho = double(this->parameters_(0));

        double abseps = 0.001, releps = 0, error = 0;
        int d = 2, nu = 0, maxpts = 25000, inform;
        std::vector<double> lower(2), upper(2);
        std::vector<int> infin(2);

        auto v = tools_stats::qnorm(u);
        for (ptrdiff_t i = 0; i < n; i++)
        {
            upper[0] = v(i,0);
            upper[1] = v(i,1);
            mvtdst_(&d, &nu, &lower[0], &upper[0], &infin[0], &rho,
                    &lower[0], &maxpts,&abseps, &releps, &error, &p(i), &inform);
        }
        return p;
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
        parameters(0) = sin(tau * boost::math::constants::pi<double>() / 2);
        return parameters;
    }
}
