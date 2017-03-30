// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#include "bicop/student.hpp"
#include "misc/tools_stats.hpp"
#include <cmath>
#ifndef M_PI
#define M_PI       3.14159265358979323846
#endif

namespace vinecopulib
{
    StudentBicop::StudentBicop()
    {
        family_ = BicopFamily::student;
        parameters_ = Eigen::VectorXd(2);
        parameters_lower_bounds_ = Eigen::VectorXd(2);
        parameters_upper_bounds_ = Eigen::VectorXd(2);
        parameters_ << 0, 50;
        parameters_lower_bounds_ << -1, 2;
        parameters_upper_bounds_ << 1, 50;
    }

    Eigen::VectorXd StudentBicop::pdf(
        const Eigen::Matrix<double, Eigen::Dynamic, 2>& u
    )
    {
        double rho = double(this->parameters_(0));
        double nu = double(this->parameters_(1));
        Eigen::VectorXd f = Eigen::VectorXd::Ones(u.rows());
        Eigen::Matrix<double, Eigen::Dynamic, 2> tmp = tools_stats::qt(u, nu);

        f = tmp.col(0).cwiseAbs2() + tmp.col(1).cwiseAbs2() - (2 * rho) * tmp.rowwise().prod();
        f /= nu * (1.0 - pow(rho, 2.0));
        f = f + Eigen::VectorXd::Ones(u.rows());
        f = f.array().pow(-(nu + 2.0) / 2.0);
        f = f.cwiseQuotient(tools_stats::dt(tmp, nu).rowwise().prod());
        //f *= StableGammaDivision((nu + 2.0) / 2.0, nu / 2.0) / (nu * M_PI * sqrt(1.0 - pow(rho, 2.0)));
        f *= boost::math::tgamma_ratio((nu + 2.0) / 2.0, nu / 2.0) / (nu * M_PI * sqrt(1.0 - pow(rho, 2.0)));

        return f;
    }

    Eigen::VectorXd StudentBicop::hfunc1(
        const Eigen::Matrix<double, Eigen::Dynamic, 2>& u
    )
    {
        double rho = double(this->parameters_(0));
        double nu = double(this->parameters_(1));
        Eigen::VectorXd h = Eigen::VectorXd::Ones(u.rows());
        Eigen::Matrix<double, Eigen::Dynamic, 2> tmp = tools_stats::qt(u, nu);
        h = nu * h + tmp.col(0).cwiseAbs2();
        h *= (1.0 - pow(rho, 2)) / (nu + 1.0);
        h = h.cwiseSqrt().cwiseInverse().cwiseProduct(tmp.col(1) - rho * tmp.col(0));
        h = tools_stats::pt(h, nu + 1.0);

        return h;
    }

    Eigen::VectorXd StudentBicop::hinv1(
        const Eigen::Matrix<double, Eigen::Dynamic, 2>& u
    )
    {
        double rho = double(this->parameters_(0));
        double nu = double(this->parameters_(1));
        Eigen::VectorXd hinv = Eigen::VectorXd::Ones(u.rows());
        Eigen::VectorXd tmp = u.col(1);
        Eigen::VectorXd tmp2 = u.col(0);
        tmp = tools_stats::qt(tmp, nu + 1.0);
        tmp2 = tools_stats::qt(tmp2, nu);

        hinv = nu * hinv + tmp2.cwiseAbs2();
        hinv *= (1.0 - pow(rho, 2)) / (nu + 1.0);
        hinv = hinv.cwiseSqrt().cwiseProduct(tmp) + rho * tmp2;
        hinv = tools_stats::pt(hinv, nu );

        return hinv;
    }

    Eigen::VectorXd StudentBicop::get_start_parameters(const double tau)
    {
        Eigen::VectorXd parameters = get_parameters();
        parameters(0) = sin(tau * M_PI / 2);;
        parameters(1) = 5;
        return parameters;
    }

    Eigen::MatrixXd StudentBicop::tau_to_parameters(const double& tau)
    {
        return vinecopulib::no_tau_to_parameters(tau);
    }
}
