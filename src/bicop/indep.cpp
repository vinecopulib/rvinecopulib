// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/bicop/indep.hpp>

namespace vinecopulib
{
    IndepBicop::IndepBicop()
    {
        family_ = BicopFamily::indep;
        parameters_ = Eigen::MatrixXd();
    }

    Eigen::VectorXd IndepBicop::pdf(
        const Eigen::Matrix<double, Eigen::Dynamic, 2>& u
    )
    {
        return Eigen::VectorXd::Ones(u.rows());
    }

    Eigen::VectorXd IndepBicop::cdf(
            const Eigen::Matrix<double, Eigen::Dynamic, 2>& u
    )
    {
        return u.rowwise().prod();
    }

    Eigen::VectorXd IndepBicop::hfunc1(
        const Eigen::Matrix<double, Eigen::Dynamic, 2>& u
    )
    {
        return u.col(1);
    }

    Eigen::VectorXd IndepBicop::hfunc2(
        const Eigen::Matrix<double, Eigen::Dynamic, 2>& u
    )
    {
        return u.col(0);
    }

    Eigen::VectorXd IndepBicop::hinv1(
        const Eigen::Matrix<double, Eigen::Dynamic, 2>& u
    )
    {
        return u.col(1);
    }

    Eigen::VectorXd IndepBicop::hinv2(
        const Eigen::Matrix<double, Eigen::Dynamic, 2>& u
    )
    {
        return u.col(0);
    }

    Eigen::MatrixXd IndepBicop::tau_to_parameters(const double &)
    {
        return Eigen::VectorXd();
    }

    double IndepBicop::parameters_to_tau(const Eigen::MatrixXd&)
    {
        return 0.0;
    }

    Eigen::VectorXd IndepBicop::get_start_parameters(const double tau)
    {
        return tau_to_parameters(tau);
    }

    void IndepBicop::flip()
    {
        // nothing to do because independence copula is radially syemmetric
    }
}
