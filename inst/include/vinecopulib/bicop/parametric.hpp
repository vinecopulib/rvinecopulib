// Copyright © 2018 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <vinecopulib/bicop/abstract.hpp>

namespace vinecopulib {
// Pre-declaration to allow ParBicop to befriend the two functions
namespace tools_optimization {
// the objective function for maximum likelihood estimation
double mle_objective(void *f_data, long n, const double *x);

// the objective function for profile maximum likelihood estimation
double pmle_objective(void *f_data, long n, const double *x);
}

//! @brief An abstract class for parametric copula families
//!
//! This class is used in the implementation underlying the Bicop class.
//! Users should not use AbstractBicop or derived classes directly, but
//! always work with the Bicop interface.
class ParBicop : public AbstractBicop
{
    friend double tools_optimization::mle_objective(void *f_data, long n,
                                                    const double *x);

    friend double tools_optimization::pmle_objective(void *f_data, long n,
                                                     const double *x);

protected:
    // Getters and setters
    Eigen::MatrixXd get_parameters() const;

    Eigen::MatrixXd get_parameters_lower_bounds() const;

    Eigen::MatrixXd get_parameters_upper_bounds() const;

    void set_parameters(const Eigen::MatrixXd &parameters);

    void flip();

    // Data members
    Eigen::MatrixXd parameters_;
    Eigen::MatrixXd parameters_lower_bounds_;
    Eigen::MatrixXd parameters_upper_bounds_;

private:
    void fit(const Eigen::Matrix<double, Eigen::Dynamic, 2> &data,
             std::string method, 
             double, 
             const Eigen::VectorXd& weights);

    double calculate_npars();

    virtual Eigen::VectorXd get_start_parameters(const double tau) = 0;

    void check_parameters(const Eigen::MatrixXd &parameters);

    void check_parameters_size(const Eigen::MatrixXd &parameters);

    void check_parameters_upper(const Eigen::MatrixXd &parameters);

    void check_parameters_lower(const Eigen::MatrixXd &parameters);
};
}

#include <vinecopulib/bicop/implementation/parametric.ipp>
