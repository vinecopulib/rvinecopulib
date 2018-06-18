// Copyright © 2018 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <vinecopulib/misc/tools_bobyqa.hpp>
#include <vinecopulib/misc/tools_eigen.hpp>
#include <vinecopulib/bicop/parametric.hpp>

namespace vinecopulib {

namespace tools_optimization {

//! @brief A helper struct for (profile) maximum likelihood estimation
typedef struct
{
    const Eigen::Matrix<double, Eigen::Dynamic, 2> &U; //!< the data.
    vinecopulib::ParBicop *bicop; //!< a pointer to the bivariate copula to optimize.
    double par0;  //!< main dependence parameter.
    unsigned int objective_calls; //!< number of evaluations of the objective.
    Eigen::VectorXd weights; //!< weights for the observations.
    double objective_min;  //!< final value of the objective function 
} ParBicopOptData;

//! @brief A class for the controls to Bobyqa
class BobyqaControls
{
public:

    BobyqaControls();

    BobyqaControls(double initial_trust_region,
                   double final_trust_region,
                   int maxeval);

    double get_initial_trust_region();

    double get_final_trust_region();

    int get_maxeval();

private:
    double initial_trust_region_; //! Initial trust region
    double final_trust_region_; //! Final trust region
    int maxeval_; //! Maximal number of evaluations of the objective

    //! Sanity checks
    //! @{
    void check_parameters(double initial_trust_region,
                          double final_trust_region,
                          int maxeval);
    //! @}
};

//! @brief A class for optimization (wrapping Bobyqa).
class Optimizer
{
public:
    Optimizer(unsigned int n_parameters,
              const Eigen::MatrixXd &lower_bounds,
              const Eigen::MatrixXd &upper_bounds);

    void set_controls(double initial_trust_region,
                      double final_trust_region,
                      int maxeval);

    Eigen::VectorXd optimize(Eigen::VectorXd initial_parameters,
                             std::function<double(void *, long,
                                                  const double *)> objective,
                             void *f_data);

private:
    unsigned int n_parameters_;
    BobyqaControls controls_;
    Eigen::MatrixXd lb_;
    Eigen::MatrixXd ub_;
};
}

}

#include <vinecopulib/misc/implementation/tools_optimization.ipp>
