// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/misc/tools_stats.hpp>

namespace vinecopulib {

//! Utilities for numerical optimization (based on Bobyqa)
namespace tools_optimization {

//! creates an Optimizer using the default controls, see BobyqaControls.
//!
//! @param n_parameters Number of parameters to optimize
//! @param lower_bounds
//! @param upper_bounds
//! @param objective The optimizer's objective function
inline Optimizer::Optimizer(unsigned int n_parameters,
                            const Eigen::MatrixXd &lower_bounds,
                            const Eigen::MatrixXd &upper_bounds)
{
    if (n_parameters < 1) {
        throw std::runtime_error("n_parameters should be larger than 0.");
    }
    n_parameters_ = n_parameters;
    controls_ = BobyqaControls();
    lb_ = lower_bounds;
    ub_ = upper_bounds;
}

//! set the optimizer's controls.
//!
//! @param initial_trust_region initial trust region.
//! @param final_trust_region final trust region.
//! @param maxeval maximal number of evaluations of the objective.
inline void Optimizer::set_controls(double initial_trust_region,
                                    double final_trust_region,
                                    int maxeval)
{
    controls_ = BobyqaControls(initial_trust_region,
                               final_trust_region,
                               maxeval);
}

//! Create controls using the default contructor
//!
//! The defaults are
//! ```
//! initial_trust_region_ = 1e-3;
//! final_trust_region_ = 1e3;
//! maxeval_ = 1000;
//! ```
inline BobyqaControls::BobyqaControls()
{
    initial_trust_region_ = 1e-3;
    final_trust_region_ = 1e3;
    maxeval_ = 1000;
}

//! Create controls by passing the arguments
//!
//! @param initial_trust_region initial trust region.
//! @param final_trust_region final trust region.
//! @param maxeval maximal number of evaluations of the objective.
inline BobyqaControls::BobyqaControls(double initial_trust_region,
                                      double final_trust_region,
                                      int maxeval)
{
    check_parameters(initial_trust_region, final_trust_region, maxeval);
    initial_trust_region_ = initial_trust_region;
    final_trust_region_ = final_trust_region;
    maxeval_ = maxeval;
}

inline void BobyqaControls::check_parameters(double initial_trust_region,
                                             double final_trust_region,
                                             int maxeval)
{
    if (initial_trust_region <= 0) {
        throw std::runtime_error(
            "initial_trust_region should be larger than 0");
    }
    if (final_trust_region <= 0) {
        throw std::runtime_error("final_trust_region should be larger than 0");
    }
    if (maxeval <= 0) {
        throw std::runtime_error("maxeval should be larger than 0");
    }
}

//! @name Getters and setters
//! @{

//! @return the initial trust region.
inline double
BobyqaControls::get_initial_trust_region()
{
    return initial_trust_region_;
};

//! @return the final trust region.
inline double
BobyqaControls::get_final_trust_region()
{
    return final_trust_region_;
};

//! @return the maximal number of evaluations of the objective.
inline int BobyqaControls::get_maxeval()
{
    return maxeval_;
};

//! @}

//! @name (Pseudo-) maximum likelihood estimation
//! @param f_data a pointer to a ParBicopOptData object.
//! @param n number of parameters.
//! @param x the parameters.
//! @{

//! evaluates the objective function for maximum likelihood estimation.
inline double mle_objective(void *f_data, long n, const double *x)
{
    ParBicopOptData *newdata = (ParBicopOptData *) f_data;
    ++newdata->objective_calls;
    Eigen::Map<const Eigen::VectorXd> par(&x[0], n);
    newdata->bicop->set_parameters(par);
    return (-1) * newdata->bicop->pdf(newdata->U).array().log().sum();
}

//! evaluates the objective function for profile maximum likelihood
//! estimation.
inline double pmle_objective(void *f_data, long n, const double *x)
{
    ParBicopOptData *newdata = (ParBicopOptData *) f_data;
    ++newdata->objective_calls;
    Eigen::VectorXd par = Eigen::VectorXd::Ones(n + 1);
    par(0) = newdata->par0;
    for (long i = 0; i < n; ++i) {
        par(i + 1) = x[i];
    }
    newdata->bicop->set_parameters(par);
    double nll = newdata->bicop->pdf(newdata->U).array().log().sum();
    nll *= -1;
    return nll;
}

//! @}

//! solve the optimization problem.
//!
//! @param initial_parameters of starting values for the optimization
//!     algorithm.
//! @return the optimal parameters.
inline Eigen::VectorXd Optimizer::optimize(Eigen::VectorXd initial_parameters,
                                           tools_bobyqa::BobyqaClosureFunction objective,
                                           void *data)
{
    if (initial_parameters.size() != n_parameters_) {
        throw std::runtime_error("The size of x should be n_parameters_.");
    }

    //const int number_interpolation_conditions = (n_parameters_ + 1) *
    //        (n_parameters_ + 2)/2;
    int number_interpolation_conditions = n_parameters_ + 3;
    std::size_t ws_size = (number_interpolation_conditions + 5) *
                          (number_interpolation_conditions
                           + n_parameters_) +
                          3 * n_parameters_ * (n_parameters_ + 5) / 2;
    double *working_space = new double[ws_size];

    double *lb = new double[n_parameters_];
    double *ub = new double[n_parameters_];
    double eps = 1e-6;
    Eigen::VectorXd::Map(lb, n_parameters_) = lb_.array() + eps;
    Eigen::VectorXd::Map(ub, n_parameters_) = ub_.array() - eps;

    double *x = new double[n_parameters_];
    Eigen::VectorXd::Map(x, n_parameters_) = initial_parameters;
    Eigen::VectorXd optimized_parameters(n_parameters_);
    std::string err_msg = "";
    try {
        auto f = [&](long n, const double *x) -> double {
            return objective(data, n, x);
        };
        tools_bobyqa::impl(f, n_parameters_,
                           number_interpolation_conditions,
                           x, lb, ub,
                           controls_.get_initial_trust_region(),
                           controls_.get_final_trust_region(),
                           controls_.get_maxeval(),
                           working_space);

        // convert result to VectorXd
        for (size_t i = 0; i < n_parameters_; i++) {
            optimized_parameters(i) = x[i];
        }
    } catch (std::invalid_argument err) {
        err_msg = std::string("Invalid arguments. ") + err.what();
    } catch (std::bad_alloc err) {
        err_msg = std::string("Ran out of memory. ") + err.what();
    } catch (std::runtime_error err) {
        err_msg = std::string("Generic failure. ") + err.what();
    } catch (...) {
        // do nothing for other errors (results are fine)
    }

    // delete dynamically allocated objects
    delete[] x;
    delete[] lb;
    delete[] ub;
    delete[] working_space;

    // throw error if optimization failed
    if (err_msg != "") {
        throw std::runtime_error(err_msg);
    }

    return optimized_parameters;
}
}

}
