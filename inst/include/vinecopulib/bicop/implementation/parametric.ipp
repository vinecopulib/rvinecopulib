// Copyright © 2018 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/misc/tools_optimization.hpp>
#include <vinecopulib/misc/tools_stl.hpp>
#include <vinecopulib/misc/tools_stats.hpp>
#include <iostream>
#include <vinecopulib/misc/tools_bobyqa.hpp>
#include <wdm/eigen.hpp>

namespace vinecopulib {
inline Eigen::MatrixXd ParBicop::get_parameters() const
{
    return parameters_;
}

inline Eigen::MatrixXd ParBicop::get_parameters_lower_bounds() const
{
    return parameters_lower_bounds_;
}

inline Eigen::MatrixXd ParBicop::get_parameters_upper_bounds() const
{
    return parameters_upper_bounds_;
}

inline void ParBicop::set_parameters(const Eigen::MatrixXd &parameters)
{
    check_parameters(parameters);
    parameters_ = parameters;
}

inline void ParBicop::flip()
{
    // Most parametric families can be flipped by changing the rotation.
    // This is done in Bicop::flip() directly. All other families need to
    // override this method.
}

// calculate number of parameters
inline double ParBicop::calculate_npars()
{
    // indepence copula has no parameters
    if (family_ == BicopFamily::indep) {
        return 0.0;
    }
    // otherwise, return length of parameter vector
    return static_cast<double>(parameters_.size());
}

// fit
inline void ParBicop::fit(const Eigen::Matrix<double, Eigen::Dynamic, 2> &data,
                          std::string method, 
                          double, 
                          const Eigen::VectorXd& weights)
{
    if (family_ == BicopFamily::indep) {
        // almost nothing to do
        set_loglik(0.0);
        return;
    }
    
    using namespace tools_optimization;
    std::vector <std::string> methods = {"itau", "mle"};
    if (!tools_stl::is_member(method, methods)) {
        throw std::runtime_error("Method not implemented.");
    }

    int npars = static_cast<int>(calculate_npars());
    if (method == "itau") {
        npars = npars - 1;
        if ((npars > 0) & (family_ != BicopFamily::student)) {
            throw std::runtime_error(
                "itau method is not available for this family.");
        }
    }

    double tau = wdm::wdm(data, "tau", weights)(0, 1);
    auto newpar = get_start_parameters(tau);
    
    if (npars > 0) {
        // Set bounds and starting values
        auto lb = get_parameters_lower_bounds();
        auto ub = get_parameters_upper_bounds();

        // ensure that starting values are sufficiently separated from bounds
        double sign = 1.0;
        if (tau < 0)
            sign = -1.0;
        if (std::abs(tau) < 0.01) {
            tau = 0.01 * sign;
        } else if (std::abs(tau) > 0.9) {
            tau = 0.9 * sign;
        }
        auto initial_parameters = get_start_parameters(tau);
        ParBicopOptData my_data = {
            data,
            this,
            initial_parameters(0), 
            0, 
            weights, 
        0};
        std::function<double(void *, long, const double *)> objective =
            mle_objective;
        if (method == "itau") {
            lb.resize(1, 1);
            lb(0) = get_parameters_lower_bounds()(1);
            ub.resize(1, 1);
            ub(0) = get_parameters_upper_bounds()(1);
            initial_parameters = newpar.tail(1);
            if (family_ == BicopFamily::student) {
                // the df parameter doesn't need to be estimated as accurately
                ub(0) = 15;
            }
            objective = pmle_objective;
        }
            
        // refine search interval for Brent algorithm        
        if (tools_stl::is_member(family_, bicop_families::one_par)) {
            auto lb2 = lb;
            auto ub2 = ub;
            if (tools_stl::is_member(family_, bicop_families::rotationless)) {
                lb = tau_to_parameters(std::max(std::fabs(tau) - 0.1, 1e-10));
                ub = tau_to_parameters(std::min(std::fabs(tau) + 0.1, 0.95));
            } else {
                lb = tau_to_parameters(std::max(tau - 0.1, -0.99));
                ub = tau_to_parameters(std::min(tau + 0.1, 0.99));                    
            }
            // make sure that parameter bounds are respected
            lb = lb2.cwiseMax(lb);
            ub = ub2.cwiseMin(ub);
        }
            
        // create optimizer
        Optimizer optimizer(npars, lb, ub);

        // optimize and store result
        auto optimized_parameters = optimizer.optimize(initial_parameters,
                                                       objective, 
                                                       &my_data);
                                                       
        // check if fit is reasonable, otherwise increase search interval 
        // and refit
        if (tools_stl::is_member(family_, bicop_families::one_par)) {
            if (my_data.objective_min > 0.1) {
                // -loglik should always be negative!
                lb = get_parameters_lower_bounds();
                ub = get_parameters_upper_bounds();
                optimizer = Optimizer(npars, lb, ub);
                optimized_parameters = optimizer.optimize(initial_parameters,
                                                          objective, 
                                                          &my_data);
            }
        }
        if (method == "itau") {
            newpar(1) = optimized_parameters(0);
        } else {
            newpar = optimized_parameters;
        }
        // set the new parameters
        set_loglik(-my_data.objective_min);
    }
    set_parameters(newpar);
    
    if (npars == 0) {
        // loglik has not been computed
        if (weights.size() > 0) {
            set_loglik((pdf(data).array().log() * weights.array()).sum());
        } else {
            set_loglik(pdf(data).array().log().sum());
        }
    }
}

//! Sanity checks
//! @{
inline void ParBicop::check_parameters(const Eigen::MatrixXd &parameters)
{
    check_parameters_size(parameters);
    check_parameters_lower(parameters);
    check_parameters_upper(parameters);
}


inline void ParBicop::check_parameters_size(const Eigen::MatrixXd &parameters)
{
    if (parameters.size() != parameters_.size()) {
        if (parameters.rows() != parameters_.rows()) {
            std::stringstream message;
            message <<
                    "parameters have has wrong number of rows " <<
                    "for " << get_family_name() << " copula; " <<
                    "expected: " << parameters_.rows() << ", " <<
                    "actual: " << parameters.rows() << std::endl;
            throw std::runtime_error(message.str().c_str());
        }
        if (parameters.cols() != parameters_.cols()) {
            std::stringstream message;
            message <<
                    "parameters have wrong number of columns " <<
                    "for " << get_family_name() << " copula; " <<
                    "expected: " << parameters_.cols() << ", " <<
                    "actual: " << parameters.cols() << std::endl;
            throw std::runtime_error(message.str().c_str());
        }
    }
}


inline void
ParBicop::check_parameters_lower(const Eigen::MatrixXd &parameters)
{
    if (parameters_lower_bounds_.size() > 0) {
        std::stringstream message;
        if ((parameters.array() < parameters_lower_bounds_.array()).any()) {
            message <<
                    "parameters exceed lower bound " <<
                    "for " << get_family_name() << " copula; " << std::endl <<
                    "bound:" << std::endl << parameters_lower_bounds_
                    << std::endl <<
                    "actual:" << std::endl << parameters << std::endl;
            throw std::runtime_error(message.str().c_str());
        }
    }
}

inline void
ParBicop::check_parameters_upper(const Eigen::MatrixXd &parameters)
{
    if (parameters_upper_bounds_.size() > 0) {
        std::stringstream message;
        if ((parameters.array() > parameters_upper_bounds_.array()).any()) {
            message <<
                    "parameters exceed upper bound " <<
                    "for " << get_family_name() << " copula; " << std::endl <<
                    "bound:" << std::endl << parameters_upper_bounds_
                    << std::endl <<
                    "actual:" << std::endl << parameters << std::endl;
            throw std::runtime_error(message.str().c_str());
        }
    }
}

//! @}

}
