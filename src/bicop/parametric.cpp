// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/bicop/parametric.hpp>
#include <vinecopulib/misc/tools_optimization.hpp>
#include <vinecopulib/misc/tools_stl.hpp>
#include <vinecopulib/misc/tools_stats.hpp>
#include <iostream>

namespace vinecopulib
{
    Eigen::MatrixXd ParBicop::get_parameters() const
    {
        return parameters_;
    }

    Eigen::MatrixXd ParBicop::get_parameters_lower_bounds() const
    {
        return parameters_lower_bounds_;
    }

    Eigen::MatrixXd ParBicop::get_parameters_upper_bounds() const
    {
        return parameters_upper_bounds_;
    }

    void ParBicop::set_parameters(const Eigen::MatrixXd& parameters)
    {
        check_parameters(parameters);
        parameters_ = parameters;
    }

    void ParBicop::flip()
    {
        // Most parametric families can be flipped by changing the rotation.
        // This is done in Bicop::flip() directly. All other families need to
        // override this method.
    }
    
    // calculate number of parameters
    double ParBicop::calculate_npars() {
        // indepence copula has no parameters
        if (family_ == BicopFamily::indep) {
            return 0.0;
        }
        // otherwise, return length of parameter vector
        return (double) parameters_.size();
    }

    // fit
    void ParBicop::fit(const Eigen::Matrix<double, Eigen::Dynamic, 2>& data,
                       std::string method, double)
    {
        if (family_ != BicopFamily::indep) {
            using namespace tools_optimization;

            std::vector<std::string> methods = {"itau", "mle"};
            if (!tools_stl::is_member(method, methods)) {
                throw std::runtime_error("Method not implemented.");
            }

            int npars = (int) calculate_npars();
            if (method == "itau") {
                npars = npars - 1;
                if ((npars > 0) & (family_ != BicopFamily::student)) {
                    throw std::runtime_error(
                        "itau method is not available for this family.");
                }
            }

            auto temp_data = data;
            double tau = tools_stats::pairwise_tau(temp_data);
            auto newpar = get_start_parameters(tau);
            if (npars > 0) {
                // Create optimizer
                Optimizer optimizer(npars);

                // Set bounds and starting values
                auto lb = get_parameters_lower_bounds();
                auto ub = get_parameters_upper_bounds();
                
                // ensure that starting values are sufficiently separated from 
                // bounds
                double sign = 1.0;
                if (tau < 0) sign = -1.0;
                if (std::abs(tau) < 0.01) {
                    tau = 0.01 * sign;
                } else if (std::abs(tau) > 0.9) {
                    tau = 0.9 * sign;
                }
                auto initial_parameters = get_start_parameters(tau);
                
                ParBicopOptData my_data = {temp_data, this, initial_parameters(0), 0};
                if (method == "itau") {
                      lb.resize(1, 1);
                      lb(0) = get_parameters_lower_bounds()(1);
                      ub.resize(1, 1);
                      ub(0) = get_parameters_upper_bounds()(1);
                      initial_parameters = newpar.tail(1);
                      if (family_ == BicopFamily::student) {
                          // the df parameter doesn't need to be estimated as
                          // accurately
                          ub(0) = 15;
                          optimizer = Optimizer(npars, 1e-2, 5e-1, 1e-4, 1e-4, 100);
                          optimizer.set_objective(pmle_objective, &my_data);
                      }
                      optimizer.set_objective(pmle_objective, &my_data);
                } else {
                      optimizer.set_objective(mle_objective, &my_data);
                }

                // optimize and store result
                optimizer.set_bounds(lb, ub);
                auto optimized_parameters = optimizer.optimize(initial_parameters);
                if (method == "itau") {
                    newpar(1) = optimized_parameters(0);
                } else {
                    newpar = optimized_parameters;
                }
            }
            
            // set the new parameters
            set_parameters(newpar);
        }
    }

    //! Sanity checks
    //! @{
    void ParBicop::check_parameters(const Eigen::MatrixXd& parameters)
    {
        check_parameters_size(parameters);
        check_parameters_lower(parameters);
        check_parameters_upper(parameters);
    }


    void ParBicop::check_parameters_size(const Eigen::MatrixXd& parameters)
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


    void ParBicop::check_parameters_lower(const Eigen::MatrixXd& parameters)
    {
        if (parameters_lower_bounds_.size() > 0) {
            std::stringstream message;
            if ((parameters.array() < parameters_lower_bounds_.array()).any()) {
                message <<
                        "parameters exceed lower bound " <<
                        "for " << get_family_name() << " copula; " << std::endl <<
                        "bound:" << std::endl << parameters_lower_bounds_ << std::endl <<
                        "actual:" << std::endl << parameters << std::endl;
                throw std::runtime_error(message.str().c_str());
            }
        }
    }

    void ParBicop::check_parameters_upper(const Eigen::MatrixXd& parameters)
    {
        if (parameters_upper_bounds_.size() > 0) {
            std::stringstream message;
            if ((parameters.array() > parameters_upper_bounds_.array()).any()) {
                message <<
                        "parameters exceed upper bound " <<
                        "for " << get_family_name() << " copula; " << std::endl <<
                        "bound:" << std::endl << parameters_upper_bounds_ << std::endl <<
                        "actual:" << std::endl << parameters << std::endl;
                throw std::runtime_error(message.str().c_str());
            }
        }
    }

    //! @}
    
}

/*void remove_row(Eigen::MatrixXd& matrix, unsigned int to_remove)
{
    unsigned int n = matrix.rows()-1;
    unsigned int m = matrix.cols();

    if(to_remove < numRows )
        matrix.block(to_remove,0,numRows-to_remove,numCols) = matrix.block(to_remove+1,0,numRows-to_remove,numCols);

    matrix.conservativeResize(numRows,numCols);
}*/
