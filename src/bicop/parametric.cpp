// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#include <vinecopulib/bicop/parametric.hpp>
#include <vinecopulib/misc/tools_optimization.hpp>
#include <vinecopulib/misc/tools_stl.hpp>
#include <vinecopulib/misc/tools_stats.hpp>

namespace vinecopulib
{
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
            double tau = tools_stats::pairwise_ktau(temp_data);
            auto newpar = get_start_parameters(tau);
            if (npars > 0) {
                // Create optimizer
                Optimizer optimizer(npars);

                // Set bounds and starting values
                auto lb = get_parameters_lower_bounds();
                auto ub = get_parameters_upper_bounds();
                auto initial_parameters = newpar;
                ParBicopOptData my_data = {temp_data, this, newpar(0), 0};
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
    
}

/*void remove_row(Eigen::MatrixXd& matrix, unsigned int to_remove)
{
    unsigned int n = matrix.rows()-1;
    unsigned int m = matrix.cols();

    if(to_remove < numRows )
        matrix.block(to_remove,0,numRows-to_remove,numCols) = matrix.block(to_remove+1,0,numRows-to_remove,numCols);

    matrix.conservativeResize(numRows,numCols);
}*/
