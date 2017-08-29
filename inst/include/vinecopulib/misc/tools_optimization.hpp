// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <nlopt.hpp>
#include <vinecopulib/misc/tools_eigen.hpp>
#include <vinecopulib/bicop/parametric.hpp>

namespace vinecopulib {

namespace tools_optimization {

    //! @brief A helper struct for (profile) maximum likelihood estimation
    typedef struct
    {
        const Eigen::Matrix<double, Eigen::Dynamic, 2>& U; //!< the data.
        vinecopulib::ParBicop* bicop; //!< a pointer to the bivariate copula to optimize.
        double par0;  //!< main dependence parameter.
        unsigned int objective_calls; //!< number of evaluations of the objective.
    } ParBicopOptData;

    //! @brief A class for the controls to NLopt
    class NLoptControls
    {
    public:

        NLoptControls();
        NLoptControls(double xtol_rel, double xtol_abs, double ftol_rel, 
            double ftol_abs, int maxeval);

        void set_controls(nlopt::opt* opt);

        double get_xtol_rel();
        double get_xtol_abs();
        double get_ftol_rel();
        double get_ftol_abs();
        double get_maxeval();

    private:
        double xtol_rel_; //! Relative tolerance on the parameter value
        double xtol_abs_; //! Absolute tolerance on the parameter value
        double ftol_rel_; //! Relative tolerance on the function value
        double ftol_abs_; //! Absolute tolerance on the function value
        int maxeval_; //! Maximal number of evaluations of the objective

        //! Sanity checks
        //! @{
        void check_parameters(double xtol_rel, double xtol_abs, double ftol_rel,
             double ftol_abs, int maxeval);
        //! @}
    };

    //! @brief A class for optimization (wrapping NLopt).
    class Optimizer {
    public:
        Optimizer(unsigned int n_parameters);
        Optimizer(unsigned int n_parameters, double xtol_rel, double xtol_abs,
                  double ftol_rel, double ftol_abs, int maxeval);
        
        void set_bounds(const Eigen::MatrixXd& lower_bounds,
            const Eigen::MatrixXd& upper_bounds);
        void set_objective(nlopt::vfunc f, void* f_data);

        Eigen::VectorXd optimize(Eigen::VectorXd initial_parameters);

    private:
        unsigned int n_parameters_;
        nlopt::opt opt_;
        NLoptControls controls_;
    };
}

}
