// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/misc/tools_optimization.hpp>
#include <vinecopulib/misc/tools_stats.hpp>

#include <cmath>

namespace vinecopulib {
    
//! Utilities for numerical optimization (based on NLopt)
namespace tools_optimization {

    //! creates an Optimizer using the default controls, see NLoptControls.
    //!
    //! @param n_parameters Number of parameters to optimize
    Optimizer::Optimizer(unsigned int n_parameters)
    {
        if (n_parameters < 1) {
            throw std::runtime_error("n_parameters should be larger than 0.");
        }
        n_parameters_ = n_parameters;
        opt_ = nlopt::opt(nlopt::LN_BOBYQA, n_parameters);
        controls_ = NLoptControls();
        controls_.set_controls(&opt_);
    }

    //! creates an optimizer using custom controls.
    //!
    //! @param n_parameters number of parameters to optimize.
    //! @param xtol_rel relative tolerance on the parameter value.
    //! @param xtol_abs absolute tolerance on the parameter value.
    //! @param ftol_rel relative tolerance on the function value.
    //! @param ftol_abs absolue tolerance on the function value.
    //! @param maxeval maximal number of evaluations of the objective.
    Optimizer::Optimizer(unsigned int n_parameters, double xtol_rel, double xtol_abs,
                         double ftol_rel, double ftol_abs, int maxeval)
    {
        if (n_parameters < 1) {
            throw std::runtime_error("n_parameters should be larger than 0.");
        }
        n_parameters_ = n_parameters;
        opt_ = nlopt::opt(nlopt::LN_BOBYQA, n_parameters);
        controls_ = NLoptControls(xtol_rel, xtol_abs, ftol_rel, ftol_abs, maxeval);
        controls_.set_controls(&opt_);
    }

    
    //! sets the optimizer's bounds.
    //! @param lower_bounds
    //! @param upper_bounds
    void Optimizer::set_bounds(const Eigen::MatrixXd& lower_bounds,
        const Eigen::MatrixXd& upper_bounds)
    {
        std::vector<double> lb(n_parameters_);
        std::vector<double> ub(n_parameters_);
        double eps = 1e-6;
        Eigen::VectorXd::Map(&lb[0], n_parameters_) = lower_bounds.array() + eps;
        Eigen::VectorXd::Map(&ub[0], n_parameters_) = upper_bounds.array() - eps;
        opt_.set_lower_bounds(lb);
        opt_.set_upper_bounds(ub);
    }

    //! Create controls using the default contructor
    //! 
    //! The defaults are 
    //! ```
    //! xtol_rel_ = 1e-6;
    //! xtol_abs_ = 1e-6;
    //! ftol_rel_ = 1e-6;
    //! ftol_abs_ = 1e-6;
    //! maxeval_ = 1000;
    //! ```
    NLoptControls::NLoptControls()
    {
        xtol_rel_ = 1e-6;
        xtol_abs_ = 1e-6;
        ftol_rel_ = 1e-6;
        ftol_abs_ = 1e-6;
        maxeval_ = 1000;
    }

    //! Create controls by passing the arguments
    //!
    //! @param xtol_rel relative tolerance on the parameter value.
    //! @param xtol_abs absolute tolerance on the parameter value.
    //! @param ftol_rel relative tolerance on the function value.
    //! @param ftol_abs absolue tolerance on the function value.
    //! @param maxeval maximal number of evaluations of the objective.
    NLoptControls::NLoptControls(double xtol_rel, double xtol_abs, 
        double ftol_rel, double ftol_abs, int maxeval)
    {
        check_parameters(xtol_rel, xtol_abs, ftol_rel, ftol_abs, maxeval);
        xtol_rel_ = xtol_rel;
        xtol_abs_ = xtol_abs;
        ftol_rel_ = ftol_rel;
        ftol_abs_ = ftol_abs;
        maxeval_ = maxeval;
    }

    void NLoptControls::check_parameters(double xtol_rel, double xtol_abs, 
                                         double ftol_rel, double ftol_abs, 
                                         int maxeval)
    {
        if (xtol_rel <= 0 || xtol_rel > 1) {
            throw std::runtime_error("xtol_rel should be in (0,1]");
        }
        if (ftol_rel <= 0 || ftol_rel > 1) {
            throw std::runtime_error("ftol_rel should be in (0,1]");
        }
        if (xtol_abs <= 0) {
            throw std::runtime_error("xtol_abs should be larger than 0");
        }
        if (ftol_abs <= 0) {
            throw std::runtime_error("ftol_abs should be larger than 0");
        }
        if (maxeval <= 0) {
            throw std::runtime_error("maxeval should be larger than 0");
        }
    }

    //! @name Getters and setters
    //! @{
    
    //! @return the relative tolerance on the parameter value. 
    double NLoptControls::get_xtol_rel() {return xtol_rel_;};
    //! @return the absolute tolerance on the parameter value. 
    double NLoptControls::get_xtol_abs() {return xtol_abs_;};
    //! @return the relative tolerance on the function value.
    double NLoptControls::get_ftol_rel() {return ftol_rel_;};
    //! @return the absolute tolerance on the function value.
    double NLoptControls::get_ftol_abs() {return ftol_abs_;};
    //! @return the maximal number of evaluations of the objective.
    double NLoptControls::get_maxeval() {return maxeval_;};
    
    //! sets controls of an Optimizer
    //! @param opt Pointer to the optimizer for control setting
    void NLoptControls::set_controls(nlopt::opt* opt)
    {
        opt->set_xtol_rel(xtol_rel_);
        opt->set_xtol_abs(xtol_abs_);
        opt->set_ftol_rel(ftol_rel_);
        opt->set_ftol_abs(ftol_abs_);
        opt->set_maxeval(maxeval_);
    };
    //! @}

    //! @name (Pseudo-) maximum likelihood estimation
    //! @param x the parameters.
    //! @param f_data a pointer to a ParBicopOptData object.
    //! @{
    
    //! evaluates the objective function for maximum likelihood estimation.
    double mle_objective(const std::vector<double>& x,
                         std::vector<double>&,
                         void* f_data)
    {
        ParBicopOptData* newdata = (ParBicopOptData*) f_data;
        ++newdata->objective_calls;
        Eigen::Map<const Eigen::VectorXd> par(&x[0], x.size());
        newdata->bicop->set_parameters(par);
        return (-1)*newdata->bicop->pdf(newdata->U).array().log().sum();
    }

    //! evaluates the objective function for profile maximum likelihood 
    //! estimation.
    double pmle_objective(const std::vector<double>& x,
                          std::vector<double>&,
                          void* f_data)
    {
        ParBicopOptData* newdata = (ParBicopOptData*) f_data;
        ++newdata->objective_calls;
        Eigen::VectorXd par = Eigen::VectorXd::Ones(x.size()+1);
        par(0) = newdata->par0;
        for (unsigned int i = 0; i < x.size(); ++i) {
            par(i + 1) = x[i];
        }
        newdata->bicop->set_parameters(par);
        double nll = newdata->bicop->pdf(newdata->U).array().log().sum();
        nll *= -1;
        return nll;
    }
    
    //! @}

    //! Set the optimizer's objective and data
    //!
    //! @param f The optimizer's objective function (see nlopt's documentation)
    //! @param f_data The optimizer's data (see nlopt's documentation)
    void Optimizer::set_objective(nlopt::vfunc f, void* f_data)
    {
        opt_.set_min_objective(f, f_data);
    }

    //! solve the optimization problem.
    //!
    //! @param initial_parameters of starting values for the optimization
    //!     algorithm.
    //! @return the optimal parameters.
    Eigen::VectorXd Optimizer::optimize(Eigen::VectorXd initial_parameters)
    {
        if (initial_parameters.size() != n_parameters_) {
            throw std::runtime_error("The size of x should be n_parameters_.");
        }

        double nll;
        std::vector<double> x(n_parameters_);
        Eigen::VectorXd::Map(&x[0], n_parameters_) = initial_parameters;
        try {
            opt_.optimize(x, nll);
        } catch (nlopt::roundoff_limited err) {
            throw std::runtime_error(std::string("Halted because roundoff errors limited progress! ") + err.what());
        } catch (nlopt::forced_stop err) {
            throw std::runtime_error(std::string("Halted because of a forced termination! ") + err.what());
        } catch (std::invalid_argument err) {
            throw std::runtime_error(std::string("Invalid arguments. ") + err.what());
        } catch (std::bad_alloc err) {
            throw std::runtime_error(std::string("Ran out of memory. ") + err.what());
        } catch (std::runtime_error err) {
            throw std::runtime_error(std::string("Generic failure. ") + err.what());
        } catch (...) {
            // do nothing for other errors (results are fine)
        }

        Eigen::Map<const Eigen::VectorXd> optimized_parameters(&x[0], x.size());
        return optimized_parameters;
    }
}

}
