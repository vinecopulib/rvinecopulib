// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#include <exception>
#include <cmath>

#include <vinecopulib/bicop/abstract.hpp>
#include <vinecopulib/misc/tools_stl.hpp>

#include <vinecopulib/bicop/bb1.hpp>
#include <vinecopulib/bicop/bb6.hpp>
#include <vinecopulib/bicop/bb7.hpp>
#include <vinecopulib/bicop/bb8.hpp>
#include <vinecopulib/bicop/clayton.hpp>
#include <vinecopulib/bicop/frank.hpp>
#include <vinecopulib/bicop/gaussian.hpp>
#include <vinecopulib/bicop/gumbel.hpp>
#include <vinecopulib/bicop/indep.hpp>
#include <vinecopulib/bicop/joe.hpp>
#include <vinecopulib/bicop/student.hpp>
#include <vinecopulib/bicop/tll0.hpp>

namespace vinecopulib
{
    //! Create a bivariate copula using the default contructor
    //!
    //! @param family the copula family.
    //! @param parameters the copula parameters (optional, must be compatible 
    //!     with family).
    //! @return A pointer to an object that inherits from AbstractBicop.
    //! @{
    BicopPtr AbstractBicop::create(BicopFamily family,
                                   const Eigen::MatrixXd& parameters)
    {
        BicopPtr new_bicop;
        switch (family) {
            case BicopFamily::indep:
                new_bicop = BicopPtr(new IndepBicop());
                break;
            case BicopFamily::gaussian:
                new_bicop = BicopPtr(new GaussianBicop());
                break;
            case BicopFamily::student:
                new_bicop = BicopPtr(new StudentBicop());
                break;
            case BicopFamily::clayton:
                new_bicop = BicopPtr(new ClaytonBicop());
                break;
            case BicopFamily::gumbel:
                new_bicop = BicopPtr(new GumbelBicop());
                break;
            case BicopFamily::frank:
                new_bicop = BicopPtr(new FrankBicop());
                break;
            case BicopFamily::joe:
                new_bicop = BicopPtr(new JoeBicop());
                break;
            case BicopFamily::bb1:
                new_bicop = BicopPtr(new Bb1Bicop());
                break;
            case BicopFamily::bb6:
                new_bicop = BicopPtr(new Bb6Bicop());
                break;
            case BicopFamily::bb7:
                new_bicop = BicopPtr(new Bb7Bicop());
                break;
            case BicopFamily::bb8:
                new_bicop = BicopPtr(new Bb8Bicop());
                break;
            case BicopFamily::tll0:
                new_bicop =  BicopPtr(new Tll0Bicop());
                break;

            default:
                throw std::runtime_error(std::string("Family not implemented"));
        }
        
        if (parameters.size() > 0) {
            new_bicop->set_parameters(parameters);
        }
            
        return new_bicop;
    }
        
    //!@}

    Eigen::VectorXd no_tau_to_parameters(const double&)
    {
        throw std::runtime_error("Method not implemented for this family");
    }

    //! Getters and setters.
    //! @{
    BicopFamily AbstractBicop::get_family() const 
    {
        return family_;
    }
    
    std::string AbstractBicop::get_family_name() const 
    {
        return vinecopulib::get_family_name(family_);
    };
    
    Eigen::MatrixXd AbstractBicop::get_parameters() const 
    {
        return parameters_;
    }
    
    Eigen::MatrixXd AbstractBicop::get_parameters_lower_bounds() const 
    {
        return parameters_lower_bounds_;
    }
    
    Eigen::MatrixXd AbstractBicop::get_parameters_upper_bounds() const 
    {
        return parameters_upper_bounds_;
    }

    void AbstractBicop::set_parameters(const Eigen::MatrixXd& parameters)
    {
        check_parameters(parameters);
        parameters_ = parameters;
    }
    //! @}
    
    //! Adjust the copula to flipped columns
    void AbstractBicop::flip()
    {
        // Most parametric families can be flipped by changing the rotation. 
        // This is done in Bicop::flip() directly. All other families need to
        // override this method.
    }
    

    //! Numerical inversion of h-functions
    //!
    //! These are generic functions to invert the hfunctions numerically.
    //! They can be used in derived classes to define \c hinv1 and \c hinv2.
    //!
    //! @param u \f$m \times 2\f$ matrix of evaluation points.
    //! @return The numerical inverse of h-functions.
    //! @{
    Eigen::VectorXd AbstractBicop::hinv1_num(const Eigen::Matrix<double, Eigen::Dynamic, 2> &u)
    {
        Eigen::Matrix<double, Eigen::Dynamic, 2> u_new = u;
        auto h1 = [&](const Eigen::VectorXd &v) {
            u_new.col(1) = v;
            return hfunc1(u_new);
        };
        return tools_eigen::invert_f(u.col(1), h1);
    }

    Eigen::VectorXd AbstractBicop::hinv2_num(const Eigen::Matrix<double, Eigen::Dynamic, 2> &u)
    {
        Eigen::Matrix<double, Eigen::Dynamic, 2> u_new = u;
        auto h1 = [&](const Eigen::VectorXd &x) {
            u_new.col(0) = x;
            return hfunc2(u_new);
        };

        return tools_eigen::invert_f(u.col(0), h1);
    }
    //! @}

    
    //! Sanity checks
    //! @{
    void AbstractBicop::check_parameters(const Eigen::MatrixXd& parameters)
    {
        check_parameters_size(parameters);
        check_parameters_lower(parameters);
        check_parameters_upper(parameters);
    }
    
    
    void AbstractBicop::check_parameters_size(const Eigen::MatrixXd& parameters)
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
    
    
    void AbstractBicop::check_parameters_lower(const Eigen::MatrixXd& parameters)
    {
        if (parameters_lower_bounds_.size() > 0) {
            std::stringstream message;
            if ((parameters.array() < parameters_lower_bounds_.array()).any()) {
                message <<
                    "parameters exceed lower bound " << 
                    " for " << get_family_name() << " copula; \n" << 
                    "bound: \n" << parameters_lower_bounds_ << "\n" <<
                    "actual: " << parameters << "\n";
                throw std::runtime_error(message.str().c_str());  
            }
        }
    }
    
    void AbstractBicop::check_parameters_upper(const Eigen::MatrixXd& parameters)
    {
        if (parameters_upper_bounds_.size() > 0) {
            std::stringstream message;
            if ((parameters.array() > parameters_upper_bounds_.array()).any()) {
                message <<
                    "parameters exceed upper bound " << 
                    " for " << get_family_name() << " copula; \n" << 
                    "bound: \n" << parameters_upper_bounds_ << "\n" <<
                    "actual: " << parameters << "\n";
                throw std::runtime_error(message.str().c_str());  
            }
        }
    }

    //! @}
}
