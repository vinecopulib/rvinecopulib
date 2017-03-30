// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include <memory>
#include <vector>

#include "misc/tools_eigen.hpp"
#include "bicop/family.hpp"

// Pre-declaration to allow AbtractBicop to befriend the two functions
namespace tools_optimization
{
    // the objective function for maximum likelihood estimation
    double mle_objective(const std::vector<double>& x,
                         std::vector<double>& grad,
                         void* data);

    // the objective function for profile maximum likelihood estimation
    double pmle_objective(const std::vector<double>& x,
                          std::vector<double> &,
                          void* data);
}

namespace vinecopulib
{
    //! @brief An abstract class for bivariate copula families
    //!
    //! This class is used in the implementation underlying the Bicop class. 
    //! Users should not use AbstractBicop or derived classes directly, but 
    //! always work with the Bicop interface.
    class AbstractBicop
    {
    friend class Bicop;
    friend double tools_optimization::mle_objective(
        const std::vector<double>& x, std::vector<double>& grad, void* data);
    friend double tools_optimization::pmle_objective(
        const std::vector<double>& x, std::vector<double>& grad, void* data);
    
    protected:
        // Factories
        static std::shared_ptr<AbstractBicop> create(
            BicopFamily family = BicopFamily::indep,
            const Eigen::MatrixXd& parameters = Eigen::MatrixXd());

        // Getters and setters
        BicopFamily get_family() const;
        std::string get_family_name() const;
        Eigen::MatrixXd get_parameters() const;
        Eigen::MatrixXd get_parameters_lower_bounds() const;
        Eigen::MatrixXd get_parameters_upper_bounds() const;
        void set_parameters(const Eigen::MatrixXd& parameters);
        void flip();

        // Virtual methods
        virtual void fit(const Eigen::Matrix<double, Eigen::Dynamic, 2> &data,
                std::string method, double mult) = 0;
        virtual double calculate_npars() = 0;
        virtual double parameters_to_tau(const Eigen::VectorXd& parameters) = 0;
        virtual Eigen::VectorXd pdf(
                const Eigen::Matrix<double, Eigen::Dynamic, 2>& u) = 0;
        virtual Eigen::VectorXd hfunc1(
                const Eigen::Matrix<double, Eigen::Dynamic, 2>& u) = 0;
        virtual Eigen::VectorXd hfunc2(
                const Eigen::Matrix<double, Eigen::Dynamic, 2>& u) = 0;
        virtual Eigen::VectorXd hinv1(
                const Eigen::Matrix<double, Eigen::Dynamic, 2>& u) = 0;
        virtual Eigen::VectorXd hinv2(
                const Eigen::Matrix<double, Eigen::Dynamic, 2>& u) = 0;
        virtual Eigen::MatrixXd tau_to_parameters(const double& tau) = 0;

        // Misc methods
        Eigen::VectorXd hinv1_num(
                const Eigen::Matrix<double, Eigen::Dynamic, 2>& u);
        Eigen::VectorXd hinv2_num(
                const Eigen::Matrix<double, Eigen::Dynamic, 2>& u);

        // Data members
        BicopFamily family_;
        Eigen::MatrixXd parameters_;
        Eigen::MatrixXd parameters_lower_bounds_;
        Eigen::MatrixXd parameters_upper_bounds_;

    private:
        void check_parameters(const Eigen::MatrixXd& parameters);
        void check_parameters_size(const Eigen::MatrixXd& parameters);
        void check_parameters_upper(const Eigen::MatrixXd& parameters);
        void check_parameters_lower(const Eigen::MatrixXd& parameters);
    };
    
    //! A shared pointer to an object of class AbstracBicop.
    typedef std::shared_ptr<AbstractBicop> BicopPtr;
    Eigen::VectorXd no_tau_to_parameters(const double&);
}
