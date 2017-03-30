// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#pragma once

#include "abstract.hpp"
#include "fit_controls.hpp"

namespace vinecopulib {
    //! @brief A class for bivariate copula models.
    //! 
    //! The copula model is fully characterized by the family, rotation,
    //! and parameters. 
    class Bicop
    {
    public:
        // Constructors
        Bicop();
        Bicop(BicopFamily family, int rotation = 0,
              const Eigen::MatrixXd& parameters = Eigen::MatrixXd());
        Bicop(Eigen::Matrix<double, Eigen::Dynamic, 2> data,
              FitControlsBicop controls = FitControlsBicop());

        // Getters and setters
        BicopFamily get_family() const;
        std::string get_family_name() const;
        int get_rotation() const;
        Eigen::MatrixXd get_parameters() const;
        void set_rotation(int rotation);
        void set_parameters(const Eigen::MatrixXd& parameters);

        // Stats methods
        Eigen::VectorXd pdf(const Eigen::Matrix<double,Eigen::Dynamic,2>& u);
        Eigen::VectorXd hfunc1(const Eigen::Matrix<double,Eigen::Dynamic,2>& u);
        Eigen::VectorXd hfunc2(const Eigen::Matrix<double,Eigen::Dynamic,2>& u);
        Eigen::VectorXd hinv1(const Eigen::Matrix<double,Eigen::Dynamic,2>& u);
        Eigen::VectorXd hinv2(const Eigen::Matrix<double,Eigen::Dynamic,2>& u);
        Eigen::Matrix<double,Eigen::Dynamic,2> simulate(const int& n);


        // Methods modifying the family/rotation/parameters
        void fit(const Eigen::Matrix<double, Eigen::Dynamic, 2> &data,
                 FitControlsBicop controls = FitControlsBicop());
        void select(Eigen::Matrix<double, Eigen::Dynamic, 2> data,
                    FitControlsBicop controls = FitControlsBicop());

        // Fit statistics
        double loglik(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u);
        double aic(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u);
        double bic(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u);

        // Misc
        std::string str();
        double calculate_npars();
        double parameters_to_tau(const Eigen::VectorXd& parameters);
        Eigen::MatrixXd tau_to_parameters(const double& tau);
        void flip();

    private:
        Eigen::MatrixXd get_parameters_lower_bounds() const;
        Eigen::MatrixXd get_parameters_upper_bounds() const;
        Eigen::Matrix<double, Eigen::Dynamic, 2> cut_and_rotate(
                const Eigen::Matrix<double, Eigen::Dynamic, 2>& u);
        void check_rotation(int rotation);

        BicopPtr get_bicop();
        BicopPtr bicop_;
        int rotation_;
    };
}
