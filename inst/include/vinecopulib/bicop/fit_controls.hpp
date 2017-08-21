// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <vector>
#include <vinecopulib/bicop/family.hpp>

namespace vinecopulib {
    //! @brief A class for controlling fit of bivariate copula models.
    //!
    class FitControlsBicop
    {
    public:
        // Constructor
        FitControlsBicop(std::vector<BicopFamily> family_set = bicop_families::all,
                         std::string parametric_method = "mle",
                         std::string nonparametric_method = "quadratic",
                         double nonparametric_mult = 1.0,
                         std::string selection_criterion = "bic",
                         bool preselect_families = true);
        FitControlsBicop(std::string parametric_method);
        FitControlsBicop(std::string nonparametric_method,
                         double nonparametric_mult);
        
        // Getters
        std::vector<BicopFamily> get_family_set() const;
        std::string get_parametric_method() const;
        std::string get_nonparametric_method() const;
        double get_nonparametric_mult() const;
        std::string get_selection_criterion() const;
        bool get_preselect_families() const;

        // Setters
        void set_family_set(std::vector<BicopFamily> family_set);
        void set_parametric_method(std::string parametric_method);
        void set_nonparametric_method(std::string nonparametric_method);
        void set_nonparametric_mult(double nonparametric_mult);
        void set_selection_criterion(std::string selection_criterion);
        void set_preselect_families(bool preselect_families);

    private:
        std::vector<BicopFamily> family_set_;
        std::string parametric_method_;
        std::string nonparametric_method_;
        double nonparametric_mult_;
        std::string selection_criterion_;
        bool preselect_families_;

        void check_parametric_method(std::string parametric_method);
        void check_nonparametric_method(std::string nonparametric_method);
        void check_nonparametric_mult(double nonparametric_mult);
        void check_selection_criterion(std::string selection_criterion);
    };
}
