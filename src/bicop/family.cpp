// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://tvatter.github.io/vinecopulib/.

#include "bicop/family.hpp"
#include <stdexcept>

namespace vinecopulib {
    //! converts a BicopFamily into a string with its name.
    //! @param family the family.
    std::string get_family_name(BicopFamily family)
    {
        std::string family_name;
        switch (family) {
            case BicopFamily::indep:
                family_name = "Independence";
                break;
            case BicopFamily::gaussian:
                family_name = "Gaussian";
                break;
            case BicopFamily::student:
                family_name = "Student";
                break;
            case BicopFamily::clayton:
                family_name = "Clayton";
                break;
            case BicopFamily::gumbel:
                family_name = "Gumbel";
                break;
            case BicopFamily::frank:
                family_name = "Frank";
                break;
            case BicopFamily::joe:
                family_name = "Joe";
                break;
            case BicopFamily::bb1:
                family_name = "BB1";
                break;
            case BicopFamily::bb6:
                family_name = "BB6";
                break;
            case BicopFamily::bb7:
                family_name = "BB7";
                break;
            case BicopFamily::bb8:
                family_name = "BB8";
                break;
            case BicopFamily::tll0:
                family_name = "TLL0";
                break;

            default:
                throw std::runtime_error(std::string("Family not implemented"));
        }
        return family_name;
    };
}
