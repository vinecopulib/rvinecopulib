// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

// interface specfifc #defines can be set here 
// (R package does: #define INTERFACED_FROM_R)

#ifndef INTERFACED_FROM_R
    #include <iostream>
#endif

namespace vinecopulib {

    namespace tools_interface {
        inline void print(std::string text)
        {
            #ifndef INTERFACED_FROM_R
                std::cout << text;
            #else
                Rcpp::Rcout << text;
            #endif
        };
    }
}
