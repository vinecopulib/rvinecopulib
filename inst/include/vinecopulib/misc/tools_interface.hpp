// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

// interface specfifc #defines can be set here 
// (R package does: #define INTERFACED_FROM_R)
#define INTERFACED_FROM_R


// interface specific headers
#ifdef INTERFACED_FROM_R
    #include <RcppEigen.h>
    // [[Rcpp::depends(RcppProgress)]]
    #include <progress.hpp>
#else
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
        }
    
        inline void check_user_interrupt(bool do_check = true)
        {
            if (do_check) {
                #ifdef INTERFACED_FROM_R
                    if (Progress::check_abort()) {
                        throw std::runtime_error("C++ call interrupted by user");
                    }
                #endif
            }
        }
    }
}
