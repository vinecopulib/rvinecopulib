#pragma once

#include "vinecopulib.hpp"

vinecopulib::Bicop bicop_wrap(const Rcpp::List& bicop_r);
Rcpp::List bicop_wrap(vinecopulib::Bicop bicop_cpp);
    