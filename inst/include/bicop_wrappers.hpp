#pragma once

#include "vinecopulib.hpp"

vinecopulib::BicopFamily to_cpp_family(const std::string& fam);
std::string to_r_family(const vinecopulib::BicopFamily& fam);
vinecopulib::Bicop bicop_wrap(const Rcpp::List& bicop_r);
Rcpp::List bicop_wrap(vinecopulib::Bicop bicop_cpp, bool is_fitted = FALSE);
