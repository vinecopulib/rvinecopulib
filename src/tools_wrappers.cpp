#include <RcppEigen.h>
#include "bicop_wrappers.hpp"

// [[Rcpp::export]]
Eigen::MatrixXd pseudo_obs_cpp(Eigen::MatrixXd x, std::string ties_method) 
{
    return vinecopulib::tools_stats::to_pseudo_obs(x, ties_method);
}