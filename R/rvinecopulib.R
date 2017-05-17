#' High Performance Algorithms for Vine Copula Modeling
#'
#' @name rvinecopulib
#' @docType package
#' @useDynLib rvinecopulib
#' @importFrom Rcpp evalCpp
#'
#' @author Thomas Nagler, Thibault Vatter
#
#' @keywords package
#'
#' @examples
#' dist <- bicop_dist("clayton", 90, 3)
#' u <- rbicop(100, dist)
#' fit <- bicop_fit(u)
#' fit
NULL
