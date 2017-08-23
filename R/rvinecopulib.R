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
#' # bivariate copula
#' dist <- bicop_dist("clayton", 90, 3)
#' u <- rbicop(100, dist)
#' fit <- bicop(u)
#' fit
#' summary(fit)
#' 
#' # vine copula
#' u <- sapply(1:3, function(i) runif(50))
#' fit <- vinecop(u, "par")
#' fit
#' summary(fit)
NULL
