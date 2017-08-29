#' High Performance Algorithms for Vine Copula Modeling
#' 
#' 'vinecopulib' is a high performance C++ library based on
#' 'Boost', 'Eigen' and 'NLopt'. It provides high-performance implementations of 
#' the core features of the popular VineCopula package, in particular
#' inference algorithms for both vine copula and bivariate copula models.
#' Advantages over VineCopula are a sleaker and more modern API, shorter runtimes, 
#' especially in high dimensions, nonparametric and multi-parameter families.
#'
#' @name rvinecopulib
#' @docType package
#' @useDynLib rvinecopulib, .registration = TRUE
#' @importFrom Rcpp evalCpp
#'
#' @author Thomas Nagler, Thibault Vatter
#
#' @keywords package
#'
#' @examples
#' ## bicop_dist objects
#' bicop_dist("gaussian", 0, 0.5)
#' str(bicop_dist("gauss", 0, 0.5))
#' bicop <- bicop_dist("clayton", 90, 3)
#' 
#' ## bicop objects
#' u <- rbicop(500, "gauss", 0, 0.5)
#' fit1 <- bicop(u, "par")
#' fit1
#' 
#' ## vinecop_dist objects
#' ## specify pair-copulas
#' bicop <- bicop_dist("bb1", 90, c(3, 2))
#' pcs <- list(
#'     list(bicop, bicop),  # pair-copulas in first tree 
#'     list(bicop)          # pair-copulas in second tree 
#' )
#' ## specify R-vine matrix
#' mat <- matrix(c(1, 2, 3, 1, 2, 0, 1, 0, 0), 3, 3) 
#' ## build the vinecop_dist object
#' vc <- vinecop_dist(pcs, mat)
#' summary(vc)
#' 
#' ## vinecop objects
#' u <- sapply(1:3, function(i) runif(50))
#' vc <- vinecop(u, "par")
#' summary(vc)
NULL
