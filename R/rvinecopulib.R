#' High Performance Algorithms for Vine Copula Modeling
#'
#' Provides an interface to 'vinecopulib', a C++ library for vine copula
#' modeling based on 'Boost' and 'Eigen'. The 'rvinecopulib' package implements
#' the core features of the popular 'VineCopula' package, in particular
#' inference algorithms for both vine copula and bivariate copula models.
#' Advantages over 'VineCopula' are a sleeker and more modern API, improved
#' performances, especially in high dimensions, nonparametric and
#' multi-parameter families. The 'rvinecopulib' package includes 'vinecopulib'
#' as header-only C++ library (currently version 0.5.2). Thus users do not need
#' to install 'vinecopulib' itself in order to use 'rvinecopulib'. Since their
#' initial releases, 'vinecopulib' is licensed under the MIT License, and
#' 'rvinecopulib' is licensed under the GNU GPL version 3.
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
#' fit1 <- bicop(u, family = "par")
#' fit1
#'
#' ## vinecop_dist objects
#' ## specify pair-copulas
#' bicop <- bicop_dist("bb1", 90, c(3, 2))
#' pcs <- list(
#'   list(bicop, bicop), # pair-copulas in first tree
#'   list(bicop) # pair-copulas in second tree
#' )
#' ## specify R-vine matrix
#' mat <- matrix(c(1, 2, 3, 1, 2, 0, 1, 0, 0), 3, 3)
#' ## build the vinecop_dist object
#' vc <- vinecop_dist(pcs, mat)
#' summary(vc)
#'
#' ## vinecop objects
#' u <- sapply(1:3, function(i) runif(50))
#' vc <- vinecop(u, family = "par")
#' summary(vc)
#'
#' ## vine_dist objects
#' vc <- vine_dist(list(distr = "norm"), pcs, mat)
#' summary(vc)
#'
#' ## vine objects
#' x <- sapply(1:3, function(i) rnorm(50))
#' vc <- vine(x, copula_controls = list(family_set = "par"))
#' summary(vc)
NULL
