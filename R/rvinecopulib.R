#' High Performance Algorithms for Vine Copula Modeling
#'
#' rvinecopulib provides high-performance tools for constructing, fitting,
#' selecting, simulating, and visualizing bivariate and vine copula models. It
#' also fits complete multivariate distributions by combining marginal models
#' with a vine copula. Continuous, integer-valued discrete, mixed, and
#' zero-inflated variables are supported.
#'
#' The package is the R interface to the header-only 'vinecopulib' C++ library,
#' which is bundled with rvinecopulib. Users do not need to install the C++
#' library separately. 'vinecopulib' is licensed under the MIT License and
#' rvinecopulib under the GNU GPL version 3.
#'
#' @section Main capabilities:
#' rvinecopulib provides:
#' \itemize{
#'   \item parametric, nonparametric, rotated, and extreme-value pair copulas;
#'   \item automatic pair-copula family, vine structure, truncation, and
#'     threshold selection;
#'   \item complete multivariate distributions with nonparametric,
#'     parametric, or custom margins;
#'   \item continuous, integer-valued discrete, mixed, and zero-inflated data;
#'   \item observation weights and parallel fitting;
#'   \item conditional simulation and Rosenblatt transforms;
#'   \item likelihood scores, Hessians, and observation-specific parameters;
#'   \item custom marginal families and custom tree-selection criteria.
#' }
#'
#' @section Modeling interfaces:
#' The framework is exposed at three levels. [vine()] models observations on
#' their original scale by fitting the marginal distributions and dependence
#' model together. [vinecop()] models uniform pseudo-observations when only the
#' dependence model is required. [bicop()] provides the corresponding
#' two-variable workflow.
#'
#' The constructors [vine_dist()], [vinecop_dist()], and [bicop_dist()] create
#' models from explicitly specified components.
#'
#' @section Learning more:
#' Start with `vignette("getting-started")`. Dedicated articles cover
#' bivariate copulas, vine structures and selection, marginal models,
#' discrete data, conditional simulation, and likelihood derivatives. The
#' [package website](https://vinecopulib.github.io/rvinecopulib/) contains the
#' rendered articles and complete function reference.
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
#' ## complete distribution on the original data scale
#' x <- data.frame(a = rnorm(50), b = rexp(50), n = rpois(50, 2))
#' fit <- vine(
#'   x,
#'   var_types = c("c", "c", "d"),
#'   copula_controls = list(family_set = "onepar")
#' )
#' rvine(3, fit)
#'
#' ## dependence model on the copula scale
#' u <- pseudo_obs(as.matrix(USArrests))
#' copula_fit <- vinecop(u, family_set = "onepar")
#' summary(copula_fit)
#'
#' ## bivariate model
#' uv <- rbicop(100, "gaussian", 0, 0.5)
#' bicop(uv, family_set = "par")
"_PACKAGE"
