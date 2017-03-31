#' Bivariate copula distributions
#'
#' A bivariate copula distribution is specified by:
#'
#' @param family the copula family, a string containing the family name (see
#' *Details* for all possible families).
#' @param rotation the rotation of the copula, one of `0`, `90`, `180` `270`.
#' @param parameters a vector or matrix of copula paramters.
#'
#' @details
#' The implemented families are:
#' \describe{
#' \item{`indep`}{Independence copula.}
#' \item{`gaussian`}{Gaussian copula.}
#' \item{`t`, `student`}{Student t copula.}
#' \item{`clayton`}{Clayton copula.}
#' \item{`gumbel`}{Gumbel copula.}
#' \item{`frank`}{Frank copula.}
#' \item{`Joe`}{Joe copula.}
#' \item{`bb1`}{BB1 copula.}
#' \item{`bb6`}{BB6 copula.}
#' \item{`bb7`}{BB7 copula.}
#' \item{`bb8`}{BB8 copula.}
#' }
#'
#' @return An object of class `bicop_dist`.
#'
#' @examples
#' bicop_dist("clayton", 0, 3)
#'
#' str(bicop_dist("clayton", 0, 3))
#'
#' @export
bicop_dist <- function(family = "indep", rotation = 0, parameters = numeric(0)) {
    stopifnot(length(family) == 1)
    if (family == "t")
        family <- "student"
    dist <- list(family = family,
                 rotation = rotation,
                 parameters = as.matrix(parameters))
    bicop_check_cpp(dist)
    structure(dist, class = "bicop_dist")
}

#' @export
print.bicop_dist <- function(x, ...) {
    cat("Bivariate copula ('bicop_dist'): ",
        "family = ", x$family,
        ", rotation = ", x$rotation,
        ", parameters = ", x$parameters,
        sep = "")
}
