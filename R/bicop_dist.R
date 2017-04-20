#' Bivariate copula distributions
#'
#' A bivariate copula distribution is specified by:
#'
#' @param family the copula family, a string containing the family name (see
#' *Details* for all possible families).
#' @param rotation the rotation of the copula, one of `0`, `90`, `180` `270`.
#' @param parameters a vector or matrix of copula paramters.
#' 
#' @note 
#' The evaluation functions can optionally be used with a `bicop_dist` object, 
#' e.g., `dbicop(c(0.1, 0.5), bicop_dist("indep"))`.
#'
#' @details
#' The implemented families listed below. Partial matching is activated, i.e., 
#' `"gauss"` is equivalent to `"gaussian"`.
#' \describe{
#' \item{`indep`}{Independence copula.}
#' \item{`gaussian`}{Gaussian copula.}
#' \item{`t`, `student`}{Student t copula.}
#' \item{`clayton`}{Clayton copula.}
#' \item{`gumbel`}{Gumbel copula.}
#' \item{`frank`}{Frank copula.}
#' \item{`joe`}{Joe copula.}
#' \item{`bb1`}{BB1 copula.}
#' \item{`bb6`}{BB6 copula.}
#' \item{`bb7`}{BB7 copula.}
#' \item{`bb8`}{BB8 copula.}
#' }
#' 
#'
#' @return An object of class `bicop_dist`.
#'
#' @examples
#' bicop_dist("gaussian", 0, 3)
#' str(bicop_dist("gauss", 0, 3))
#' 
#' bicop <- bicop_dist("clayton", 90, 3)
#' @export
bicop_dist <- function(family = "indep", rotation = 0, parameters = numeric(0)) {
    stopifnot(length(family) == 1)
    if (family == "t")
        family <- "student"
    family <- family_set_all[pmatch(family, family_set_all)]
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

#' @rdname bicop_dist
#' @examples
#' # evaluate the copula density
#' dbicop(c(0.1, 0.2), "clay", 90, 3)
#' dbicop(c(0.1, 0.2), bicop)
#' 
#' @export
dbicop <- function(u, family, rotation, parameters) {
    bicop <- args2bicop(as.list(environment()))
    bicop_pdf_cpp(if_vec_to_matrix(u), bicop)
}

#' @rdname bicop_dist
#' @param n number of observations. If `length(n) > 1``, the length is taken to 
#'   be the number required.
#' @examples
#' # evaluate the copula density
#' plot(rbicop(500, "clay", 90, 3))
#' plot(rbicop(500, bicop))
#' 
#' @export
rbicop <- function(n, family, rotation, parameters) {
    if (length(n) > 1)
        n <- length(n)
    bicop <- args2bicop(as.list(environment()))
    bicop_simulate_cpp(n, bicop)
}

#' @rdname bicop_dist
#' @examples
#' # evaluate the copula density
#' h1bicop(c(0.1, 0.2), "frank", 0, 5)
#' h2bicop(c(0.1, 0.2), bicop)
#' 
#' @export
h1bicop <- function(u, family, rotation, parameters) {
    bicop <- args2bicop(as.list(environment()))
    bicop_hfunc1_cpp(if_vec_to_matrix(u), bicop)
}

#' @rdname bicop_dist
#' @export
h2bicop <- function(u, family, rotation, parameters) {
    bicop <- args2bicop(as.list(environment()))
    bicop_hfunc2_cpp(if_vec_to_matrix(u), bicop)
}

#' @rdname bicop_dist
#' @examples
#' # evaluate the copula density
#' hi1bicop(c(0.1, 0.2), "bb6", 180, 5)
#' hi2bicop(c(0.1, 0.2), bicop)
#' 
#' @export
hi1bicop <- function(u, family, rotation, parameters) {
    bicop <- args2bicop(as.list(environment()))
    bicop_hinv1_cpp(if_vec_to_matrix(u), bicop)
}

#' @rdname bicop_dist
#' @export
hi2bicop <- function(u, family, rotation, parameters) {
    bicop <- args2bicop(as.list(environment()))
    bicop_hinv2_cpp(if_vec_to_matrix(u), bicop)
}

