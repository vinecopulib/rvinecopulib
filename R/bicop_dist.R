#' Bivariate copula distributions
#'
#' A bivariate copula distribution is specified by:
#'
#' @param family the copula family, a string containing the family name (see
#' *Details* for all possible families).
#' @param rotation the rotation of the copula, one of `0`, `90`, `180`, `270`.
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
#' \item{`t`}{Student t copula.}
#' \item{`clayton`}{Clayton copula.}
#' \item{`gumbel`}{Gumbel copula.}
#' \item{`frank`}{Frank copula.}
#' \item{`joe`}{Joe copula.}
#' \item{`bb1`}{BB1 copula.}
#' \item{`bb6`}{BB6 copula.}
#' \item{`bb7`}{BB7 copula.}
#' \item{`bb8`}{BB8 copula.}
#' \item{`tll0`}{local constant transformation kernel, should only be used with
#' data, see `bicop()`.}
#' }
#'
#'
#' @return An object of class `bicop_dist`.
#'
#' @examples
#' bicop_dist("gaussian", 0, 0.5)
#' str(bicop_dist("gauss", 0, 0.5))
#'
#' bicop <- bicop_dist("clayton", 90, 3)
#' @export
bicop_dist <- function(family = "indep", rotation = 0, parameters = numeric(0)) {
    stopifnot(length(family) == 1)
    if (family %in% setdiff(family_set_nonparametric, "indep"))
        stop("bicop_dist should not be used directly with nonparametric families.")
    family <- family_set_all[pmatch(family, family_set_all)]
    dist <- list(family     = family,
                 rotation   = rotation,
                 parameters = as.matrix(parameters),
                 npars      = length(parameters))
    bicop_check_cpp(dist)
    structure(dist, class = "bicop_dist")
}

#' @export
print.bicop_dist <- function(x, ...) {
    if (x$family %in% setdiff(family_set_nonparametric, "indep")) {
        x$parameters <- "[30x30 grid]"
    }
    cat("Bivariate copula ('bicop_dist'): ",
        "family = ", x$family,
        ", rotation = ", x$rotation,
        ", parameters = ", x$parameters,
        sep = "")
}

#' @rdname bicop_dist
#' @param u evaluation points, either a length 2 vector or a two-column matrix.
#' @examples
#' # evaluate the copula density
#' dbicop(c(0.1, 0.2), "clay", 90, 3)
#' dbicop(c(0.1, 0.2), bicop)
#' @export
dbicop <- function(u, family, rotation, parameters) {
    bicop <- args2bicop(family, rotation, parameters)
    bicop_pdf_cpp(if_vec_to_matrix(u), bicop)
}
#' @rdname bicop_dist
#' @examples
#' # evaluate the copula cdf
#' pbicop(c(0.1, 0.2), "clay", 90, 3)
#' @export
pbicop <- function(u, family, rotation, parameters) {
    bicop <- args2bicop(family, rotation, parameters)
    bicop_cdf_cpp(if_vec_to_matrix(u), bicop)
}

#' @rdname bicop_dist
#' @param n number of observations. If `length(n) > 1``, the length is taken to
#'   be the number required.
#' @examples
#' # simulate data
#' plot(rbicop(500, "clay", 90, 3))
#' plot(rbicop(500, bicop))
#'
#' @export
rbicop <- function(n, family, rotation, parameters) {
    if (length(n) > 1)
        n <- length(n)
    bicop <- args2bicop(family, rotation, parameters)
    bicop_simulate_cpp(n, bicop)
}


#' H-functions and their inverses for bivariate copula distributions
#'
#' @param u evaluation points, either a length 2 vector or a two-column matrix.
#' @param cond_var either `1` or `2`; `cond_var = 1` conditions on the first
#'    variable, `cond_var = 2` on the second.
#' @param family the copula family, a string containing the family name (see
#' [`bicop_dist()`]).
#' @param rotation the rotation of the copula, one of `0`, `90`, `180`, `270`.
#' @param parameters a vector or matrix of copula paramters.
#' @param inverse whether to compute the h-function or its inverse.
#'
#' @details H-functions are conditional distributions derived from a copula.
#' If \eqn{C(u, v) = P(U \le u, V \le v)} is a copula, then
#' \deqn{h_1(v | u) = P(U \le u | V = v),}
#' \deqn{h_2(u | v) = P(V \le v | U = u).}
#'
#' @return A numeric vector containing the value of the (inverse) h-function.
#' @export
#'
#' @examples
#' joe_cop <- bicop_dist("joe", 0, 3)
#'
#' # h_1(0.2 | 0.1)
#' hbicop(c(0.1, 0.2), 1, "bb8", 0, c(2, 0.5))
#'
#' # h_2(0.1 | 0.2)
#' hbicop(c(0.1, 0.2), 2, joe_cop)
#'
#' # h_1^{-1}(0.2 | 0.1)
#' hbicop(c(0.1, 0.2), 1, "bb8", 0, c(2, 0.5), inverse = TRUE)
#'
#' # h_2^{-1}(0.1 | 0.2)
#' hbicop(c(0.1, 0.2), 2, joe_cop, inverse = TRUE)
hbicop <- function(u, cond_var, family, rotation, parameters, inverse = FALSE) {
    stopifnot(length(cond_var) == 1)
    stopifnot(cond_var %in% c(1, 2))
    stopifnot(is.logical(inverse))
    bicop <- args2bicop(family, rotation, parameters)

    if (!inverse) {
        if (cond_var == 1) {
            return(bicop_hfunc1_cpp(if_vec_to_matrix(u), bicop))
        } else {
            return(bicop_hfunc2_cpp(if_vec_to_matrix(u), bicop))
        }
    } else {
        if (cond_var == 1) {
            return(bicop_hinv1_cpp(if_vec_to_matrix(u), bicop))
        } else {
            return(bicop_hinv2_cpp(if_vec_to_matrix(u), bicop))
        }
    }
}

#' Conversion between Kendall's tau and parameters
#'
#' @param family a copula family (see `bicop_dist()`) or a `bicop_dist` object.
#' @param rotation the rotation of the copula, one of `0`, `90`, `180`, `270`.
#' @param parameters vector or matrix of copula parameters, not used when 
#'   `family` is a `bicop_dist` object.
#'
#' @export
#'
#' @examples
#' # the following are equivalent
#' par_to_tau(bicop_dist("clayton", 0, 3))
#' par_to_tau("clayton", 0, 3)
#' 
#' tau_to_par("clayton", 0.5)
#' tau_to_par(bicop_dist("clayton", 0, 3), 0.5)
par_to_tau <- function(family, rotation, parameters) {
    bicop <- args2bicop(family, rotation, parameters)
    bicop_par_to_tau_cpp(bicop)
}

#' @rdname par_to_tau
#' @export
tau_to_par <- function(family, tau) {
    bicop <- args2bicop(family)
    if (!(bicop$family %in% c(family_set_elliptical, family_set_nonparametric)))
        bicop$rotation <- ifelse(tau > 0, 0, 90)
    bicop_tau_to_par_cpp(bicop, tau)
}
