#' @param family the copula family, a string containing the family name (see
#' *Details* for all possible families).
#' @param rotation the rotation of the copula, one of `0`, `90`, `180`, `270`.
#' @param parameters a vector or matrix of copula paramters.
#' @rdname bicop
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

#' Bivariate copula distributions
#' 
#' Density, distribution function, random generation and h-functions (with 
#' their inverses) for the bivariate copula distribution.
#' 
#' @aliases dbicop pbicop rbicop hbicop dbicop_dist pbicop_dist rbicop_dist hbicop_dist
#' 
#' @param u evaluation points, either a length 2 vector or a two-column matrix.
#' @param family the copula family, a string containing the family name (see 
#'   \code{\link{bicop}} for all possible families).
#' @param rotation the rotation of the copula, one of `0`, `90`, `180`, `270`.
#' @param parameters a vector or matrix of copula paramters.
#'   
#' @note The functions can optionally be used with a `bicop_dist`
#' object, e.g., `dbicop(c(0.1, 0.5), bicop_dist("indep"))`.
#' 
#' @details 
#' See \code{\link{bicop}} for the various implemented copula families. 
#' H-functions (`hbicop()`) are conditional distributions derived
#' from a copula. If \eqn{C(u, v) = P(U \le u, V \le v)} is a copula, then 
#' \deqn{h_1(v | u) = P(U \le u | V = v),} \deqn{h_2(u | v) = P(V \le v | U =
#' u).}
#' 
#' @return 
#' `dbicop` gives the density, `pbicop` gives the distribution function, 
#' `rbicop` generates random deviates, and `hbicop` gives the h-functions 
#' (and their inverses).
#' 
#' The length of the result is determined by `n` for `rbicop`, and 
#' the number of rows in `u` for the other functions.
#' 
#' The numerical arguments other than `n` are recycled to the length of the 
#' result.
# 
#' @examples
#' ## evaluate the copula density
#' dbicop(c(0.1, 0.2), "clay", 90, 3)
#' dbicop(c(0.1, 0.2), bicop_dist("clay", 90, 3))
#' 
#' ## evaluate the copula cdf
#' pbicop(c(0.1, 0.2), "clay", 90, 3)
#' 
#' ## simulate data
#' plot(rbicop(500, "clay", 90, 3))
#' 
#' ## h-functions
#' joe_cop <- bicop_dist("joe", 0, 3)
#' # h_1(0.2 | 0.1)
#' hbicop(c(0.1, 0.2), 1, "bb8", 0, c(2, 0.5))
#' # h_2(0.1 | 0.2)
#' hbicop(c(0.1, 0.2), 2, joe_cop)
#' # h_1^{-1}(0.2 | 0.1)
#' hbicop(c(0.1, 0.2), 1, "bb8", 0, c(2, 0.5), inverse = TRUE)
#' # h_2^{-1}(0.1 | 0.2)
#' hbicop(c(0.1, 0.2), 2, joe_cop, inverse = TRUE)
#' @rdname bicop_methods
#' @export
dbicop <- function(u, family, rotation, parameters) {
    bicop <- args2bicop(family, rotation, parameters)
    bicop_pdf_cpp(if_vec_to_matrix(u), bicop)
}
#' @rdname bicop_methods
#' @export
pbicop <- function(u, family, rotation, parameters) {
    bicop <- args2bicop(family, rotation, parameters)
    bicop_cdf_cpp(if_vec_to_matrix(u), bicop)
}

#' @param n number of observations. If `length(n) > 1``, the length is taken to
#'   be the number required.
#' @param U optionally, an \eqn{n \times 2} matrix of values in \eqn{(0,1)}.
#'    The result is then the inverse Rosenblatt transform of `U`; if `U` is a
#'    matrix of independent \eqn{U(0, 1)} variebls, this simulates data 
#'    from `vinecop`.
#' @examples
#' 
#' ## simulate data
#' plot(rbicop(500, "clay", 90, 3))
#' plot(rbicop(500, bicop))
#' @rdname bicop_methods
#' @export
rbicop <- function(n, family, rotation, parameters, U = NULL) {
    if (length(n) > 1)
        n <- length(n)
    bicop <- args2bicop(family, rotation, parameters)
    U <- prep_uniform_data(n, 2, U)
    U <- cbind(U[, 1], bicop_hinv1_cpp(U, bicop))
    if (!is.null(bicop$names)) 
        colnames(U) <- bicop$names
    
    U
}


#' @rdname bicop_methods
#' @param cond_var either `1` or `2`; `cond_var = 1` conditions on the first
#'    variable, `cond_var = 2` on the second.
#' @param inverse whether to compute the h-function or its inverse.
#' @export
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
#' @param tau Kendall's \eqn{\tau}.
#' @export
tau_to_par <- function(family, tau) {
    bicop <- args2bicop(family)
    if (!(bicop$family %in% c(family_set_elliptical, family_set_nonparametric)))
        bicop$rotation <- ifelse(tau > 0, 0, 90)
    bicop_tau_to_par_cpp(bicop, tau)
}
