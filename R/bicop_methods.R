#' Bivariate copula distributions
#' 
#' Density, distribution function, random generation and h-functions (with 
#' their inverses) for the bivariate copula distribution.
#' 
#' @name bicop_distributions
#' @aliases dbicop pbicop rbicop hbicop dbicop_dist pbicop_dist rbicop_dist hbicop_dist
#' 
#' @param u evaluation points, either a length 2 vector or a two-column matrix.
#' @param family the copula family, a string containing the family name (see 
#'   \code{\link{bicop}} for all possible families).
#' @param rotation the rotation of the copula, one of `0`, `90`, `180`, `270`.
#' @param parameters a vector or matrix of copula parameters.
#'   
#' @note The functions can optionally be used with a [bicop_dist]
#' object, e.g., `dbicop(c(0.1, 0.5), bicop_dist("indep"))`.
#' 
#' @details 
#' See [bicop] for the various implemented copula families. 
#' 
#' H-functions (`hbicop()`) are conditional distributions derived
#' from a copula. If \eqn{C(u, v) = P(U \le u, V \le v)} is a copula, then 
#' \deqn{h_1(v | u) = P(U \le u | V = v),} \deqn{h_2(u | v) = P(V \le v | U =
#' u).}
#' 
#' @return 
#' `dbicop()` gives the density, `pbicop()` gives the distribution function, 
#' `rbicop()` generates random deviates, and `hbicop()` gives the h-functions 
#' (and their inverses).
#' 
#' The length of the result is determined by `n` for `rbicop()`, and 
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
#'    matrix of independent \eqn{U(0, 1)} variables, this simulates data 
#'    from `vinecop`.
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
#' @param family a copula family (see [bicop_dist()]) or a [bicop_dist] object.
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
#' @name par_to_tau
#' @rdname par_to_tau
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

#' Predictions and fitted values for a bivariate copula model
#' 
#' Predictions of the density, distribution function,  
#' h-functions (with their inverses) for a bivariate copula model.
#'
#' @name bicop_predict_and_fitted
#' @aliases predict.bicop fitted.bicop
#' @param object a `bicop` object.
#' @param newdata points where the fit shall be evaluated.
#' @param what what to predict, one of `"pdf"`, `"cdf"`, `"hfunc1"`, `"hfunc2"`, 
#'    `"hinv1"`, `"hinv2"`.
#' @param ... unused.
#' @return 
#' `fitted()` and `logLik()` have return values similar to [dbicop()], 
#' [pbicop()], and [hbicop()].
#' @examples
#' # Simulate and fit a bivariate copula model
#' u <- rbicop(500, "gauss", 0, 0.5)
#' fit <- bicop(u, "par")
#' 
#' # Predictions
#' all.equal(predict(fit, u, "hfunc1"), fitted(fit, "hfunc1"))
#' @rdname predict_bicop
#' @export
predict.bicop <- function(object, newdata, what = "pdf", ...) {
    stopifnot(what %in% c("pdf", "cdf", "hfunc1", "hfunc2", "hinv1", "hinv2"))
    newdata <- if_vec_to_matrix(newdata)
    switch(
        what,
        "pdf"    = bicop_pdf_cpp(newdata, object),
        "cdf"    = bicop_cdf_cpp(newdata, object),
        "hfunc1" = bicop_hfunc1_cpp(newdata, object),
        "hfunc2" = bicop_hfunc2_cpp(newdata, object),
        "hinv1"  = bicop_hinv1_cpp(newdata, object),
        "hinv2"  = bicop_hinv2_cpp(newdata, object)
    )
}

#' @rdname predict_bicop
#' @export
fitted.bicop <- function(object, what = "pdf", ...) {
    if (is.null(object$data))
        stop("data have not been stored, use keep_data = TRUE when fitting.")
    stopifnot(what %in% c("pdf", "cdf", "hfunc1", "hfunc2", "hinv1", "hinv2"))
    switch(
        what,
        "pdf"    = bicop_pdf_cpp(object$data, object),
        "cdf"    = bicop_cdf_cpp(object$data, object),
        "hfunc1" = bicop_hfunc1_cpp(object$data, object),
        "hfunc2" = bicop_hfunc2_cpp(object$data, object),
        "hinv1"  = bicop_hinv1_cpp(object$data, object),
        "hinv2"  = bicop_hinv2_cpp(object$data, object)
    )
}

#' @importFrom stats logLik
#' @export
logLik.bicop <- function(object, ...) {
    if (is.null(object$data))
        stop("data have not been stored, use keep_data = TRUE when fitting.")
    structure(bicop_loglik_cpp(object$data, object), "df" = object$npars)
}

#' @export
print.bicop_dist <- function(x, ...) {
    if (x$family %in% setdiff(family_set_nonparametric, "indep")) {
        x$parameters <- paste0(round(x$npars, 2), sep = " d.f.")
    }
    cat("Bivariate copula ('bicop_dist'): ",
        "family = ", x$family,
        ", rotation = ", x$rotation,
        ", parameters = ", ifelse(length(x$parameters) > 1, 
                                  paste(round(x$parameters, 2), 
                                        collapse = ", "),
                                  x$parameters),
        sep = "")
    cat("\n")
}

#' @export
summary.bicop_dist <- function(object, ...) {
    print.bicop_dist(object, ...)
}

#' @export
print.bicop <- function(x, ...) {
    if (x$family %in% setdiff(family_set_nonparametric, "indep")) {
        pars_formatted <- paste0(round(x$npars, 2), sep = " d.f.")
    } else {
        pars_formatted <- paste(round(x$parameters, 2), collapse = ", ")
    }
    cat("Bivariate copula fit ('bicop'): ",
        "family = ", x$family,
        ", rotation = ", x$rotation,
        ", parameters = ", pars_formatted,
        "\n",
        sep = "")
    
    invisible(x)
}

#' @export
summary.bicop <- function(object, ...) {
    print.bicop(object, ...)
    cat("nobs =", object$nobs, "  ")
    if (!is.null(object$dat))
        info <- bicop_fit_info(object)
    if (!is.null(object$dat)) {
        cat("logLik =", round(info$logLik, 2), "  ")
        cat("npars =", round(info$npars, 2), "  ")
        cat("AIC =", round(info$AIC, 2), "  ")
        cat("BIC =", round(info$BIC, 2), "  ")
        attr(object, "info") <- info
    } else {
        cat("(for mor information, fit model with keep_data = TRUE)")
    }
    cat("\n")
    invisible(object)
}

bicop_fit_info <- function(bc) {
    ll <- logLik(bc)
    list(
        nobs   = bc$nobs,
        logLik = ll[1],
        npars  = attr(ll, "df"),
        AIC    = -2 * ll[1] + 2 * attr(ll, "df"),
        BIC    = -2 * ll[1] + log(bc$nobs) * attr(ll, "df")
    )
}