#' Bivariate copula distributions
#'
#' Density, distribution function, random generation and h-functions (with their
#' inverses) for the bivariate copula distribution.
#'
#' @name bicop_distributions
#' @aliases dbicop pbicop rbicop hbicop dbicop_dist pbicop_dist rbicop_dist
#'   hbicop_dist
#'
#' @param u evaluation points, a matrix with at least two columns, see
#'   *Details*.
#' @param family the copula family, a string containing the family name (see
#'   \code{\link{bicop}} for all possible families).
#' @param rotation the rotation of the copula, one of `0`, `90`, `180`, `270`.
#' @param parameters a vector or matrix of copula parameters.
#' @param var_types variable types, a length 2 vector; e.g., `c("c", "c")` for
#'   both continuous (default), or `c("c", "d")` for first variable continuous
#'   and second discrete.
#'
#' @note The functions can optionally be used with a [bicop_dist] object, e.g.,
#'   `dbicop(c(0.1, 0.5), bicop_dist("indep"))`.
#'
#' @details
#' See [bicop] for the various implemented copula families.
#'
#' H-functions (`hbicop()`) are conditional distributions derived
#' from a copula. If \eqn{C(u, v) = P(U \le u, V \le v)} is a copula, then
#' \deqn{h_1(u, v) = P(V \le v | U = u) = \partial C(u, v) / \partial u,}
#' \deqn{h_2(u, v) = P(U \le u | V = v) = \partial C(u, v) / \partial v.}
#' In other words, the H-function number refers to the conditioning variable.
#' When inverting H-functions, the inverse is then taken with respect to the
#' other variable, that is `v` when `cond_var = 1` and `u` when `cond_var = 2`.
#'
#' ## Discrete variables
#' When at least one variable is discrete, mote than two columns are required
#' for `u`: the first \eqn{n \times 2} block contains realizations of
#' \eqn{F_{X_1}(x_1), F_{X_2}(x_2)}. The second \eqn{n \times 2} block contains
#' realizations of \eqn{F_{X_1}(x_1^-), F_{X_1}(x_1^-)}. The minus indicates a
#' left-sided limit of the cdf. For, e.g., an integer-valued variable, it holds
#' \eqn{F_{X_1}(x_1^-) = F_{X_1}(x_1 - 1)}. For continuous variables the left
#' limit and the cdf itself coincide. Respective columns can be omitted in the
#' second block.
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
#'
#' @seealso [bicop_dist()], [bicop()]
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
#' # h_1(0.1, 0.2)
#' hbicop(c(0.1, 0.2), 1, "bb8", 0, c(2, 0.5))
#' # h_2^{-1}(0.1, 0.2)
#' hbicop(c(0.1, 0.2), 2, joe_cop, inverse = TRUE)
#'
#' ## mixed discrete and continuous data
#' x <- cbind(rpois(10, 1), rnorm(10, 1))
#' u <- cbind(ppois(x[, 1], 1), pnorm(x[, 2]), ppois(x[, 1] - 1, 1))
#' pbicop(u, "clay", 90, 3, var_types = c("d", "c"))
#'
#' @rdname bicop_methods
#' @export
dbicop <- function(u, family, rotation, parameters, var_types = c("c", "c")) {
  bicop <- args2bicop(family, rotation, parameters, var_types)
  bicop_pdf_cpp(if_vec_to_matrix(u), bicop)
}
#' @rdname bicop_methods
#' @export
pbicop <- function(u, family, rotation, parameters, var_types = c("c", "c")) {
  bicop <- args2bicop(family, rotation, parameters, var_types)
  bicop_cdf_cpp(if_vec_to_matrix(u), bicop)
}

#' @param n number of observations. If `length(n) > 1``, the length is taken to
#'   be the number required.
#' @param qrng if `TRUE`, generates quasi-random numbers using the bivariate
#' Generalized Halton sequence (default `qrng = FALSE`).
#' @rdname bicop_methods
#' @export
rbicop <- function(n, family, rotation, parameters, qrng = FALSE) {
  if (length(n) > 1) {
    n <- length(n)
  }
  if (inherits(family, "bicop_dist") & !missing(rotation)) {
    qrng <- rotation
  }
  assert_that(is.flag(qrng))

  bicop <- args2bicop(family, rotation, parameters)
  U <- bicop_sim_cpp(bicop, n, qrng, get_seeds())
  if (!is.null(bicop$names)) {
    colnames(U) <- bicop$names
  }

  U
}


#' @rdname bicop_methods
#' @param cond_var either `1` or `2`; `cond_var = 1` conditions on the first
#'    variable, `cond_var = 2` on the second.
#' @param inverse whether to compute the h-function or its inverse.
#' @export
hbicop <- function(u, cond_var, family, rotation, parameters, inverse = FALSE,
                   var_types = c("c", "c")) {
  assert_that(in_set(cond_var, 1:2), is.flag(inverse))
  bicop <- args2bicop(family, rotation, parameters, var_types)
  u <- if_vec_to_matrix(u)

  if (!inverse) {
    if (cond_var == 1) {
      return(bicop_hfunc1_cpp(u, bicop))
    } else {
      return(bicop_hfunc2_cpp(u, bicop))
    }
  } else {
    if (cond_var == 1) {
      return(bicop_hinv1_cpp(u, bicop))
    } else {
      return(bicop_hinv2_cpp(u, bicop))
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
#' par_to_ktau(bicop_dist("clayton", 0, 3))
#' par_to_ktau("clayton", 0, 3)
#'
#' ktau_to_par("clayton", 0.5)
#' ktau_to_par(bicop_dist("clayton", 0, 3), 0.5)
#' @name par_to_ktau
#' @rdname par_to_ktau
par_to_ktau <- function(family, rotation, parameters) {
  bicop <- args2bicop(family, rotation, parameters)
  bicop_par_to_tau_cpp(bicop)
}

#' @rdname par_to_ktau
#' @param tau Kendall's \eqn{\tau}.
#' @export
ktau_to_par <- function(family, tau) {
  bicop <- args2bicop(family)
  if (!(bicop$family %in% c(family_set_elliptical, family_set_nonparametric))) {
    bicop$rotation <- ifelse(tau > 0, 0, 90)
  }
  bicop_tau_to_par_cpp(bicop, tau)
}

#' Predictions and fitted values for a bivariate copula model
#'
#' Predictions of the density, distribution function, h-functions (with their
#' inverses) for a bivariate copula model.
#'
#' @name bicop_predict_and_fitted
#' @aliases predict.bicop fitted.bicop
#' @param object a `bicop` object.
#' @param newdata points where the fit shall be evaluated.
#' @param what what to predict, one of `"pdf"`, `"cdf"`, `"hfunc1"`, `"hfunc2"`,
#'   `"hinv1"`, `"hinv2"`.
#' @param ... unused.
#' @return `fitted()` and `logLik()` have return values similar to [dbicop()],
#' [pbicop()], and [hbicop()].
#'
#' @details `fitted()` can only be called if the model was fit with the
#'   `keep_data = TRUE` option.
#'
#'
#' ## Discrete variables
#' When at least one variable is discrete, mote than two columns are required
#' for `newdata`: the first \eqn{n \times 2} block contains realizations of
#' \eqn{F_{X_1}(x_1), F_{X_2}(x_2)}. The second \eqn{n \times 2} block contains
#' realizations of \eqn{F_{X_1}(x_1^-), F_{X_1}(x_1^-)}. The minus indicates a
#' left-sided limit of the cdf. For, e.g., an integer-valued variable, it holds
#' \eqn{F_{X_1}(x_1^-) = F_{X_1}(x_1 - 1)}. For continuous variables the left
#' limit and the cdf itself coincide. Respective columns can be omitted in the
#' second block.
#'
#'
#' @examples
#' # Simulate and fit a bivariate copula model
#' u <- rbicop(500, "gauss", 0, 0.5)
#' fit <- bicop(u, family = "par", keep_data = TRUE)
#'
#' # Predictions
#' all.equal(predict(fit, u, "hfunc1"), fitted(fit, "hfunc1"))
#' @rdname predict_bicop
#' @export
predict.bicop_dist <- function(object, newdata, what = "pdf", ...) {
  assert_that(in_set(what, what_allowed))
  newdata <- if_vec_to_matrix(newdata)
  switch(
    what,
    "pdf" = bicop_pdf_cpp(newdata, object),
    "cdf" = bicop_cdf_cpp(newdata, object),
    "hfunc1" = bicop_hfunc1_cpp(newdata, object),
    "hfunc2" = bicop_hfunc2_cpp(newdata, object),
    "hinv1" = bicop_hinv1_cpp(newdata, object),
    "hinv2" = bicop_hinv2_cpp(newdata, object)
  )
}

what_allowed <- c("pdf", "cdf", "hfunc1", "hfunc2", "hinv1", "hinv2")

#' @rdname predict_bicop
#' @export
fitted.bicop <- function(object, what = "pdf", ...) {
  if (is.null(object$data)) {
    stop("data have not been stored, use keep_data = TRUE when fitting.")
  }
  assert_that(in_set(what, what_allowed))
  switch(
    what,
    "pdf" = bicop_pdf_cpp(object$data, object),
    "cdf" = bicop_cdf_cpp(object$data, object),
    "hfunc1" = bicop_hfunc1_cpp(object$data, object),
    "hfunc2" = bicop_hfunc2_cpp(object$data, object),
    "hinv1" = bicop_hinv1_cpp(object$data, object),
    "hinv2" = bicop_hinv2_cpp(object$data, object)
  )
}

#' @importFrom stats logLik
#' @export
logLik.bicop <- function(object, ...) {
  structure(object$loglik, "df" = object$npars)
}

#' @export
print.bicop_dist <- function(x, ...) {
  x0 <- x
  if (x$family %in% setdiff(family_set_nonparametric, "indep")) {
    x$parameters <- paste0(round(x$npars, 2), sep = " d.f.")
  }
  cat("Bivariate copula ('bicop_dist'): ",
    "family = ", x$family,
    ", rotation = ", x$rotation,
    ", parameters = ", ifelse(length(x$parameters) > 1,
      paste(round(x$parameters, 2),
        collapse = ", "
      ),
      x$parameters
    ),
    ", var_types = ", paste(x$var_types, collapse = ","),
    sep = ""
  )
  cat("\n")
  invisible(x0)
}

#' @export
summary.bicop_dist <- function(object, ...) {
  print.bicop_dist(object, ...)
}

#' @export
print.bicop <- function(x, ...) {
  x0 <- x
  if (x$family %in% setdiff(family_set_nonparametric, "indep")) {
    pars_formatted <- paste0(round(x$npars, 2), sep = " d.f.")
  } else {
    pars_formatted <- paste(round(x$parameters, 2), collapse = ", ")
  }
  cat("Bivariate copula fit ('bicop'): ",
    "family = ", x$family,
    ", rotation = ", x$rotation,
    ", parameters = ", pars_formatted,
    ", var_types = ", paste(x$var_types, collapse = ","),
    "\n",
    sep = ""
  )

  invisible(x0)
}

#' @export
summary.bicop <- function(object, ...) {
  print.bicop(object, ...)
  cat("nobs =", object$nobs, "  ")

  info <- bicop_fit_info(object)
  cat("logLik =", round(info$logLik, 2), "  ")
  cat("npars =", round(info$npars, 2), "  ")
  cat("AIC =", round(info$AIC, 2), "  ")
  cat("BIC =", round(info$BIC, 2), "  ")
  attr(object, "info") <- info
  cat("\n")

  invisible(object)
}

#' @importFrom stats coef
#' @export
coef.bicop_dist <- function(object, ...) {
  object$parameters
}

bicop_fit_info <- function(bc) {
  ll <- logLik(bc)
  list(
    nobs = bc$nobs,
    logLik = ll[1],
    npars = attr(ll, "df"),
    AIC = -2 * ll[1] + 2 * attr(ll, "df"),
    BIC = -2 * ll[1] + log(bc$nobs) * attr(ll, "df")
  )
}
