#' Parameter uncertainty
#'
#' Covariance matrices and Wald confidence intervals for the estimated
#' pair-copula parameters.
#'
#' @name parameter_uncertainty
#' @aliases vcov.bicop_dist vcov.vinecop_dist confint.bicop_dist
#'   confint.vinecop_dist
#' @param object an object of class `"bicop"` or `"vinecop"`, or of class
#'   `"bicop_dist"` or `"vinecop_dist"` together with `newdata`.
#' @param newdata copula data to evaluate the derivatives at. Defaults to the
#'   data stored in `object`, which requires `keep_data = TRUE` at fitting
#'   time.
#' @param step_wise if `TRUE` (default), derivatives of the step-wise
#'   estimating equations; if `FALSE`, of the joint log-likelihood. See
#'   [scores()].
#' @param cores number of cores to use; if larger than one, computations are
#'   done in parallel on `cores` batches.
#' @param parm a specification of which parameters to give intervals for: a
#'   vector of names or of positions. Defaults to all of them.
#' @param level the confidence level.
#' @param ... passed on to `vcov()`.
#'
#' @details
#' The estimator solves \eqn{\bar s(\hat\theta) = 0} with
#' \eqn{\bar s(\theta) = n^{-1} \sum_i s_i(\theta)}, so it is an M-estimator.
#' Writing
#' \deqn{A = -\partial \bar s(\theta) / \partial \theta^\top, \qquad
#'       B = n^{-1} \sum_i s_i s_i^\top,}
#' `vcov()` returns the sandwich \eqn{A^{-1} B A^{-\top} / n} of White (1982).
#' The simpler \eqn{A^{-1}/n} would need the information equality \eqn{A = B},
#' which holds only at a maximum-likelihood estimator and so not for the
#' default step-wise objective below; the sandwich is correct either way, so it
#' is the only form offered.
#'
#' ## Which objective is differentiated
#'
#' `step_wise` selects the estimating equation, and it changes the nature of
#' \eqn{A}.
#'
#' With `step_wise = FALSE`, \eqn{\bar s} is the gradient of the joint
#' log-likelihood and \eqn{A} is minus its averaged Hessian: a symmetric
#' matrix, and the classical setting for both formulas.
#'
#' With `step_wise = TRUE` (the default), the pseudo-observations entering each
#' pair copula are held fixed. This is the estimating equation that sequential,
#' tree-by-tree estimation actually solves, but it is not the gradient of any
#' objective: \eqn{A} is a Jacobian rather than a Hessian, and it is block
#' triangular rather than symmetric. The sandwich formula still applies --
#' which is why it is written with \eqn{A^{-\top}} on the right.
#'
#' ## Restrictions
#'
#' No derivatives are implemented for nonparametric pair copulas, so a model
#' containing one is refused. For models with discrete variables the backend
#' differentiates by central finite differences rather than in closed form; the
#' result is correct but carries the accuracy of a finite-difference
#' approximation.
#'
#' ## Why there is no method for vine distributions
#'
#' These functions treat the copula data as observed, so they quantify
#' uncertainty in the pair-copula parameters *given* the margins. A model
#' fitted by [vine()] estimates its margins too, and the second stage inherits
#' the first stage's sampling variability: intervals that ignore it are too
#' short. In a small experiment (three-dimensional Gaussian D-vine,
#' equicorrelation 0.5, `n = 1000`, 1000 replications) nominal 95% intervals
#' attained 0.956 when the true margins were used and 0.929 when they were
#' estimated by ranks.
#'
#' Correcting for it needs the influence function of the marginal estimator,
#' and [fit_margin()] is a user-supplied black box, so there is nothing to
#' differentiate. We therefore provide no method for `"vine"` objects rather
#' than a silently anticonservative one. To quantify uncertainty for a complete
#' vine distribution, resample:
#'
#' ```
#' boot <- replicate(B, {
#'   idx <- sample(nrow(x), replace = TRUE)
#'   fit <- vine(x[idx, ], margins_controls = mc, copula_controls = cc)
#'   unlist(get_all_parameters(fit$copula))
#' })
#' apply(boot, 1, quantile, c(0.025, 0.975))
#' ```
#'
#' @return `vcov()` returns a covariance matrix with one row and column per
#' estimated pair-copula parameter, ordered by tree, then edge, then parameter
#' within the pair copula, with dimnames giving that position. `confint()`
#' returns a matrix with columns for the lower and upper endpoints.
#'
#' @seealso [scores()], [hessian()], [vinecop()], [vine()]
#'
#' @references
#' White H (1982). "Maximum Likelihood Estimation of Misspecified Models."
#' *Econometrica*, **50**(1), 1--25. \doi{10.2307/1912526}
#'
#' Stoeber J, Schepsmeier U (2013). "Estimating Standard Errors in Regular Vine
#' Copula Models." *Computational Statistics*, **28**(6), 2679--2707.
#' \doi{10.1007/s00180-013-0423-8}
#'
#' @examples
#' # bivariate copula ------------------------------------------
#' u <- rbicop(500, "gaussian", 0, 0.5)
#' fit <- bicop(u, family_set = "parametric", keep_data = TRUE)
#' sqrt(diag(vcov(fit)))
#' confint(fit)
#'
#' # vine copula -----------------------------------------------
#' u <- pseudo_obs(matrix(rnorm(500 * 3), 500, 3))
#' fit <- vinecop(u, family_set = "parametric", keep_data = TRUE)
#' vcov(fit)
#' confint(fit, level = 0.9)
#'
#' # derivatives of the joint log-likelihood instead of the step-wise ones
#' sqrt(diag(vcov(fit, step_wise = FALSE)))
#'
#' @importFrom stats vcov confint qnorm
NULL

## ---------------------------------------------------------------------------
## helpers
## ---------------------------------------------------------------------------

#' Names for the columns of `scores()`, in (tree, edge, parameter) order.
#' @noRd
par_labels <- function(object) {
  if (inherits(object, "bicop_dist")) {
    np <- NROW(object$parameters)
    if (np <= 1) {
      return(object$family)
    }
    return(paste0(object$family, ".par", seq_len(np)))
  }

  pcs <- object$pair_copulas
  mat <- get_matrix(object)
  d <- dim(object)[1]
  labs <- character(0)
  for (t in seq_along(pcs)) {
    for (e in seq_along(pcs[[t]])) {
      pc <- pcs[[t]][[e]]
      np <- NROW(pc$parameters)
      if (np == 0) {
        next
      }
      pair <- paste0(mat[d - e + 1, e], ",", mat[t, e])
      cond <- if (t > 1) {
        paste0(";", paste(mat[seq_len(t - 1), e], collapse = ","))
      } else {
        ""
      }
      base <- paste0("T", t, ".", pair, cond)
      labs <- c(labs, if (np == 1) base else paste0(base, ".par", seq_len(np)))
    }
  }
  labs
}

#' @noRd
vcov_check_supported <- function(object) {
  fams <- if (inherits(object, "vinecop_dist")) {
    unlist(get_all_families(object))
  } else {
    object$family
  }
  if (any(fams == "tll")) {
    stop(
      "vcov() requires analytic derivatives, which are unavailable for ",
      "nonparametric ('tll') pair copulas. Refit with a parametric ",
      "`family_set`, or use scores()/hessian() directly.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' @noRd
vcov_data <- function(object, newdata) {
  if (!is.null(newdata)) {
    return(if_vec_to_matrix(newdata, dim(object)[1] == 1))
  }
  if (is.null(object$data)) {
    stop(
      "no data available: either fit with `keep_data = TRUE` or supply ",
      "`newdata`.",
      call. = FALSE
    )
  }
  object$data
}

#' @noRd
vcov_assemble <- function(u, object, step_wise, cores, is_bicop) {
  s <- if (is_bicop) {
    scores(u, object, cores = cores)
  } else {
    scores(u, object, step_wise = step_wise, cores = cores)
  }
  H <- if (is_bicop) {
    hessian(u, object, cores = cores)
  } else {
    hessian(u, object, step_wise = step_wise, cores = cores)
  }

  n <- NROW(s)
  if (NCOL(s) == 0) {
    return(matrix(numeric(), 0, 0))
  }

  ## hessian() returns d sbar / d theta, so the M-estimation bread is its
  ## negative transpose. Under step_wise = TRUE this matrix is block
  ## triangular, not symmetric: symmetrizing it would discard the cross-tree
  ## propagation terms, so we must not.
  A <- -t(H)

  Ainv <- tryCatch(solve(A), error = function(e) {
    stop(
      "the matrix of derivatives of the estimating equations is numerically ",
      "singular, so no covariance matrix can be formed. This usually means a ",
      "parameter sits at a boundary of its family's parameter space.",
      call. = FALSE
    )
  })

  B <- crossprod(s) / n
  V <- Ainv %*% B %*% t(Ainv) / n
  ## the sandwich is symmetric in exact arithmetic
  V <- (V + t(V)) / 2

  labs <- par_labels(object)
  if (length(labs) == ncol(V)) {
    dimnames(V) <- list(labs, labs)
  }
  V
}

## ---------------------------------------------------------------------------
## methods
## ---------------------------------------------------------------------------

#' @rdname parameter_uncertainty
#' @export
vcov.vinecop_dist <- function(
  object,
  newdata = NULL,
  step_wise = TRUE,
  cores = 1,
  ...
) {
  vcov_check_supported(object)
  u <- vcov_data(object, newdata)
  vcov_assemble(u, object, step_wise, cores, is_bicop = FALSE)
}

#' @rdname parameter_uncertainty
#' @export
vcov.bicop_dist <- function(
  object,
  newdata = NULL,
  step_wise = TRUE,
  cores = 1,
  ...
) {
  vcov_check_supported(object)
  u <- vcov_data(object, newdata)
  vcov_assemble(u, object, step_wise, cores, is_bicop = TRUE)
}

#' @rdname parameter_uncertainty
#' @export
confint.vinecop_dist <- function(object, parm, level = 0.95, ...) {
  confint_from_vcov(object, parm, level, ...)
}

#' @rdname parameter_uncertainty
#' @export
confint.bicop_dist <- function(object, parm, level = 0.95, ...) {
  confint_from_vcov(object, parm, level, ...)
}

#' @noRd
confint_from_vcov <- function(object, parm, level, ...) {
  assert_that(is.number(level), level > 0, level < 1)
  V <- vcov(object, ...)
  est <- unlist(get_all_parameters_flat(object))
  se <- sqrt(diag(V))
  if (length(est) != length(se)) {
    stop(
      "parameter vector and covariance matrix have inconsistent lengths; ",
      "this is a bug, please report it.",
      call. = FALSE
    )
  }
  names(est) <- rownames(V)

  a <- (1 - level) / 2
  z <- qnorm(1 - a)
  ci <- cbind(est - z * se, est + z * se)
  colnames(ci) <- paste(format(100 * c(a, 1 - a), trim = TRUE, digits = 3), "%")
  rownames(ci) <- names(est)

  if (!missing(parm)) {
    ci <- ci[parm, , drop = FALSE]
  }
  ci
}

#' Estimated parameters in the order used by `scores()`.
#' @noRd
get_all_parameters_flat <- function(object) {
  if (inherits(object, "bicop_dist")) {
    return(as.vector(object$parameters))
  }
  out <- numeric(0)
  for (tree in object$pair_copulas) {
    for (pc in tree) {
      out <- c(out, as.vector(pc$parameters))
    }
  }
  out
}
