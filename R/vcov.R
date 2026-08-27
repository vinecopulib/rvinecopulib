#' Parameter uncertainty
#'
#' Covariance matrices and Wald confidence intervals for the estimated
#' pair-copula parameters.
#'
#' @name parameter_uncertainty
#' @aliases vcov.bicop_dist vcov.vinecop_dist vcov.vine_dist confint.bicop_dist
#'   confint.vinecop_dist confint.vine_dist
#' @param object an object of class `"bicop"`, `"vinecop"` or `"vine"`; or of
#'   class `"bicop_dist"`, `"vinecop_dist"` or `"vine_dist"` together with
#'   `newdata`.
#' @param newdata data to evaluate the derivatives at. Copula data for the
#'   copula classes, original-scale data for `"vine_dist"`. Defaults to the data
#'   stored in `object`, which requires `keep_data = TRUE` at fitting time.
#' @param type `"sandwich"` (the default) or `"model"`; see *Details*.
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
#' `type = "sandwich"` returns \eqn{A^{-1} B A^{-\top} / n}, the
#' misspecification-robust form of White (1982), and `type = "model"` returns
#' \eqn{A^{-1} / n}.
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
#' which is why it is written with \eqn{A^{-\top}} on the right -- but the
#' information equality does not, so `type = "model"` is refused.
#'
#' ## Restrictions
#'
#' No derivatives are implemented for nonparametric pair copulas, so a model
#' containing one is refused. For models with discrete variables the backend
#' differentiates by central finite differences rather than in closed form; the
#' result is correct but carries the accuracy of a finite-difference
#' approximation.
#'
#' Neither variant accounts for estimation of the margins: both treat the
#' copula data as observed. For vine distributions fitted with [vine()] this is
#' the usual two-stage approximation, and it makes the intervals somewhat too
#' short.
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
#' # the model-based form requires the joint objective
#' vcov(fit, type = "model", step_wise = FALSE)
#'
#' # vine distribution on the data scale -----------------------
#' x <- data.frame(a = rnorm(500), b = rnorm(500))
#' fit <- vine(x, copula_controls = list(family_set = "parametric"),
#'             keep_data = TRUE)
#' confint(fit)
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
vcov_assemble <- function(u, object, type, step_wise, cores, is_bicop) {
  type <- match.arg(type, c("sandwich", "model"))
  if (type == "model" && !is_bicop && isTRUE(step_wise)) {
    stop(
      "type = \"model\" relies on the information equality, which holds for a ",
      "maximum-likelihood estimator but not for the step-wise one. Use ",
      "type = \"sandwich\", or step_wise = FALSE to differentiate the joint ",
      "log-likelihood.",
      call. = FALSE
    )
  }
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

  V <- if (type == "model") {
    Ainv / n
  } else {
    B <- crossprod(s) / n
    Ainv %*% B %*% t(Ainv) / n
  }
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
  type = c("sandwich", "model"),
  step_wise = TRUE,
  cores = 1,
  ...
) {
  vcov_check_supported(object)
  u <- vcov_data(object, newdata)
  vcov_assemble(u, object, type, step_wise, cores, is_bicop = FALSE)
}

#' @rdname parameter_uncertainty
#' @export
vcov.bicop_dist <- function(
  object,
  newdata = NULL,
  type = c("sandwich", "model"),
  step_wise = TRUE,
  cores = 1,
  ...
) {
  vcov_check_supported(object)
  u <- vcov_data(object, newdata)
  vcov_assemble(u, object, type, step_wise, cores, is_bicop = TRUE)
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

#' @rdname parameter_uncertainty
#' @export
vcov.vine_dist <- function(
  object,
  newdata = NULL,
  type = c("sandwich", "model"),
  step_wise = TRUE,
  cores = 1,
  ...
) {
  vcov_check_supported(object$copula)
  u <- vine_copula_data(object, newdata)
  vcov_assemble(u, object$copula, type, step_wise, cores, is_bicop = FALSE)
}

#' @rdname parameter_uncertainty
#' @export
confint.vine_dist <- function(object, parm, level = 0.95, ...) {
  assert_that(is.number(level), level > 0, level < 1)
  V <- vcov(object, ...)
  est <- get_all_parameters_flat(object$copula)
  se <- sqrt(diag(V))
  names(est) <- rownames(V)
  a <- (1 - level) / 2
  ci <- cbind(est - qnorm(1 - a) * se, est + qnorm(1 - a) * se)
  colnames(ci) <- paste(format(100 * c(a, 1 - a), trim = TRUE, digits = 3), "%")
  rownames(ci) <- names(est)
  if (!missing(parm)) {
    ci <- ci[parm, , drop = FALSE]
  }
  ci
}

#' Copula-scale data for a fitted vine distribution.
#'
#' A `vine` stores data on the original scale, so the marginal transformation
#' has to be reapplied before the copula derivatives can be evaluated.
#' @noRd
vine_copula_data <- function(object, newdata) {
  dat <- if (!is.null(newdata)) newdata else object$data
  if (is.null(dat)) {
    stop(
      "no data available: either fit with `keep_data = TRUE` or supply ",
      "`newdata`.",
      call. = FALSE
    )
  }
  compute_pseudo_obs(expand_factors(as.data.frame(dat)), object)
}
