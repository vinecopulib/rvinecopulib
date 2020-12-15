#' Vine based distributions
#'
#' Density, distribution function and random generation
#' for the vine based distribution.
#'
#' @name vine_distributions
#' @aliases dvine pvine rvine dvine_dist pvine_dist rvine_dist
#' @param x evaluation points, either a length d vector or a d-column matrix,
#'   where d is the number of variables in the vine.
#' @param vine an object of class `"vine_dist"`.
#' @param cores number of cores to use; if larger than one, computations are
#'   done in parallel on `cores` batches .
#' @details
#' See [vine] for the estimation and construction of vine models.
#' Here, the density, distribution function and random generation
#' for the vine distributions are standard.
#'
#' The functions are based on [dvinecop()], [pvinecop()] and [rvinecop()] for
#' [vinecop] objects, and either [kde1d::dkde1d()], [kde1d::pkde1d()] and
#' [kde1d::qkde1d()] for estimated vines (i.e., output of [vine()]), or the
#' standard *d/p/q-xxx* from [stats::Distributions] for custom vines
#' (i.e., output of [vine_dist()]).
#' @return
#' `dvine()` gives the density, `pvine()` gives the distribution function,
#' and `rvine()` generates random deviates.
#'
#' The length of the result is determined by `n` for `rvine()`, and
#' the number of rows in `u` for the other functions.
#'
#' The `vine` object is recycled to the length of the
#' result.
#' @examples
#' # specify pair-copulas
#' bicop <- bicop_dist("bb1", 90, c(3, 2))
#' pcs <- list(
#'   list(bicop, bicop), # pair-copulas in first tree
#'   list(bicop) # pair-copulas in second tree
#' )
#'
#' # set up vine copula model
#' mat <- rvine_matrix_sim(3)
#' vc <- vine_dist(list(distr = "norm"), pcs, mat)
#'
#' # simulate from the model
#' x <- rvine(200, vc)
#' pairs(x)
#'
#' # evaluate the density and cdf
#' dvine(x[1, ], vc)
#' pvine(x[1, ], vc)
#' @rdname vine_methods
#' @export
dvine <- function(x, vine, cores = 1) {
  stopifnot(inherits(vine, "vine_dist"))
  if (NCOL(x) == 1) {
    x <- t(x)
  }

  x <- expand_factors(x)
  if (!is.null(vine$names)) {
    x <- x[, vine$names, drop = FALSE]
  }

  # prepare marginals if only one is specified
  d <- ncol(x)
  if (!inherits(vine, "vine") & depth(vine$margins) == 1) {
    vine$margins <- replicate(d, vine$margins, simplify = FALSE)
  }

  ## evaluate marginal densities
  margvals <- dpq_marg(x, vine, "d")

  if (!is.null(vine$copula)) {
    u <- compute_pseudo_obs(x, vine)
    vinevals <- dvinecop(u, vine$copula, cores)
  } else {
    vinevals <- rep(1, nrow(x))
  }

  ## final density estimate is product of marginals and copula density
  apply(cbind(margvals, vinevals), 1, prod)
}

#' @rdname vine_methods
#' @param n_mc number of samples used for quasi Monte Carlo integration.
#' @export
pvine <- function(x, vine, n_mc = 10^4, cores = 1) {
  stopifnot(inherits(vine, "vine_dist"))

  if (NCOL(x) == 1) {
    x <- t(x)
  }
  x <- expand_factors(x)
  if (!is.null(vine$names)) {
    x <- x[, vine$names, drop = FALSE]
  }

  # prepare marginals if only one is specified
  if (!inherits(vine, "vine") & depth(vine$margins) == 1) {
    vine$margins <- replicate(ncol(x), vine$margins, simplify = FALSE)
  }

  # PIT to copula data
  u <- compute_pseudo_obs(x, vine)

  # Evaluate copula if needed
  if (!is.null(vine$copula)) {
    vals <- pvinecop(u, vine$copula, n_mc, cores)
  } else {
    vals <- apply(u, 1, prod)
  }

  vals
}

#' @rdname vine_methods
#' @param n number of observations.
#' @param qrng if `TRUE`, generates quasi-random numbers using the multivariate
#' Generalized Halton sequence up to dimension 300 and the Generalized Sobol
#' sequence in higher dimensions (default `qrng = FALSE`).
#' @export
rvine <- function(n, vine, qrng = FALSE, cores = 1) {
  assert_that(inherits(vine, "vine_dist"), is.flag(qrng))

  # simulate copula data
  U <- rvinecop(n, vine$copula, qrng, cores)

  # prepare marginals if only one is specified
  if (!inherits(vine, "vine") & depth(vine$margins) == 1) {
    vine$margins <- replicate(dim(vine)[1], vine$margins, simplify = FALSE)
  }

  # use quantile transformation for marginals
  X <- dpq_marg(U, vine, "q")
  colnames(X) <- vine$names
  X
}

#' @export
print.vine_dist <- function(x, ...) {
  cat(dim(x)[1], "-dimensional vine distribution model ('vine_dist')", sep = "")
  print_truncation_info(x$copula)
  invisible(x)
}

#' @export
summary.vine_dist <- function(object, ...) {
  list(
    margins = get_vine_dist_margin_summary(object),
    copula = summary(object$copula)
  )
}

get_vine_dist_margin_summary <- function(vd) {
  margins <- vd$margins
  if (length(margins) == 1) {
    margins <- rep(list(margins), dim(vd$copula)[1])
  }
  df <- data.frame(
    margin = seq_along(margins),
    distr = sapply(margins, function(x) x$distr)
  )
  class(df) <- c("summary_df", class(df))
  df
}

#' Predictions and fitted values for a vine copula model
#'
#' Predictions of the density and distribution function
#' for a vine copula model.
#'
#' @name vine_predict_and_fitted
#' @aliases fitted.vine predict.vine
#' @param object a `vine` object.
#' @param newdata points where the fit shall be evaluated.
#' @param what what to predict, either `"pdf"` or `"cdf"`.
#' @param n_mc number of samples used for quasi Monte Carlo integration when
#'    `what = "cdf"`.
#' @param cores number of cores to use; if larger than one, computations are
#'   done in parallel on `cores` batches .
#' @param ... unused.
#'
#' @return
#' `fitted()` and `predict()` have return values similar to [dvine()]
#' and [pvine()].
#' @export
#' @rdname predict_vine
#' @examples
#' x <- sapply(1:5, function(i) rnorm(50))
#' fit <- vine(x, copula_controls = list(family_set = "par"), keep_data = TRUE)
#' all.equal(predict(fit, x), fitted(fit), check.environment = FALSE)
predict.vine <- function(object, newdata, what = "pdf", n_mc = 10^4,
                         cores = 1, ...) {
  stopifnot(what %in% c("pdf", "cdf"))
  switch(
    what,
    "pdf" = dvine(newdata, object, cores),
    "cdf" = pvine(newdata, object, n_mc, cores)
  )
}

#' @rdname predict_vine
#' @export
fitted.vine <- function(object, what = "pdf", n_mc = 10^4, cores = 1, ...) {
  if (all(is.na(object$data))) {
    stop("data have not been stored, use keep_data = TRUE when fitting.")
  }
  stopifnot(what %in% c("pdf", "cdf"))
  switch(
    what,
    "pdf" = dvine(object$data, object, cores),
    "cdf" = pvine(object$data, object, n_mc, cores)
  )
}

#' @export
logLik.vine <- function(object, ...) {
  structure(object$loglik, "df" = object$npars)
}

#' @export
print.vine <- function(x, ...) {
  cat(dim(x)[1], "-dimensional vine distribution fit ('vine')", sep = "")
  print_truncation_info(x$copula)
  print_fit_info(x)
  invisible(x)
}

#' @export
summary.vine <- function(object, ...) {
  list(
    margins = get_vine_margin_summary(object),
    copula = summary(object$copula)
  )
}

get_vine_margin_summary <- function(object) {
  capture.output(info <- sapply(object$margins, summary))
  info <- as.data.frame(t(info))
  info <- cbind(
    data.frame(margin = seq_len(nrow(info)), name = object$names),
    info
  )
  class(info) <- c("summary_df", "data.frame")
  info
}


dpq_marg <- function(x, vine, what = "p") {
  d <- ncol(x)
  res <- lapply(seq_len(d),
                function(i) eval_one_dpq(x[, i], vine$margins[[i]], what))
  do.call(cbind, res)
}

eval_one_dpq <- function(x, margin, what = "p") {
  if (inherits(margin, "kde1d")) {
    dpq <- switch(what,
                  p = pkde1d(x, margin),
                  d = dkde1d(x, margin),
                  q = qkde1d(x, margin))
  } else {
    par <- margin[names(margin) != "distr"]
    par[[length(par) + 1]] <- x
    names(par)[[length(par)]] <- switch(what,
                                        p = "q",
                                        d = "x",
                                        q = "p")
    dpq <- do.call(get(paste0(what, margin$distr)), par)
  }
  dpq
}

collate_u <- function(x) {
  if (!all(is.na(x$data))) {
    u <- dpq_marg(x$data, x)
    x$copula$data <- u
  }
  x
}

#' @export
dim.vine_dist <- function(x) {
  dim(x$copula)
}
