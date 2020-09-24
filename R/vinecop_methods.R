#' Vine copula distributions
#'
#' Density, distribution function and random generation for the vine copula
#' distribution.
#'
#' @name vinecop_distributions
#' @aliases dvinecop pvinecop rvinecop dvinecop_dist pvinecop_dist rvinecop_dist
#' @param u matrix of evaluation points; must contain at least d columns, where
#'   d is the number of variables in the vine. More columns are required for
#'   discrete models, see *Details*.
#' @param vinecop an object of class `"vinecop_dist"`.
#' @param cores number of cores to use; if larger than one, computations are
#'   done in parallel on `cores` batches .
#' @details See [vinecop()] for the estimation and construction of vine copula
#' models.
#'
#' The copula density is defined as joint density divided by marginal
#' densities, irrespective of variable types.
#'
#' ## Discrete variables
#'
#' When at least one variable is discrete, two types of
#' "observations" are required in `u`: the first \eqn{n \; x \; d} block
#' contains realizations of \eqn{F_{X_j}(X_j)}. The second \eqn{n \; x \; d}
#' block contains realizations of \eqn{F_{X_j}(X_j^-)}. The minus indicates a
#' left-sided limit of the cdf. For, e.g., an integer-valued variable, it holds
#' \eqn{F_{X_j}(X_j^-) = F_{X_j}(X_j - 1)}. For continuous variables the left
#' limit and the cdf itself coincide. Respective columns can be omitted in the
#' second block.
#'
#' @return
#' `dvinecop()` gives the density, `pvinecop()` gives the distribution function,
#' and `rvinecop()` generates random deviates.
#'
#' The length of the result is determined by `n` for `rvinecop()`, and
#' the number of rows in `u` for the other functions.
#'
#' The `vinecop` object is recycled to the length of the
#' result.
#'
#' @seealso [vinecop_dist()], [vinecop()], [plot.vinecop()], [contour.vinecop()]
#'
#' @examples
#' ## simulate dummy data
#' x <- rnorm(30) * matrix(1, 30, 5) + 0.5 * matrix(rnorm(30 * 5), 30, 5)
#' u <- pseudo_obs(x)
#'
#' ## fit a model
#' vc <- vinecop(u, family = "clayton")
#'
#' # simulate from the model
#' u <- rvinecop(100, vc)
#' pairs(u)
#'
#' # evaluate the density and cdf
#' dvinecop(u[1, ], vc)
#' pvinecop(u[1, ], vc)
#'
#' ## Discrete models
#' vc$var_types <- rep("d", 5)  # convert model to discrete
#'
#' # with discrete data we need two types of observations (see Details)
#' x <- qpois(u, 1)  # transform to Poisson margins
#' u_disc <- cbind(ppois(x, 1), ppois(x - 1, 1))
#'
#' dvinecop(u_disc[1:5, ], vc)
#' pvinecop(u_disc[1:5, ], vc)
#'
#' # simulated data always has uniform margins
#' pairs(rvinecop(200, vc))
#' @rdname vinecop_methods
#' @export
dvinecop <- function(u, vinecop, cores = 1) {
  assert_that(inherits(vinecop, "vinecop_dist"))
  u <- if_vec_to_matrix(u, dim(vinecop)[1] == 1)
  vinecop_pdf_cpp(u, vinecop, cores)
}

#' @rdname vinecop_methods
#' @param n_mc number of samples used for quasi Monte Carlo integration.
#' @importFrom assertthat is.count
#' @export
pvinecop <- function(u, vinecop, n_mc = 10^4, cores = 1) {
  assert_that(
    inherits(vinecop, "vinecop_dist"),
    is.number(n_mc), is.count(cores)
  )
  u <- if_vec_to_matrix(u, dim(vinecop)[1] == 1)
  vinecop_cdf_cpp(as.matrix(u), vinecop, n_mc, cores, get_seeds())
}

#' @rdname vinecop_methods
#' @param n number of observations.
#' @param qrng if `TRUE`, generates quasi-random numbers using the multivariate
#' Generalized Halton sequence up to dimension 300 and the Generalized Sobol
#' sequence in higher dimensions (default `qrng = FALSE`).
#' @export
rvinecop <- function(n, vinecop, qrng = FALSE, cores = 1) {
  assert_that(
    is.number(n),
    inherits(vinecop, "vinecop_dist"),
    is.flag(qrng),
    is.number(cores)
  )

  U <- vinecop_sim_cpp(vinecop, n, qrng, cores, get_seeds())
  if (!is.null(vinecop$names)) {
    colnames(U) <- vinecop$names
  }

  U
}

#' @export
print.vinecop_dist <- function(x, ...) {
  cat(dim(x)[1], "-dimensional vine copula model ('vinecop_dist')", sep = "")
  print_truncation_info(x)
  invisible(x)
}

#' @importFrom utils capture.output
#' @export
summary.vinecop_dist <- function(object, ...) {
  mat <- as_rvine_matrix(get_structure(object))
  d <- dim(object)[1]
  n_trees <- dim(object)[2]
  n_pcs <- length(unlist(object$pair_copulas, recursive = FALSE))
  mdf <- as.data.frame(matrix(NA, n_pcs, 10))
  names(mdf) <- c(
    "tree", "edge",
    "conditioned", "conditioning", "var_types",
    "family", "rotation", "parameters", "df", "tau"
  )
  k <- 1
  for (t in seq_len(n_trees)) {
    for (e in seq_len(d - t)) {
      mdf$tree[k] <- t
      mdf$edge[k] <- e
      mdf$conditioned[k] <- list(c(mat[d - e + 1, e], mat[t, e]))
      mdf$conditioning[k] <- list(mat[rev(seq_len(t - 1)), e])
      pc <- object$pair_copulas[[t]][[e]]
      mdf$var_types[k] <- paste(pc$var_types, collapse = ",")
      mdf$family[k] <- pc$family
      mdf$rotation[k] <- pc$rotation
      mdf$parameters[k] <- list(pc$parameters)
      if (pc$family %in% setdiff(family_set_nonparametric, "indep"))
        mdf$parameters[k] <- list("[30x30 grid]")
      mdf$df[k] <- pc$npars
      mdf$tau[k] <- par_to_ktau(pc)
      k <- k + 1
    }
  }
  class(mdf) <- c("summary_df", class(mdf))
  mdf
}

#' Predictions and fitted values for a vine copula model
#'
#' Predictions of the density and distribution function
#' for a vine copula model.
#'
#' @name vinecop_predict_and_fitted
#' @aliases fitted.vinecop predict.vinecop
#' @param object a `vinecop` object.
#' @param newdata points where the fit shall be evaluated.
#' @param what what to predict, either `"pdf"` or `"cdf"`.
#' @param n_mc number of samples used for quasi Monte Carlo integration when
#'    `what = "cdf"`.
#' @param cores number of cores to use; if larger than one, computations are
#'   done in parallel on `cores` batches.
#' @param ... unused.
#'
#' @details `fitted()` can only be called if the model was fit with the
#'    `keep_data = TRUE` option.
#'
#' ## Discrete variables
#'
#' When at least one variable is discrete, two types of
#' "observations" are required in `newdata`: the first \eqn{n \; x \; d} block
#' contains realizations of \eqn{F_{X_j}(X_j)}. The second \eqn{n \; x \; d}
#' block contains realizations of \eqn{F_{X_j}(X_j^-)}. The minus indicates a
#' left-sided limit of the cdf. For, e.g., an integer-valued variable, it holds
#' \eqn{F_{X_j}(X_j^-) = F_{X_j}(X_j - 1)}. For continuous variables the left
#' limit and the cdf itself coincide. Respective columns can be omitted in the
#' second block.
#'
#' @return
#' `fitted()` and `predict()` have return values similar to [dvinecop()]
#' and [pvinecop()].
#' @export
#' @rdname predict_vinecop
#' @examples
#' u <- sapply(1:5, function(i) runif(50))
#' fit <- vinecop(u, family = "par", keep_data = TRUE)
#' all.equal(predict(fit, u), fitted(fit))
predict.vinecop <- function(object, newdata, what = "pdf", n_mc = 10^4,
                            cores = 1, ...) {
  assert_that(
    in_set(what, c("pdf", "cdf")),
    is.number(n_mc),
    is.number(cores), cores > 0
  )
  newdata <- if_vec_to_matrix(newdata, dim(object)[1] == 1)
  switch(
    what,
    "pdf" = vinecop_pdf_cpp(newdata, object, cores),
    "cdf" = vinecop_cdf_cpp(newdata, object, n_mc, cores, get_seeds())
  )
}

#' @rdname predict_vinecop
#' @export
fitted.vinecop <- function(object, what = "pdf", n_mc = 10^4, cores = 1, ...) {
  if (is.null(object$data)) {
    stop("data have not been stored, use keep_data = TRUE when fitting.")
  }
  assert_that(
    in_set(what, c("pdf", "cdf")),
    is.number(n_mc),
    is.number(cores), cores > 0
  )
  switch(
    what,
    "pdf" = vinecop_pdf_cpp(object$data, object, cores),
    "cdf" = vinecop_cdf_cpp(object$data, object, n_mc, cores, get_seeds())
  )
}

#' @export
logLik.vinecop <- function(object, ...) {
  structure(object$loglik, "df" = object$npars)
}

#' Modified vine copula Bayesian information criterion (mBICv)
#'
#' Calculates the modified vine copula Bayesian information criterion.
#'
#' The modified vine copula Bayesian information criterion (mBICv) is defined as
#'
#' \deqn{BIC = -2 loglik +  \nu log(n) - 2 \sum_{t=1}^{d - 1} (q_t log(\psi_0^t)
#' - (d - t - q_t) log(1 - \psi_0^t)) }
#'
#' where \eqn{\mathrm{loglik}} is the log-likelihood and \eqn{\nu} is the
#' (effective) number of parameters of the model, \eqn{t} is the tree level
#' \eqn{\psi_0} is the prior probability of having a non-independence copula and
#' \eqn{q_t} is the number of non-independence copulas in tree \eqn{t}. The
#' mBICv is a consistent model selection criterion for parametric sparse vine
#' copula models.
#'
#' @param object a fitted `vinecop` object.
#' @param psi0 baseline prior probability of a non-independence copula.
#' @param newdata optional; a new data set.
#'
#' @references Nagler, T., Bumann, C., Czado, C. (2019). Model selection for
#'   sparse high-dimensional vine copulas with application to portfolio risk.
#'   *Journal of Multivariate Analysis, in press*
#'   (\url{https://arxiv.org/pdf/1801.09739.pdf})
#'
#' @export mBICV
#' @examples
#' u <- sapply(1:5, function(i) runif(50))
#' fit <- vinecop(u, family = "par", keep_data = TRUE)
#' mBICV(fit, 0.9) # with a 0.9 prior probability of a non-independence copula
#' mBICV(fit, 0.1) # with a 0.1 prior probability of a non-independence copula
mBICV <- function(object, psi0 = 0.9, newdata = NULL) {
  assert_that(inherits(object, "vinecop_dist"), is.number(psi0))
  ll <- ifelse(is.null(newdata),
    object$loglik,
    sum(log(dvinecop(newdata, object)))
  )
  - 2 * ll + compute_mBICV_penalty(object, psi0)
}

compute_mBICV_penalty <- function(object, psi0) {
  d <- dim(object)[1]
  smr <- summary(object)
  q_m <- tapply(smr$family, smr$tree, function(x) sum(x == "indep"))
  q_m <- c(q_m, rep(0, d - 1 - length(q_m)))
  m_seq <- seq_len(d - 1)
  pen <- object$npars * log(object$nobs)
  pen - 2 * sum(
    q_m * log(psi0^m_seq) + (d - seq_len(d - 1) - q_m) * log(1 - psi0^m_seq)
  )
}

#' @export
print.vinecop <- function(x, ...) {
  cat(dim(x)[1], "-dimensional vine copula fit ('vinecop')", sep = "")
  print_truncation_info(x)
  print_fit_info(x)
  invisible(x)
}


#' @export
summary.vinecop <- function(object, ...) {
  mdf <- summary.vinecop_dist(object)

  d <- dim(object)[1]
  trunc_lvl <- dim(object)[2]
  k <- 1
  for (t in seq_len(trunc_lvl)) {
    for (e in seq_len(d - t)) {
      mdf$loglik[k] <- object$pair_copulas[[t]][[e]]$loglik
      k <- k + 1
    }
  }

  mdf
}

#' @export
dim.vinecop_dist <- function(x) {
  dim(x$structure)
}
