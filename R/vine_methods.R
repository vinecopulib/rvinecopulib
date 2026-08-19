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
#' @param x_cond optional conditioning values for `rvine()` on the original
#'   data scale. A vector or one-row object is repeated `n` times;
#'   alternatively, supply an `n`-row matrix or data frame for
#'   observation-specific conditioning values. If `NULL`, `rvine()` performs
#'   unconditional simulation.
#' @param conditioning_set variable indices or names corresponding to the
#'   columns of `x_cond`. When `NULL`, the columns correspond to the last
#'   variables of the current copula order. Discrete left limits are computed
#'   internally from the fitted margins.
#' @details
#' See [vine] for the estimation and construction of vine models.
#' Here, the density, distribution function and random generation
#' for the vine distributions are standard.
#'
#' The functions are based on [dvinecop()], [pvinecop()] and [rvinecop()] for
#' [vinecop] objects. Margins are evaluated through [dmargin()], [pmargin()],
#' and [qmargin()]. Methods are provided for margins fitted by [vine()] and for
#' the fixed [stats::Distributions] specifications accepted by [vine_dist()].
#' @return
#' `dvine()` gives the density, `pvine()` gives the distribution function,
#' and `rvine()` generates unconditional or conditional random deviates.
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
#' vc <- vine_dist(list(stats_margin("norm")), pcs, mat)
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
  cores <- as_count(cores, "cores")
  if (NCOL(x) == 1) {
    x <- t(x)
  }

  x <- expand_factors(x)
  if (!is.null(vine$names)) {
    x <- x[, vine$names, drop = FALSE]
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
  n_mc <- as_count(n_mc, "n_mc")
  cores <- as_count(cores, "cores")

  if (NCOL(x) == 1) {
    x <- t(x)
  }
  x <- expand_factors(x)
  if (!is.null(vine$names)) {
    x <- x[, vine$names, drop = FALSE]
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
rvine <- function(
  n,
  vine,
  qrng = FALSE,
  cores = 1,
  x_cond = NULL,
  conditioning_set = NULL
) {
  n <- as_count(n, "n")
  cores <- as_count(cores, "cores")
  assert_that(
    inherits(vine, "vine_dist"),
    is.flag(qrng)
  )

  if (is.null(x_cond)) {
    if (!is.null(conditioning_set) && length(conditioning_set) > 0) {
      stop("'conditioning_set' requires 'x_cond'.", call. = FALSE)
    }
    U <- rvinecop(n, vine$copula, qrng, cores)
    conditioned_variables <- integer()
  } else {
    x_cond <- process_conditioning_values(x_cond, n, "x_cond")
    d <- dim(vine)[1]
    conditioning_set_supplied <-
      !is.null(conditioning_set) && length(conditioning_set) > 0

    if (conditioning_set_supplied) {
      conditioned_variables <- process_conditioning_set(
        conditioning_set,
        vine$names,
        d
      )
      if (ncol(x_cond) != length(conditioned_variables)) {
        stop(
          "'x_cond' must have one column per conditioning variable.",
          call. = FALSE
        )
      }
    } else {
      k <- ncol(x_cond)
      if (k < 1L || k >= d) {
        stop("'x_cond' must have between 1 and d - 1 columns.", call. = FALSE)
      }
      conditioned_variables <- utils::tail(vine$copula$structure$order, k)
    }

    get_conditioning_column <- function(i) {
      if (is.data.frame(x_cond)) x_cond[[i]] else x_cond[, i]
    }
    u_values <- lapply(seq_along(conditioned_variables), function(i) {
      eval_one_dpq(
        get_conditioning_column(i),
        vine$margins[[conditioned_variables[i]]],
        "p"
      )
    })
    discrete_positions <- which(
      vine$copula$var_types[conditioned_variables] == "d"
    )
    u_left_limits <- lapply(discrete_positions, function(i) {
      eval_one_dpq(
        get_conditioning_column(i),
        vine$margins[[conditioned_variables[i]]],
        "p_sub"
      )
    })
    u_cond <- do.call(cbind, c(u_values, u_left_limits))

    U <- rvinecop(
      n,
      vine$copula,
      qrng,
      cores,
      u_cond = u_cond,
      conditioning_set = if (conditioning_set_supplied) {
        conditioned_variables
      }
    )
  }

  # use quantile transformation for marginals
  X <- dpq_marg(U, vine, "q")
  if (!is.null(vine$var_levels)) {
    for (k in which(lengths(vine$var_levels) > 0L)) {
      values <- if (is.data.frame(X)) X[[k]] else X[, k]
      if (!is.ordered(values)) {
        restored <- ordered(
          vine$var_levels[[k]][as.integer(values) + 1L],
          levels = vine$var_levels[[k]]
        )
        if (is.data.frame(X)) {
          X[[k]] <- restored
        } else {
          X <- as.data.frame(X)
          X[[k]] <- restored
        }
      }
    }
  }
  if (!is.null(x_cond)) {
    for (i in seq_along(conditioned_variables)) {
      if (is.data.frame(X)) {
        X[[conditioned_variables[i]]] <- get_conditioning_column(i)
      } else {
        X[, conditioned_variables[i]] <- get_conditioning_column(i)
      }
    }
  }
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
    copula = summary(object$copula, ...)
  )
}

get_vine_dist_margin_summary <- function(vd) {
  margins <- vd$margins
  if (length(margins) == 1) {
    margins <- rep(list(margins), dim(vd$copula)[1])
  }
  df <- data.frame(
    margin = seq_along(margins),
    distr = vapply(
      margins,
      function(margin) margin_info(margin)$family_name,
      character(1)
    )
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
predict.vine <- function(
  object,
  newdata,
  what = "pdf",
  n_mc = 10^4,
  cores = 1,
  ...
) {
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
  infos <- lapply(object$margins, margin_info)
  support <- t(vapply(infos, `[[`, numeric(2), "support"))
  info <- data.frame(
    margin = seq_along(object$margins),
    name = object$names,
    family = vapply(infos, `[[`, character(1), "family_name"),
    type = vapply(infos, `[[`, character(1), "type"),
    xmin = support[, 1L],
    xmax = support[, 2L],
    npars = vapply(infos, `[[`, numeric(1), "npars"),
    loglik = vapply(infos, `[[`, numeric(1), "loglik")
  )
  class(info) <- c("summary_df", "data.frame")
  info
}


dpq_marg <- function(x, vine, what = "p") {
  d <- ncol(x)
  res <- lapply(
    seq_len(d),
    function(i) eval_one_dpq(x[, i], vine$margins[[i]], what)
  )
  do.call(cbind, res)
}

eval_one_dpq <- function(x, margin, what = "p") {
  if (is.ordered(x) && what != "q") {
    x <- as.numeric(x) - 1L
  }
  dpq <- switch(
    what,
    p = pmargin(x, margin),
    d = dmargin(x, margin),
    q = qmargin(x, margin),
    p_sub = pmargin_sub(x, margin)
  )
  if (is.factor(dpq)) {
    dpq <- as.data.frame(dpq)
  }
  if (what == "p_sub") {
    dpq[is.nan(dpq) & !is.nan(x)] <- 0
  }
  dpq
}

#' @export
dim.vine_dist <- function(x) {
  dim(x$copula)
}
