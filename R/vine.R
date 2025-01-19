#' Vine copula models
#'
#' Automated fitting or creation of custom vine copula models
#'
#' @aliases vine_dist
#' @param data a matrix or data.frame. Discrete variables have to be declared as
#' `ordered()`.
#' @param margins_controls a list with arguments to be passed to
#' [kde1d::kde1d()]. Currently, there can be
#'   * `mult` numeric vector of length one or d; all bandwidths for marginal
#'   kernel density estimation are multiplied with `mult`. Defaults to
#'   `log(1 + d)` where `d` is the number of variables after applying
#'   `rvinecopulib:::expand_factors()`.
#'   * `xmin` numeric vector of length d; see [kde1d::kde1d()].
#'   * `xmax` numeric vector of length d; see [kde1d::kde1d()].
#'   * `bw` numeric vector of length d; see [kde1d::kde1d()].
#'   * `deg` numeric vector of length one or d; [kde1d::kde1d()].
#' @param copula_controls a list with arguments to be passed to [vinecop()].
#' @param weights optional vector of weights for each observation.
#' @param keep_data whether the original data should be stored; if you want to
#'   store the pseudo-observations used for fitting the copula, use the
#'   `copula_controls` argument.
#' @param cores the number of cores to use for parallel computations.
#' @details
#' `vine_dist()` creates a vine copula by specifying the margins, a nested list
#' of `bicop_dist` objects and a quadratic structure matrix.
#'
#' `vine()` provides automated fitting for vine copula models.
#' `margins_controls` is a list with the same parameters as
#' [kde1d::kde1d()] (except for `x`). `copula_controls` is a list
#' with the same parameters as [vinecop()] (except for `data`).
#'
#' @return Objects inheriting from `vine_dist` for [vine_dist()], and
#' `vine` and `vine_dist` for [vine()].
#'
#' Objects from the `vine_dist` class are lists containing:
#'
#' * `margins`, a list of marginals (see below).
#' * `copula`, an object of the class `vinecop_dist`, see [vinecop_dist()].
#'
#' For objects from the `vine` class, `copula` is also an object of the class
#' `vine`, see [vinecop()]. Additionally, objects from the `vine` class contain:
#'
#' * `margins_controls`, a `list` with the set of fit controls that was passed
#' to [kde1d::kde1d()] when estimating the margins.
#' * `copula_controls`, a `list` with the set of fit controls that was passed
#' to [vinecop()] when estimating the copula.
#' * `data` (optionally, if `keep_data = TRUE` was used), the dataset that was
#' passed to [vine()].
#' * `nobs`, an `integer` containing the number of observations that was used
#' to fit the model.
#'
#' Concerning `margins`:
#'
#' * For objects created with [vine_dist()], it simply corresponds to the `margins`
#' argument.
#' * For objects created with [vine()], it is a list of objects of class `kde1d`,
#' see [kde1d::kde1d()].
#'
#' @examples
#' # specify pair-copulas
#' bicop <- bicop_dist("bb1", 90, c(3, 2))
#' pcs <- list(
#'   list(bicop, bicop), # pair-copulas in first tree
#'   list(bicop) # pair-copulas in second tree
#' )
#'
#' # specify R-vine matrix
#' mat <- matrix(c(1, 2, 3, 1, 2, 0, 1, 0, 0), 3, 3)
#'
#' # set up vine copula model with Gaussian margins
#' vc <- vine_dist(list(distr = "norm"), pcs, mat)
#'
#' # show model
#' summary(vc)
#'
#' # simulate some data
#' x <- rvine(50, vc)
#'
#' # estimate a vine copula model
#' fit <- vine(x, copula_controls = list(family_set = "par"))
#' summary(fit)
#'
#' ## model for discrete data
#' x <- as.data.frame(x)
#' x[, 1] <- ordered(round(x[, 1]), levels = seq.int(-5, 5))
#' fit_disc <- vine(x, copula_controls = list(family_set = "par"))
#' summary(fit_disc)
#'
#' @importFrom kde1d dkde1d pkde1d qkde1d
#' @export
vine <- function(data,
                 margins_controls = list(
                   mult = NULL,
                   xmin = NaN,
                   xmax = NaN,
                   bw = NA,
                   deg = 2
                 ),
                 copula_controls = list(
                   family_set = "all",
                   structure = NA,
                   par_method = "mle",
                   nonpar_method = "constant",
                   mult = 1,
                   selcrit = "aic",
                   psi0 = 0.9,
                   presel = TRUE,
                   allow_rotations = TRUE,
                   trunc_lvl = Inf,
                   tree_crit = "tau",
                   threshold = 0,
                   keep_data = FALSE,
                   show_trace = FALSE,
                   cores = 1
                 ),
                 weights = numeric(),
                 keep_data = FALSE,
                 cores = 1) {
  ## basic sanity checks (copula_controls are checked by vinecop)
  data <- expand_factors(data)

  d <- ncol(data)
  var_types <- rep("c", d)
  var_types[sapply(data, is.ordered)] <- "d"

  assert_that(is.list(margins_controls))
  allowed_margins_controls <- c("xmin", "xmax", "type", "mult", "bw", "deg")
  assert_that(in_set(names(margins_controls), allowed_margins_controls))

  assert_that(is.list(copula_controls))
  if (is.null(copula_controls$keep_data)) {
    copula_controls$keep_data <- TRUE
  }
  copula_controls$cores <- cores

  ## expand the required arguments and compute default mult if needed
  margins_controls <- expand_margin_controls(margins_controls, d, data)
  for (k in which(sapply(data, is.ordered))) {
    margins_controls$type[k] <- "d"
    margins_controls$xmin[k] <- 0
    margins_controls$xmax[k] <- nlevels(data[[k]]) - 1
  }

  ## estimation of the marginals
  vine <- list()
  vine$margins <- fit_margins_cpp(prep_for_margins(data),
                                  xmin = margins_controls$xmin,
                                  xmax = margins_controls$xmax,
                                  type = margins_controls$type,
                                  mult = margins_controls$mult,
                                  bw = margins_controls$bw,
                                  deg = margins_controls$deg,
                                  weights = weights,
                                  cores)
  vine$margins_controls <- margins_controls
  vine$margins <- finalize_margins(data, vine$margins)

  ## estimation of the R-vine copula --------------
  copula_controls$data <- compute_pseudo_obs(data, vine)
  copula_controls$var_types <- simplify_var_types(margins_controls$type)
  copula_controls$weights <- weights
  vine$copula <- do.call(vinecop, copula_controls)
  vine$copula_controls <- copula_controls[-which(names(copula_controls) == "data")]

  finalize_vine(vine, data, weights, keep_data)
}

prep_for_margins <- function(data) {
  data <- lapply(data, function(x) if (is.ordered(x)) as.numeric(x) - 1 else x)
  do.call(cbind, data)
}

compute_pseudo_obs <- function(data, vine) {
  d <- ncol(data)
  u <- dpq_marg(data, vine)
  if (any(sapply(data, is.factor))) {
    u_sub <- u
    for (k in seq_len(d)) {
      if (is.factor(data[, k])) {
        lv <- as.numeric(data[, k]) - 1
        lv0 <- which(lv == 0)
        lv[lv0] <- 1
        xlv <- ordered(levels(data[, k])[lv], levels = levels(data[, k]))
        u_sub[, k] <- eval_one_dpq(xlv, vine$margins[[k]])
        u_sub[lv0, k] <- 0
      }
    }
  } else {
    u_sub <- NULL
  }
  cbind(u, u_sub)
}

#' @importFrom stats model.matrix
#' @noRd
expand_factors <- function(data) {
  if (is.data.frame(data)) {
    data <- lapply(data, function(x) {
      if (is.numeric(x) | is.ordered(x))
        return(x)
      x <- model.matrix(~ x)[, -1, drop = FALSE]
      x <- as.data.frame(x)
      x <- lapply(x, function(y) ordered(y, levels = 0:1))
    })
  }
  as.data.frame(data)
}

#' @param margins A list with with each element containing the specification of a
#' marginal [stats::Distributions]. Each marginal specification
#' should be a list with containing at least the distribution family (`"distr"`)
#' and optionally the parameters, e.g.
#' `list(list(distr = "norm"), list(distr = "norm", mu = 1), list(distr = "beta", shape1 = 1, shape2 = 1))`.
#' Note that parameters that have no default values have to be provided.
#' Furthermore, if `margins` has length one, it will be recycled for every component.
#' @param pair_copulas A nested list of 'bicop_dist' objects, where
#'    \code{pair_copulas[[t]][[e]]} corresponds to the pair-copula at edge `e` in
#'    tree `t`.
#' @param structure an `rvine_structure` object, namely a compressed
#' representation of the vine structure, or an object that can be coerced
#' into one (see [rvine_structure()] and [as_rvine_structure()]).
#' The dimension must be `length(pair_copulas[[1]]) + 1`.
#' @rdname vine
#' @export
vine_dist <- function(margins, pair_copulas, structure) {
  structure <- as_rvine_structure(structure)

  # sanity checks for the marg
  if (!(length(margins) %in% c(1, dim(structure)[1]))) {
    stop("marg should have length 1 or dim(structure)[1]")
  }
  stopifnot(is.list(margins))
  if (depth(margins) == 1) {
    check_marg <- check_distr(margins)
    try(npars_marg <- ncol(matrix) * get_npars_distr(margins), silent = TRUE)
  } else {
    check_marg <- lapply(margins, check_distr)
    try(npars_marg <- sum(sapply(margins, get_npars_distr)), silent = TRUE)
  }
  is_ok <- sapply(check_marg, isTRUE)
  if (!all(is_ok)) {
    msg <- "Some objects in marg aren't properly defined.\n"
    msg <- c(msg, paste0("margin ", seq_along(check_marg)[!is_ok], " : ",
                         unlist(check_marg[!is_ok]), ".",
                         sep = "\n"
    ))
    stop(msg)
  }

  # create the vinecop object
  copula <- vinecop_dist(pair_copulas, structure)

  # create object
  structure(list(
    margins = margins,
    copula = copula,
    npars = copula$npars + npars_marg,
    loglik = NA
  ), class = "vine_dist")
}

expand_margin_controls <- function(controls, d, data) {
  default_controls <-
    list(xmin = NaN, xmax = NaN, type = "c", mult = NULL, bw = NA, deg = 2)
  controls <- modifyList(default_controls, controls)
  if (is.null(controls[["mult"]])) {
    controls[["mult"]] <- log(1 + d)
  }
  for (par in names(controls)) {
    if (length(controls[[par]]) != ncol(data))
      controls[[par]] <- rep(controls[[par]], ncol(data))
  }
  controls
}

finalize_margins <- function(data, margins) {
  for (k in seq_along(margins)) {
    margins[[k]]$x <- data[[k]]
    margins[[k]]$nobs <- nrow(data)
  }
  margins
}

finalize_vine <- function(vine, data, weights, keep_data) {
  ## compute npars/loglik
  npars <- loglik <- 0
  for (k in seq_len(ncol(data))) {
    npars <- npars + vine$margins[[k]]$edf
    loglik <- loglik + vine$margins[[k]]$loglik
  }

  ## add the npars/loglik of the copulas
  vine$npars <- npars + vine$copula$npars
  vine$loglik <- loglik + vine$copula$loglik

  ## add data
  if (keep_data) {
    vine$data <- data
    vine$weights <- weights
  } else {
    vine$data <- matrix(NA, ncol = ncol(data))
    colnames(vine$data) <- colnames(data)
  }

  ## add number of observations
  vine$nobs <- nrow(data)
  vine$names <- vine$copula$names <- colnames(data)

  ## store levels for discrete variables
  vine$var_levels <- lapply(data, levels)

  ## create and return object
  structure(vine, class = c("vine", "vine_dist"))
}

simplify_var_types <- function(x) {
  x[x %in% c("cont", "continuous")] <- "c"
  x[x %in% c("disc", "discrete", "zinf", "zero-inflated")] <- "d"
  x
}
