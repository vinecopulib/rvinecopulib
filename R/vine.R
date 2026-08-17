#' Vine copula models
#'
#' Automated fitting or creation of custom vine copula models
#'
#' @aliases vine_dist
#' @param data a matrix or data.frame. Discrete variables have to be declared as
#' `ordered()`.
#' @param margins_controls a list with arguments to be passed to
#' marginal fitting. Currently, there can be
#'   * `family_set` candidate families used for every variable, or a list with
#'   one entry per variable. The default `"kde1d"` preserves nonparametric
#'   margin fitting. Other character names refer to `univariateML` families;
#'   custom candidates can be defined with [margin_family()].
#'   * `selcrit` selection criterion, one of `"loglik"`, `"aic"`, or `"bic"`.
#'   * `mult` numeric vector of length one or d; all bandwidths for marginal
#'   kernel density estimation are multiplied with `mult`. Defaults to
#'   `log(1 + d)` where `d` is the number of variables after applying
#'   `rvinecopulib:::expand_factors()`.
#'   * `xmin` numeric vector of length d; see [kde1d::kde1d()].
#'   * `xmax` numeric vector of length d; see [kde1d::kde1d()].
#'   * `type` legacy variable-type control, a character vector of length one or
#'   d. Prefer the top-level `var_types` argument for new code.
#'   * `bw` numeric vector of length d; see [kde1d::kde1d()].
#'   * `deg` numeric vector of length one or d; [kde1d::kde1d()].
#' @param copula_controls a list with arguments to be passed to [vinecop()].
#' @param weights optional vector of weights for each observation.
#' @param keep_data whether the original data should be stored; if you want to
#'   store the pseudo-observations used for fitting the copula, use the
#'   `copula_controls` argument.
#' @param cores the number of cores to use for parallel computations.
#' @param var_types optional variable types, one for each variable after factor
#'   expansion: `"c"` for continuous, `"d"` for integer-valued discrete, or
#'   `"zi"` for continuous with an atom at zero. Types are inferred from
#'   [ordered()] and [zero_inflated()] columns when omitted.
#' @details
#' `vine_dist()` creates a vine copula by specifying the margins, a nested list
#' of `bicop_dist` objects and a quadratic structure matrix.
#'
#' `vine()` provides automated fitting for vine copula models.
#' `margins_controls` controls marginal family selection and contains the same
#' smoothing parameters as [kde1d::kde1d()] (except for `x`). `copula_controls`
#' is a list with the same parameters as [vinecop()] (except for `data`).
#'
#' A character `margins_controls$family_set` is used for every variable. Its
#' entries can be `"kde1d"` or names from
#' [univariateML::univariateML_models]. A list supplies one candidate set per
#' variable, for example `list(c("norm", "cauchy"), c("pois", "geom"))`.
#' Each entry may itself be a character vector, a [margin_family()] object, or a
#' list combining both. An unnamed list is matched by position; a named list
#' must contain every expanded variable name exactly once and is reordered to
#' match the data.
#'
#' The aliases `"par"` and `"parametric"` expand to all `univariateML`
#' families, `"nonpar"` and `"nonparametric"` select `"kde1d"`, and `"all"`
#' combines both. Parametric names and aliases require the suggested
#' `univariateML` package; the default `"kde1d"` and custom [margin_family()]
#' objects do not. Duplicate candidate names are fitted only once.
#'
#' Candidates incompatible with a variable's type are removed before fitting.
#' rvinecopulib fits every remaining candidate and selects the first minimum of
#' `-logLik`, AIC, or BIC according to `selcrit`. Competing custom candidates
#' therefore need a finite [logLik()] value. Observation weights are currently
#' supported only when every candidate is `"kde1d"`.
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
#' * `margins_controls`, a `list` with the controls used to fit and select the
#' marginal families.
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
#' * For objects created with [vine()], it is a list of fitted objects
#' implementing the [margin protocol][margin_protocol].
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
#' vc <- vine_dist(list(list(distr = "norm")), pcs, mat)
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
#' # parametric margin selection (requires the suggested univariateML package)
#' if (requireNamespace("univariateML", quietly = TRUE)) {
#'   fit_par_margins <- vine(
#'     x,
#'     margins_controls = list(family_set = c("norm", "cauchy"))
#'   )
#' }
#'
#' ## model for discrete data
#' x <- as.data.frame(x)
#' x[, 1] <- ordered(round(x[, 1]), levels = seq.int(-5, 5))
#' fit_disc <- vine(x, copula_controls = list(family_set = "par"))
#' summary(fit_disc)
#'
#' @importFrom kde1d dkde1d pkde1d qkde1d
#' @export
vine <- function(
  data,
  margins_controls = list(
    family_set = "kde1d",
    selcrit = "aic",
    mult = NULL,
    xmin = NaN,
    xmax = NaN,
    type = "c",
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
    cores = 1,
    tree_algorithm = "mst_prim",
    conditioning_set = NULL
  ),
  weights = numeric(),
  keep_data = FALSE,
  cores = 1,
  var_types = NULL
) {
  ## basic sanity checks (copula_controls are checked by vinecop)
  margins_controls_supplied <- !missing(margins_controls)
  data <- expand_factors(data)

  d <- ncol(data)

  assert_that(is.list(margins_controls))
  allowed_margins_controls <- c(
    "family_set",
    "selcrit",
    "xmin",
    "xmax",
    "type",
    "mult",
    "bw",
    "deg"
  )
  assert_that(in_set(names(margins_controls), allowed_margins_controls))
  var_types <- resolve_margin_types(
    data,
    var_types,
    margins_controls,
    margins_controls_supplied
  )

  assert_that(is.list(copula_controls))
  if (is.null(copula_controls$keep_data)) {
    copula_controls$keep_data <- TRUE
  }
  copula_controls$cores <- cores

  ## expand the required arguments and compute default mult if needed
  family_set <- margins_controls$family_set
  if (is.null(family_set)) {
    family_set <- "kde1d"
  }
  family_set <- expand_margin_family_set(family_set, d, colnames(data))
  selcrit <- margins_controls$selcrit
  if (is.null(selcrit)) {
    selcrit <- "aic"
  }
  if (!is.character(selcrit) || length(selcrit) != 1L || is.na(selcrit)) {
    stop("'margins_controls$selcrit' must be a single string.", call. = FALSE)
  }
  selcrit <- tolower(selcrit)
  if (!selcrit %in% c("loglik", "aic", "bic")) {
    stop(
      "'margins_controls$selcrit' must be 'loglik', 'aic', or 'bic'.",
      call. = FALSE
    )
  }
  kde_controls <- margins_controls[
    intersect(
      names(margins_controls),
      c("xmin", "xmax", "type", "mult", "bw", "deg")
    )
  ]
  kde_controls <- expand_margin_controls(kde_controls, d, data)
  margins_controls <- kde_controls
  margins_controls$type <- var_types
  for (k in which(sapply(data, is.ordered))) {
    margins_controls$type[k] <- "d"
    margins_controls$xmin[k] <- 0
    margins_controls$xmax[k] <- nlevels(data[[k]]) - 1
  }

  ## estimation of the marginals
  vine <- list()
  margin_data <- prep_for_margins(data)
  only_kde1d <- all(vapply(
    family_set,
    function(candidates) {
      length(candidates) == 1L && identical(candidates[[1L]], "kde1d")
    },
    logical(1)
  ))
  if (only_kde1d) {
    vine$margins <- fit_margins_cpp(
      margin_data,
      xmin = margins_controls$xmin,
      xmax = margins_controls$xmax,
      type = margins_controls$type,
      mult = margins_controls$mult,
      bw = margins_controls$bw,
      deg = margins_controls$deg,
      weights = weights,
      cores
    )
  } else {
    vine$margins <- lapply(seq_len(d), function(j) {
      controls <- lapply(margins_controls, `[[`, j)
      select_margin(
        margin_data[, j],
        family_set[[j]],
        var_types[j],
        controls,
        weights,
        selcrit,
        j
      )
    })
  }
  margins_controls$family_set <- family_set
  margins_controls$selcrit <- selcrit
  vine$margins_controls <- margins_controls
  vine$margins <- finalize_margins(data, vine$margins)

  ## estimation of the R-vine copula --------------
  vine$copula <- list(var_types = simplify_var_types(margins_controls$type))
  copula_controls$var_types <- vine$copula$var_types
  copula_controls$data <- compute_pseudo_obs(data, vine)
  copula_controls$weights <- weights
  vine$copula <- do.call(vinecop, copula_controls)
  vine$copula_controls <- copula_controls[
    -which(names(copula_controls) == "data")
  ]

  finalize_vine(vine, data, weights, keep_data)
}

resolve_margin_types <- function(
  data,
  var_types,
  margins_controls,
  margins_controls_supplied
) {
  d <- ncol(data)
  inferred <- rep("c", d)
  inferred[vapply(data, is.ordered, logical(1))] <- "d"
  inferred[vapply(data, inherits, logical(1), "zero_inflated")] <- "zi"

  controls_type <- if (margins_controls_supplied) {
    margins_controls$type
  } else {
    NULL
  }
  if (!is.null(controls_type)) {
    controls_type <- expand_margin_types(
      controls_type,
      d,
      "margins_controls$type"
    )
  }
  if (!is.null(var_types)) {
    var_types <- expand_margin_types(var_types, d, "var_types")
  }
  if (
    !is.null(var_types) &&
      !is.null(controls_type) &&
      !identical(var_types, controls_type)
  ) {
    stop("'var_types' and 'margins_controls$type' disagree.", call. = FALSE)
  }

  explicit <- if (!is.null(var_types)) var_types else controls_type
  resolved <- if (is.null(explicit)) inferred else explicit
  marked <- inferred != "c"
  if (any(marked & resolved != inferred)) {
    stop(
      "explicit variable types disagree with ordered or zero_inflated columns.",
      call. = FALSE
    )
  }

  for (j in which(resolved == "d")) {
    if (!is.ordered(data[[j]])) {
      x <- data[[j]]
      if (any(!is.na(x) & (!is.finite(x) | x != floor(x)))) {
        stop(
          sprintf("discrete variable %d must be integer-valued.", j),
          call. = FALSE
        )
      }
    }
  }
  resolved
}

expand_margin_types <- function(type, d, arg) {
  type <- normalize_margin_types(type, arg)
  if (!length(type) %in% c(1L, d)) {
    stop(sprintf("'%s' must have length one or %d.", arg, d), call. = FALSE)
  }
  rep(type, length.out = d)
}

prep_for_margins <- function(data) {
  data <- lapply(data, function(x) if (is.ordered(x)) as.numeric(x) - 1 else x)
  do.call(cbind, data)
}

compute_pseudo_obs <- function(data, vine) {
  d <- ncol(data)
  u <- dpq_marg(data, vine)
  colnames(u) <- names(data)
  if (any(vine$copula$var_types != "c")) {
    u_sub <- u
    for (k in which(vine$copula$var_types != "c")) {
      u_sub[, k] <- eval_one_dpq(data[[k]], vine$margins[[k]], "p_sub")
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
      if (is.numeric(x) | is.ordered(x)) {
        return(x)
      }
      x <- model.matrix(~x)[, -1, drop = FALSE]
      x <- as.data.frame(x)
      x <- lapply(x, function(y) ordered(y, levels = 0:1))
    })
  }
  as.data.frame(data)
}

#' @param margins a list containing one marginal distribution per variable.
#' Each margin can be a [margin_dist()] object, another fitted object
#' implementing the [margin protocol][margin_protocol], a `kde1d` object, or a
#' fixed [stats::Distributions] specification. Fixed specifications must be a
#' list containing at least the distribution family (`"distr"`) and optionally
#' the parameters, e.g.
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
  if (length(margins) == 1) {
    margins <- replicate(dim(structure)[1], margins[[1]], simplify = FALSE)
  }
  stopifnot(length(margins) == dim(structure)[1])
  check_marg <- lapply(margins, check_distr)

  is_ok <- sapply(check_marg, isTRUE)
  if (!all(is_ok)) {
    msg <- "Some objects in marg aren't properly defined.\n"
    msg <- c(
      msg,
      paste0(
        "margin ",
        seq_along(check_marg)[!is_ok],
        " : ",
        unlist(check_marg[!is_ok]),
        ".",
        sep = "\n"
      )
    )
    stop(msg)
  }
  npars_marg <- sum(vapply(margins, margin_npars, numeric(1)))

  # create the vinecop object
  copula <- vinecop_dist(pair_copulas, structure)

  # create object
  structure(
    list(
      margins = margins,
      copula = copula,
      npars = copula$npars + npars_marg,
      loglik = NA
    ),
    class = "vine_dist"
  )
}

expand_margin_controls <- function(controls, d, data) {
  default_controls <-
    list(xmin = NaN, xmax = NaN, type = "c", mult = NULL, bw = NA, deg = 2)
  controls <- modifyList(default_controls, controls)
  if (is.null(controls[["mult"]])) {
    controls[["mult"]] <- log(1 + d)
  }
  for (par in names(controls)) {
    if (length(controls[[par]]) != ncol(data)) {
      controls[[par]] <- rep(controls[[par]], ncol(data))
    }
  }
  controls
}

finalize_margins <- function(data, margins) {
  for (k in seq_along(margins)) {
    if (inherits(margins[[k]], "kde1d")) {
      margins[[k]]$x <- data[[k]]
      margins[[k]]$nobs <- nrow(data)
    } else if (is.ordered(data[[k]])) {
      attr(margins[[k]], "levels") <- levels(data[[k]])
    }
  }
  margins
}

finalize_vine <- function(vine, data, weights, keep_data) {
  ## compute npars/loglik
  npars <- loglik <- 0
  for (k in seq_len(ncol(data))) {
    npars <- npars + margin_npars(vine$margins[[k]])
    loglik <- loglik + margin_loglik(vine$margins[[k]])
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
  x[x %in% c("disc", "discrete", "zi", "zinf", "zero-inflated")] <- "d"
  x
}
