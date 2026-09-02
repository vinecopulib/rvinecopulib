#' Vine copula models
#'
#' Automated fitting or creation of custom vine copula models
#'
#' @aliases vine_dist
#' @param data a matrix or data.frame. Discrete variables have to be declared as
#' `ordered()`.
#' @param margins_controls a list with arguments to be passed to
#' marginal fitting. The supported entries are
#'   * `family_set` candidate families used for every variable, or a list with
#'   one entry per variable. The default `"kde1d"` preserves nonparametric
#'   margin fitting. Other character names refer to `univariateML` families;
#'   configured built-in and custom candidates are created with
#'   [kde1d_family()], [univariateML_family()], and [margin_family()].
#'   * `selcrit` selection criterion, one of `"loglik"`, `"aic"`, or `"bic"`.
#'   * `cores` number of cores used for marginal fitting. `NULL` inherits the
#'   top-level `cores` argument. Use one core for stochastic custom fitters when
#'   results must be invariant to the number of cores.
#' @param copula_controls a list with arguments to be passed to [vinecop()].
#' @param weights optional numeric vector of weights for each observation,
#'   used for both marginal and copula estimation. Every margin family receives
#'   these weights and is responsible for using or explicitly rejecting them.
#' @param keep_data whether the original data should be stored; if you want to
#'   store the pseudo-observations used for fitting the copula, use the
#'   `copula_controls` argument.
#' @param cores the number of cores to use for parallel computations. Unless
#'   overridden by `margins_controls$cores`, margins are fitted in forked
#'   processes when `cores > 1` on non-Windows systems and serially on Windows.
#'   Stochastic custom margin fits may depend on the number of forked processes.
#' @param var_types optional variable types, one for each variable after factor
#'   expansion: `"c"` for continuous, `"d"` for integer-valued discrete, or
#'   `"zi"` for continuous with an atom at zero. Types are inferred from
#'   [ordered()] and [zero_inflated()] columns when omitted.
#' @details
#' `vine_dist()` creates a vine copula by specifying the margins, a nested list
#' of `bicop_dist` objects and a quadratic structure matrix.
#'
#' `vine()` provides automated fitting for vine copula models. Margin-specific
#' fitting options are properties of the family specification, not controls
#' interpreted by `vine()`. For example, use
#' `kde1d_family(xmin = 0, mult = 2)` in `family_set`. `copula_controls` is a
#' list with the same parameters as [vinecop()] (except for `data`).
#'
#' A character `margins_controls$family_set` is used for every variable. Its
#' entries can be `"kde1d"` or names of models available in univariateML. A
#' list supplies one candidate set per
#' variable, for example `list(c("norm", "cauchy"), c("pois", "geom"))`.
#' Each entry may itself be a character vector, a margin-family object, or a
#' list combining both. An unnamed list is matched by position; a named list
#' must contain every expanded variable name exactly once and is reordered to
#' match the data.
#'
#' The aliases `"par"` and `"parametric"` expand to all `univariateML`
#' families, `"nonpar"` and `"nonparametric"` select `"kde1d"`, and `"all"`
#' combines both. Parametric names and aliases require the suggested
#' `univariateML` package; the default `"kde1d"`, [kde1d_family()], and custom
#' [margin_family()] objects do not. Identical candidate specifications are
#' fitted only once.
#'
#' Candidates incompatible with a variable's type are removed before fitting.
#' rvinecopulib fits every remaining candidate and selects the first minimum of
#' negative log-likelihood, AIC, or BIC according to `selcrit`. Competing
#' candidates therefore need a finite `margin_info()$loglik` value, and AIC or
#' BIC additionally requires a finite `margin_info()$npars`. A sole candidate
#' without a parameter count is retained with a warning, but model AIC and BIC
#' are then unavailable. An unsupported candidate may fail while another
#' candidate succeeds; such failures and all candidate warnings are reported.
#' Protocol violations always stop fitting. The built-in univariateML family
#' explicitly rejects observation weights.
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
#' vc <- vine_dist(list(stats_margin("norm")), pcs, mat)
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
  margins_controls = list(),
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
  cores <- as_count(cores, "cores")
  data <- expand_factors(data)

  d <- ncol(data)

  assert_that(is.list(margins_controls))
  allowed_margins_controls <- c("family_set", "selcrit", "cores")

  # BEGIN legacy margins_controls compatibility ---------
  legacy_args <- c("xmin", "xmax", "type", "mult", "bw", "deg")
  allowed_margins_controls <- c(allowed_margins_controls, legacy_args)
  # END legacy margins_controls compatibility -----------

  assert_that(in_set(names(margins_controls), allowed_margins_controls))

  # BEGIN legacy margins_controls compatibility ---------
  legacy_controls <- upgrade_legacy_margin_controls(
    margins_controls,
    data,
    var_types
  )
  margins_controls <- legacy_controls$margins_controls
  var_types <- legacy_controls$var_types
  # END legacy margins_controls compatibility -----------

  margins_controls <- margins_controls[
    !vapply(margins_controls, is.null, logical(1))
  ]
  margins_controls <- utils::modifyList(
    list(family_set = "kde1d", selcrit = "aic", cores = NULL),
    margins_controls
  )
  var_types <- resolve_margin_types(data, var_types)
  validate_vine_weights(weights, nrow(data))
  marg_cores <- margins_controls$cores
  # not ifelse(): it is vectorised over the condition and would silently
  # truncate a longer `cores` to its first element before as_count() sees it
  if (is.null(marg_cores)) {
    marg_cores <- cores
  }
  marg_cores <- as_count(marg_cores, "margins_controls$cores")

  assert_that(is.list(copula_controls))
  if (is.null(copula_controls$keep_data)) {
    copula_controls$keep_data <- FALSE
  }
  copula_controls$cores <- cores

  ## expand the required arguments and compute default mult if needed
  family_set <- margins_controls$family_set
  family_set <- expand_margin_family_set(
    family_set,
    d,
    colnames(data)
  )

  # BEGIN legacy margins_controls compatibility -----------
  if (!is.null(legacy_controls$kde)) {
    family_set <- apply_legacy_kde_controls(
      family_set,
      legacy_controls$kde,
      data
    )
  }
  # END legacy margins_controls compatibility ---------------

  selcrit <- margins_controls$selcrit
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
  ## estimation of the marginals
  vine <- list()
  margin_data <- prep_for_margins(data)
  vine$margins <- fit_vine_margins(
    margin_data,
    family_set,
    var_types,
    weights,
    selcrit,
    marg_cores
  )
  vine$margins_controls <- list(
    family_set = family_set,
    selcrit = selcrit,
    cores = marg_cores
  )

  ## estimation of the R-vine copula --------------
  vine$copula <- list(var_types = simplify_var_types(var_types))
  copula_controls$var_types <- vine$copula$var_types
  copula_controls$data <- compute_pseudo_obs(data, vine)
  copula_controls$weights <- weights
  vine$copula <- do.call(vinecop, copula_controls)
  vine$copula_controls <- copula_controls[
    -which(names(copula_controls) == "data")
  ]

  finalize_vine(vine, data, weights, keep_data)
}

fit_vine_margins <- function(
  margin_data,
  family_set,
  var_types,
  weights,
  selcrit,
  cores
) {
  cores <- min(cores, length(margin_data))
  if (.Platform$OS.type == "windows") {
    cores <- 1L
  }
  results <- parallel::mclapply(
    seq_along(margin_data),
    fit_one_vine_margin,
    margin_data = margin_data,
    family_set = family_set,
    var_types = var_types,
    weights = weights,
    selcrit = selcrit,
    mc.cores = cores,
    mc.set.seed = TRUE
  )
  collect_vine_margin_results(results)
}

collect_vine_margin_results <- function(results) {
  failed <- which(vapply(
    results,
    function(result) {
      is.null(result) ||
        inherits(result, c("error", "try-error")) ||
        inherits(result$fit, c("error", "try-error"))
    },
    logical(1)
  ))
  if (length(failed)) {
    failure <- results[[failed[1L]]]
    if (is.null(failure)) {
      stop(
        sprintf(
          paste0(
            "forked margin fitting failed for variable %d; retry with ",
            "'margins_controls = list(cores = 1)'."
          ),
          failed[1L]
        ),
        call. = FALSE
      )
    }
    if (inherits(failure, c("error", "try-error"))) {
      stop(failure)
    }
    stop(failure$fit)
  }
  warning_messages <- unique(unlist(
    lapply(results, `[[`, "warnings"),
    use.names = FALSE
  ))
  if (length(warning_messages)) {
    warning(paste(warning_messages, collapse = "\n"), call. = FALSE)
  }
  lapply(results, `[[`, "fit")
}

fit_one_vine_margin <- function(
  j,
  margin_data,
  family_set,
  var_types,
  weights,
  selcrit
) {
  warning_messages <- character()
  fit <- tryCatch(
    withCallingHandlers(
      select_margin(
        margin_data[[j]],
        family_set[[j]],
        var_types[j],
        weights,
        selcrit,
        j
      ),
      warning = function(condition) {
        warning_messages <<- c(warning_messages, conditionMessage(condition))
        invokeRestart("muffleWarning")
      }
    ),
    error = identity
  )
  list(fit = fit, warnings = warning_messages)
}

# BEGIN legacy margins_controls compatibility ----------------------------------
legacy_margin_controls_state <- new.env(parent = emptyenv())

upgrade_legacy_margin_controls <- function(margins_controls, data, var_types) {
  legacy_names <- intersect(
    names(margins_controls),
    c("xmin", "xmax", "type", "mult", "bw", "deg")
  )
  if (!length(legacy_names)) {
    return(list(
      margins_controls = margins_controls,
      var_types = var_types,
      kde = NULL
    ))
  }

  if (!isTRUE(legacy_margin_controls_state$warning_issued)) {
    legacy_margin_controls_state$warning_issued <- TRUE
    warning(
      paste0(
        "Supplying 'xmin', 'xmax', 'mult', 'bw', 'deg', or 'type' in ",
        "'margins_controls' is deprecated; configure KDE options with ",
        "'kde1d_family()' in 'margins_controls$family_set' and supply ",
        "'type' through 'var_types'."
      ),
      call. = FALSE
    )
  }

  legacy_type <- margins_controls$type
  if (!is.null(legacy_type)) {
    legacy_type <- expand_margin_types(
      legacy_type,
      ncol(data),
      "margins_controls$type"
    )
    if (
      !is.null(var_types) &&
        !identical(
          expand_margin_types(var_types, ncol(data), "var_types"),
          legacy_type
        )
    ) {
      stop("'var_types' and 'margins_controls$type' disagree.", call. = FALSE)
    }
    var_types <- legacy_type
  }

  smoothing_names <- intersect(
    legacy_names,
    c("xmin", "xmax", "mult", "bw", "deg")
  )
  explicit_kde <- contains_kde1d_family(margins_controls$family_set)
  if (length(smoothing_names) && explicit_kde) {
    stop(
      paste0(
        "legacy KDE controls cannot be combined with an explicit ",
        "'kde1d_family()' specification."
      ),
      call. = FALSE
    )
  }
  kde <- NULL
  if (!explicit_kde) {
    kde <- expand_legacy_kde_controls(
      margins_controls[legacy_names],
      ncol(data)
    )
  }

  list(
    margins_controls = margins_controls[
      !names(margins_controls) %in% legacy_names
    ],
    var_types = var_types,
    kde = kde
  )
}

contains_kde1d_family <- function(family_set) {
  if (is_kde1d_family(family_set)) {
    return(TRUE)
  }
  if (has_margin_family_protocol(family_set) || !is.list(family_set)) {
    return(FALSE)
  }
  any(vapply(family_set, contains_kde1d_family, logical(1)))
}

is_kde1d_family <- function(family) {
  inherits(family, "custom_margin_family") &&
    identical(family$fit, fit_kde1d_margin)
}

expand_legacy_kde_controls <- function(controls, d) {
  controls <- modifyList(
    list(xmin = NaN, xmax = NaN, type = "c", mult = NULL, bw = NA, deg = 2),
    controls
  )
  if (is.null(controls$mult)) {
    controls$mult <- log(1 + d)
  }
  controls <- lapply(controls, rep_len, length.out = d)
  controls[c("xmin", "xmax", "mult", "bw", "deg")]
}

apply_legacy_kde_controls <- function(family_set, controls, data) {
  structure(
    lapply(seq_along(family_set), function(j) {
      lapply(family_set[[j]], function(family) {
        if (!is_kde1d_family(family)) {
          return(family)
        }
        xmin <- controls$xmin[j]
        xmax <- controls$xmax[j]
        if (is.ordered(data[[j]])) {
          xmin <- 0
          xmax <- nlevels(data[[j]]) - 1
        }
        kde1d_family(
          xmin = xmin,
          xmax = xmax,
          mult = controls$mult[j],
          bw = controls$bw[j],
          deg = controls$deg[j]
        )
      })
    }),
    names = names(family_set)
  )
}
# END legacy margins_controls compatibility ------------------------------------

resolve_margin_types <- function(data, var_types) {
  d <- ncol(data)
  inferred <- rep("c", d)
  inferred[vapply(data, is.ordered, logical(1))] <- "d"
  inferred[vapply(data, inherits, logical(1), "zero_inflated")] <- "zi"

  if (!is.null(var_types)) {
    var_types <- expand_margin_types(var_types, d, "var_types")
  }

  resolved <- if (is.null(var_types)) inferred else var_types
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

validate_vine_weights <- function(weights, n) {
  if (
    !is.numeric(weights) ||
      !length(weights) %in% c(0L, n) ||
      anyNA(weights) ||
      any(!is.finite(weights)) ||
      any(weights < 0) ||
      (length(weights) && sum(weights) <= 0)
  ) {
    stop(
      sprintf(
        "'weights' must be empty or contain %d finite non-negative values with positive sum.",
        n
      ),
      call. = FALSE
    )
  }
  invisible(weights)
}

expand_margin_types <- function(type, d, arg) {
  type <- normalize_margin_types(type, arg)
  if (!length(type) %in% c(1L, d)) {
    stop(sprintf("'%s' must have length one or %d.", arg, d), call. = FALSE)
  }
  rep(type, length.out = d)
}

prep_for_margins <- function(data) {
  lapply(data, function(x) {
    if (!is.ordered(x)) {
      return(x)
    }
    structure(as.numeric(x) - 1, margin_levels = levels(x))
  })
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
      frame <- stats::model.frame(~x, na.action = stats::na.pass)
      x <- model.matrix(attr(frame, "terms"), frame)[, -1, drop = FALSE]
      x <- as.data.frame(x)
      x <- lapply(x, function(y) ordered(y, levels = 0:1))
    })
  }
  as.data.frame(data)
}

#' @param margins a list containing one marginal distribution per variable.
#' Each margin can be a [margin_dist()] object, another fitted object
#' implementing the [margin protocol][margin_protocol], a `kde1d` object, or a
#' [stats_margin()] object. Legacy `list(distr = ...)` specifications are
#' converted through [stats_margin()]. If `margins` has length one, it is
#' recycled for every component.
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
  if (!is.list(margins) || !length(margins) %in% c(1, dim(structure)[1])) {
    stop("marg should have length 1 or dim(structure)[1]")
  }
  if (length(margins) == 1) {
    margins <- replicate(dim(structure)[1], margins[[1]], simplify = FALSE)
  }
  margins <- lapply(margins, as_margin)
  info <- lapply(margins, margin_info)
  npars_marg <- sum(vapply(info, `[[`, numeric(1), "npars"))

  # create the vinecop object
  copula <- vinecop_dist(
    pair_copulas,
    structure,
    var_types = simplify_var_types(vapply(info, `[[`, character(1), "type"))
  )

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

finalize_vine <- function(vine, data, weights, keep_data) {
  ## compute npars/loglik
  info <- lapply(vine$margins, margin_info)

  ## add the npars/loglik of the copulas
  vine$npars <- sum(vapply(info, `[[`, numeric(1), "npars")) +
    vine$copula$npars
  vine$loglik <- sum(vapply(info, `[[`, numeric(1), "loglik")) +
    vine$copula$loglik

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
