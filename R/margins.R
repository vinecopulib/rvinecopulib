#' Marginal distribution protocol
#'
#' Evaluate a fitted marginal distribution. These generics dispatch on
#' `margin`, their second argument. This lets another package integrate its
#' fitted distribution class without changing the usual data-first calling
#' convention.
#'
#' A fitted margin class must provide methods for all three generics and a
#' [logLik()] method whose return value has a finite `df` attribute. It must
#' declare its variable type through `attr(margin, "type")`, using `"c"` for a
#' continuous distribution, `"d"` for an integer-supported distribution, or
#' `"zi"` for a continuous distribution with an atom at zero. For compatibility
#' with `univariateML`, a logical `attr(margin, "continuous")` is also
#' recognized. An optional `"family"` or `"model"` attribute supplies the name
#' shown in summaries.
#'
#' Methods must be vectorized and return one value for every element of `x` or
#' `p`. For discrete margins, `dmargin()` returns probability mass and
#' rvinecopulib computes the left limit as `pmargin(x - 1, margin)`. For a
#' `"zi"` margin, `dmargin(0, margin)` is the mass at zero and the left limit at
#' zero is `pmargin(0, margin) - dmargin(0, margin)`. At nonzero values its left
#' limit equals its CDF.
#'
#' @param x vector of evaluation points.
#' @param p vector of probabilities.
#' @param margin a fitted marginal distribution object.
#'
#' @return A numeric vector containing density or probability mass values for
#'   `dmargin()`, distribution function values for `pmargin()`, and quantiles
#'   for `qmargin()`.
#'
#' @name margin_protocol
#'
#' @examples
#' fitted_normal <- structure(
#'   c(mean = 1, sd = 2),
#'   class = "example_normal_margin",
#'   type = "c",
#'   family = "normal"
#' )
#' dmargin.example_normal_margin <- function(x, margin) {
#'   dnorm(x, margin[["mean"]], margin[["sd"]])
#' }
#' pmargin.example_normal_margin <- function(x, margin) {
#'   pnorm(x, margin[["mean"]], margin[["sd"]])
#' }
#' qmargin.example_normal_margin <- function(p, margin) {
#'   qnorm(p, margin[["mean"]], margin[["sd"]])
#' }
#' logLik.example_normal_margin <- function(object, ...) {
#'   structure(-10, class = "logLik", df = 2)
#' }
#' pmargin(1, fitted_normal)
NULL

#' Create a custom marginal distribution
#'
#' `margin_dist()` creates a fitted or fixed margin from three minimal
#' evaluation functions. The functions must each accept one vector argument
#' and return a vector of the same length. They are stored in the resulting
#' object, so closures can capture fitted parameter values.
#'
#' A finite `loglik` is needed when the object competes with other candidates
#' during [vine()] fitting. It may be `NA` for a fixed margin supplied directly
#' to [vine_dist()] or when it is the only fitted candidate. `npars` is used in
#' AIC and BIC calculations and may be non-integer for effective degrees of
#' freedom.
#'
#' @param d density or probability mass function.
#' @param p distribution function.
#' @param q quantile function.
#' @param family descriptive family name.
#' @param type margin type: `"c"` for continuous, `"d"` for integer-valued
#'   discrete, or `"zi"` for a continuous distribution with an atom at zero.
#' @param npars effective number of fitted parameters.
#' @param loglik maximized marginal log-likelihood, or `NA` for a fixed margin.
#'
#' @return An object of class `margin_dist` implementing the
#'   [margin protocol][margin_protocol].
#' @export
#'
#' @examples
#' normal_margin <- margin_dist(
#'   d = function(x) dnorm(x, mean = 1, sd = 2),
#'   p = function(x) pnorm(x, mean = 1, sd = 2),
#'   q = function(p) qnorm(p, mean = 1, sd = 2),
#'   family = "norm",
#'   npars = 2
#' )
#' dmargin(1, normal_margin)
#' pmargin(1, normal_margin)
#' qmargin(0.5, normal_margin)
margin_dist <- function(
  d,
  p,
  q,
  family = "custom",
  type = "c",
  npars = 0,
  loglik = NA_real_
) {
  if (!is.function(d) || !is.function(p) || !is.function(q)) {
    stop("'d', 'p', and 'q' must be functions.", call. = FALSE)
  }
  if (
    !is.character(family) ||
      length(family) != 1L ||
      is.na(family) ||
      !nzchar(family)
  ) {
    stop("'family' must be a non-empty character string.", call. = FALSE)
  }
  if (
    !is.character(type) || length(type) != 1L || !type %in% c("c", "d", "zi")
  ) {
    stop("'type' must be one of 'c', 'd', or 'zi'.", call. = FALSE)
  }
  if (
    !is.numeric(npars) || length(npars) != 1L || !is.finite(npars) || npars < 0
  ) {
    stop("'npars' must be a finite non-negative number.", call. = FALSE)
  }
  if (
    !is.numeric(loglik) ||
      length(loglik) != 1L ||
      !(is.finite(loglik) || (is.na(loglik) && !is.nan(loglik)))
  ) {
    stop("'loglik' must be a finite number or NA.", call. = FALSE)
  }
  structure(
    list(
      d = d,
      p = p,
      q = q,
      family = family,
      type = type,
      npars = npars,
      loglik = loglik
    ),
    class = c("margin_dist", "list")
  )
}

#' Define a candidate marginal family
#'
#' `margin_family()` defines how [vine()] fits one candidate family. The fitting
#' function has a deliberately small interface: it receives the observed
#' vector and, when supplied to [vine()], a named `weights` argument. It returns
#' an object implementing the
#' [margin protocol][margin_protocol]. Additional settings can be captured in
#' the fitting function's environment. For backward compatibility, `fit` is
#' called with only the observed vector when no observation weights are used.
#'
#' The fitted object must declare the same type as the data being fitted. If
#' several candidates are supplied, it must also have a finite log-likelihood
#' and parameter count. Candidate fitting errors are ignored while another
#' candidate succeeds; if all candidates fail, [vine()] reports their error
#' messages together. The `"kde1d"` candidate supports observation weights.
#' Custom fitters may support them by accepting `weights`; otherwise that
#' candidate fails normally and another candidate may still be selected. The
#' built-in univariateML adapter does not support weights.
#'
#' A `margin_family` object can be used for every variable directly. In a
#' per-variable `family_set`, wrap multiple candidates in a nested list; see
#' the example below and [vine()] for the full family-set syntax.
#'
#' @param fit a function whose first argument is the observed data vector. When
#'   observation weights are supplied, it is also called with a named `weights`
#'   argument containing the corresponding numeric vector.
#' @param family a descriptive family name.
#' @param type variable types supported by the family: any of `"c"`, `"d"`,
#'   and `"zi"`.
#'
#' @return An object of class `margin_family`, suitable for
#'   `margins_controls$family_set` in [vine()].
#' @export
#'
#' @examples
#' normal_family <- margin_family(
#'   fit = function(x) {
#'     mu <- mean(x)
#'     sigma <- sd(x)
#'     margin_dist(
#'       d = function(y) dnorm(y, mu, sigma),
#'       p = function(y) pnorm(y, mu, sigma),
#'       q = function(p) qnorm(p, mu, sigma),
#'       family = "normal",
#'       npars = 2,
#'       loglik = sum(dnorm(x, mu, sigma, log = TRUE))
#'     )
#'   },
#'   family = "normal"
#' )
#' x <- cbind(rnorm(30), rnorm(30, 2))
#' fit <- vine(
#'   x,
#'   margins_controls = list(
#'     family_set = list(list(normal_family), list(normal_family))
#'   ),
#'   copula_controls = list(family_set = "indep")
#' )
#' summary(fit)$margins
margin_family <- function(fit, family = "custom", type = "c") {
  if (!is.function(fit)) {
    stop("'fit' must be a function.", call. = FALSE)
  }
  if (
    !is.character(family) ||
      length(family) != 1L ||
      is.na(family) ||
      !nzchar(family)
  ) {
    stop("'family' must be a non-empty character string.", call. = FALSE)
  }
  type <- unique(normalize_margin_types(type))
  if (!length(type)) {
    stop("'type' must contain at least one variable type.", call. = FALSE)
  }
  structure(
    list(fit = fit, family = family, type = type),
    class = "margin_family"
  )
}

#' Declare zero-inflated data
#'
#' Marks a numeric vector as having a continuous distribution with an atom at
#' zero. The atom is always assumed to be located at zero. The class is
#' preserved when the vector is stored in or subset from a data frame and is
#' detected automatically by [vine()]. Alternatively, use `var_types = "zi"`.
#'
#' For a fitted zero-inflated margin, rvinecopulib uses the ordinary CDF away
#' from zero. At zero, its left limit is computed as
#' `pmargin(0, margin) - dmargin(0, margin)`, so `dmargin(0, margin)` must return
#' the atom's probability mass rather than a continuous density value.
#'
#' @param x a numeric vector.
#'
#' @return `x` with an additional `zero_inflated` class.
#' @export
#'
#' @examples
#' x <- zero_inflated(c(0, 0, rexp(8)))
#' inherits(data.frame(x = x)$x, "zero_inflated")
zero_inflated <- function(x) {
  if (!is.numeric(x)) {
    stop("'x' must be numeric.", call. = FALSE)
  }
  structure(x, class = unique(c("zero_inflated", class(x))))
}

#' @export
`[.zero_inflated` <- function(x, ...) {
  out <- NextMethod("[")
  class(out) <- class(x)
  out
}

#' @rdname margin_protocol
#' @export
dmargin <- function(x, margin) {
  UseMethod("dmargin", margin)
}

#' @rdname margin_protocol
#' @export
pmargin <- function(x, margin) {
  UseMethod("pmargin", margin)
}

#' @rdname margin_protocol
#' @export
qmargin <- function(p, margin) {
  UseMethod("qmargin", margin)
}

#' @export
dmargin.default <- function(x, margin) {
  stop_margin_method("dmargin", margin)
}

#' @export
pmargin.default <- function(x, margin) {
  stop_margin_method("pmargin", margin)
}

#' @export
qmargin.default <- function(p, margin) {
  stop_margin_method("qmargin", margin)
}

stop_margin_method <- function(generic, margin) {
  stop(
    sprintf(
      "no %s() method is available for margin class <%s>",
      generic,
      paste(class(margin), collapse = "/")
    ),
    call. = FALSE
  )
}

#' @export
dmargin.kde1d <- function(x, margin) {
  dkde1d(x, margin)
}

#' @export
pmargin.kde1d <- function(x, margin) {
  pkde1d(x, margin)
}

#' @export
qmargin.kde1d <- function(p, margin) {
  qkde1d(p, margin)
}

#' @export
dmargin.univariateML <- function(x, margin) {
  check_univariateML()
  univariateML::dml(x, margin)
}

#' @export
pmargin.univariateML <- function(x, margin) {
  check_univariateML()
  univariateML::pml(x, margin)
}

#' @export
qmargin.univariateML <- function(p, margin) {
  check_univariateML()
  univariateML::qml(p, margin)
}

#' @export
dmargin.margin_dist <- function(x, margin) {
  margin$d(x)
}

#' @export
pmargin.margin_dist <- function(x, margin) {
  margin$p(x)
}

#' @export
qmargin.margin_dist <- function(p, margin) {
  margin$q(p)
}

#' @export
logLik.margin_dist <- function(object, ...) {
  structure(object$loglik, df = object$npars, class = "logLik")
}

#' @export
dmargin.list <- function(x, margin) {
  eval_legacy_margin("d", x, "x", margin)
}

#' @export
pmargin.list <- function(x, margin) {
  eval_legacy_margin("p", x, "q", margin)
}

#' @export
qmargin.list <- function(p, margin) {
  eval_legacy_margin("q", p, "p", margin)
}

eval_legacy_margin <- function(prefix, values, value_name, margin) {
  if (!is_legacy_margin(margin)) {
    stop_margin_method(paste0(prefix, "margin"), margin)
  }
  args <- margin[names(margin) != "distr"]
  args[[length(args) + 1L]] <- values
  names(args)[length(args)] <- value_name
  do.call(get(paste0(prefix, margin$distr)), args)
}

is_legacy_margin <- function(margin) {
  is.list(margin) &&
    is.character(margin$distr) &&
    length(margin$distr) == 1L &&
    margin$distr %in% supported_distributions
}

has_margin_protocol <- function(margin) {
  classes <- class(margin)
  if (identical(classes, "list")) {
    return(is_legacy_margin(margin))
  }
  generics <- c("dmargin", "pmargin", "qmargin")
  all(vapply(
    generics,
    function(generic) {
      any(vapply(
        classes,
        function(class) {
          !is.null(utils::getS3method(generic, class, optional = TRUE)) &&
            class != "default"
        },
        logical(1)
      ))
    },
    logical(1)
  ))
}

margin_npars <- function(margin) {
  if (inherits(margin, "kde1d")) {
    return(margin$edf)
  }
  if (is_legacy_margin(margin)) {
    return(get_npars_distr(margin))
  }
  ll <- stats::logLik(margin)
  npars <- attr(ll, "df")
  if (is.null(npars) || length(npars) != 1L || !is.finite(npars)) {
    stop("logLik() for a fitted margin must have a finite 'df' attribute.")
  }
  unname(npars)
}

margin_loglik <- function(margin) {
  if (inherits(margin, "kde1d")) {
    return(margin$loglik)
  }
  if (is_legacy_margin(margin)) {
    return(NA_real_)
  }
  as.numeric(stats::logLik(margin))
}

margin_type <- function(margin) {
  type <- NULL
  if (inherits(margin, "margin_dist")) {
    type <- margin$type
  } else if (inherits(margin, "kde1d")) {
    type <- margin$type
  } else if (is_legacy_margin(margin)) {
    type <- "c"
  } else {
    type <- attr(margin, "type", exact = TRUE)
    if (is.null(type)) {
      continuous <- attr(margin, "continuous", exact = TRUE)
      if (isTRUE(continuous)) {
        type <- "c"
      }
      if (isFALSE(continuous)) type <- "d"
    }
  }
  if (is.null(type) || length(type) != 1L) {
    stop("a fitted margin must declare its type.", call. = FALSE)
  }
  normalize_margin_types(type, "margin type")
}

normalize_margin_types <- function(type, arg = "type") {
  if (!is.character(type) || anyNA(type)) {
    stop(sprintf("'%s' must be a character vector.", arg), call. = FALSE)
  }
  normalized <- type
  normalized[type %in% c("cont", "continuous")] <- "c"
  normalized[type %in% c("disc", "discrete")] <- "d"
  normalized[type %in% c("zinf", "zero-inflated")] <- "zi"
  if (any(!normalized %in% c("c", "d", "zi"))) {
    stop(
      sprintf("'%s' must contain only 'c', 'd', or 'zi'.", arg),
      call. = FALSE
    )
  }
  normalized
}

pmargin_sub <- function(x, margin) {
  type <- margin_type(margin)
  if (type == "c") {
    return(pmargin(x, margin))
  }
  if (is.ordered(x) && !inherits(margin, "kde1d")) {
    x <- as.numeric(x) - 1L
  }
  if (type == "d") {
    if (inherits(margin, "kde1d") && is.ordered(margin$x)) {
      xnum <- as.numeric(x)
      levels <- levels(margin$x)
      previous <- ordered(
        levels[ifelse(xnum > 1, xnum - 1, NA)],
        levels = levels
      )
      return(pmargin(previous, margin))
    }
    return(pmargin(x - 1, margin))
  }

  out <- pmargin(x, margin)
  at_zero <- !is.na(x) & x == 0
  out[at_zero] <- pmax(
    out[at_zero] - dmargin(x[at_zero], margin),
    0
  )
  out
}

margin_family_name <- function(margin) {
  if (is_legacy_margin(margin)) {
    return(margin$distr)
  }
  if (inherits(margin, "margin_dist")) {
    return(margin$family)
  }
  family <- attr(margin, "family", exact = TRUE)
  if (is.null(family)) {
    family <- attr(margin, "model", exact = TRUE)
  }
  if (is.null(family)) {
    family <- class(margin)[1L]
  }
  as.character(family)[1L]
}

check_univariateML <- function() {
  if (!requireNamespace("univariateML", quietly = TRUE)) {
    stop(
      "package 'univariateML' is required for parametric margin families.",
      call. = FALSE
    )
  }
}

expand_margin_family_set <- function(family_set, d, variable_names = NULL) {
  if (!is.list(family_set) || inherits(family_set, "margin_family")) {
    family_set <- rep(list(family_set), d)
  } else {
    if (length(family_set) != d) {
      stop(
        sprintf(
          "a list 'family_set' must have one entry for each of %d variables.",
          d
        ),
        call. = FALSE
      )
    }
    if (!is.null(names(family_set))) {
      if (
        is.null(variable_names) ||
          any(!nzchar(names(family_set))) ||
          anyDuplicated(names(family_set)) ||
          !setequal(names(family_set), variable_names)
      ) {
        stop(
          "a named 'family_set' must contain each variable name exactly once.",
          call. = FALSE
        )
      }
      family_set <- family_set[variable_names]
    }
  }
  lapply(family_set, normalize_margin_candidates)
}

normalize_margin_candidates <- function(candidates) {
  if (inherits(candidates, "margin_family")) {
    candidates <- list(candidates)
  } else if (is.character(candidates)) {
    if (!length(candidates) || anyNA(candidates)) {
      stop("margin family sets must be non-empty and cannot contain NA.")
    }
    candidates <- as.list(candidates)
  } else if (is.list(candidates)) {
    candidates <- unlist(
      lapply(candidates, function(candidate) {
        if (is.character(candidate)) as.list(candidate) else list(candidate)
      }),
      recursive = FALSE
    )
  } else {
    stop(
      "margin family candidates must be names or margin_family() objects.",
      call. = FALSE
    )
  }

  if (!length(candidates)) {
    stop("margin family sets must be non-empty.", call. = FALSE)
  }
  valid <- vapply(
    candidates,
    function(candidate) {
      inherits(candidate, "margin_family") ||
        (is.character(candidate) &&
          length(candidate) == 1L &&
          !is.na(candidate))
    },
    logical(1)
  )
  if (!all(valid)) {
    stop(
      "margin family candidates must be names or margin_family() objects.",
      call. = FALSE
    )
  }
  expand_margin_family_aliases(candidates)
}

expand_margin_family_aliases <- function(candidates) {
  aliases <- vapply(
    candidates,
    function(candidate) {
      if (is.character(candidate)) candidate else ""
    },
    character(1)
  )
  aliases[aliases %in% c("nonpar", "nonparametric")] <- "kde1d"
  candidates[aliases == "kde1d"] <- "kde1d"

  expand <- aliases %in% c("all", "par", "parametric")
  if (any(expand)) {
    check_univariateML()
    replacements <- lapply(aliases[expand], function(alias) {
      families <- as.list(univariateML::univariateML_models)
      if (alias == "all") c(list("kde1d"), families) else families
    })
    candidates[which(expand)] <- replacements
    candidates <- unlist(candidates, recursive = FALSE)
  }
  candidates[
    !duplicated(vapply(candidates, margin_candidate_name, character(1)))
  ]
}

margin_candidate_name <- function(candidate) {
  if (inherits(candidate, "margin_family")) candidate$family else candidate
}

margin_candidate_supports <- function(candidate, type) {
  if (inherits(candidate, "margin_family")) {
    return(type %in% candidate$type)
  }
  if (candidate == "kde1d") {
    return(TRUE)
  }
  if (type == "zi") {
    return(FALSE)
  }
  check_univariateML()
  metadata <- univariateML::univariateML_metadata[[paste0("ml", candidate)]]
  if (is.null(metadata)) {
    stop(sprintf("unknown marginal family: \"%s\".", candidate), call. = FALSE)
  }
  support_type <- attr(metadata$support, "type")
  identical(support_type, if (type == "d") "Z" else "R")
}

fit_margin_candidate <- function(x, candidate, type, controls, weights) {
  if (inherits(candidate, "margin_family")) {
    if (length(weights)) {
      fit <- candidate$fit(x, weights = weights)
    } else {
      fit <- candidate$fit(x)
    }
  } else if (candidate == "kde1d") {
    fit <- kde1d::kde1d(
      x,
      xmin = controls$xmin,
      xmax = controls$xmax,
      type = type,
      mult = controls$mult,
      bw = controls$bw,
      deg = controls$deg,
      weights = weights
    )
  } else {
    check_univariateML()
    if (length(weights)) {
      stop(
        "univariateML margin families do not support observation weights.",
        call. = FALSE
      )
    }
    fitter <- getExportedValue("univariateML", paste0("ml", candidate))
    fit <- fitter(x)
    attr(fit, "family") <- candidate
  }

  check <- check_distr(fit)
  if (!isTRUE(check)) {
    stop(check, call. = FALSE)
  }
  fitted_type <- margin_type(fit)
  if (fitted_type != type) {
    stop(
      sprintf(
        "fitted family \"%s\" declared type \"%s\", expected \"%s\".",
        margin_candidate_name(candidate),
        fitted_type,
        type
      ),
      call. = FALSE
    )
  }
  fit
}

select_margin <- function(
  x,
  candidates,
  type,
  controls,
  weights,
  selcrit,
  variable
) {
  compatible <- vapply(
    candidates,
    margin_candidate_supports,
    logical(1),
    type = type
  )
  candidates <- candidates[compatible]
  if (!length(candidates)) {
    stop(
      sprintf(
        "variable %d has no candidate margin supporting type \"%s\".",
        variable,
        type
      ),
      call. = FALSE
    )
  }
  errors <- character()
  fits <- lapply(candidates, function(candidate) {
    tryCatch(
      suppressWarnings(
        fit_margin_candidate(x, candidate, type, controls, weights)
      ),
      error = function(error) {
        errors[[margin_candidate_name(candidate)]] <<- conditionMessage(error)
        NULL
      }
    )
  })
  ok <- !vapply(fits, is.null, logical(1))
  fits <- fits[ok]
  candidates <- candidates[ok]
  if (!length(fits)) {
    details <- paste(names(errors), errors, sep = ": ", collapse = "; ")
    stop(
      sprintf(
        "could not fit a margin for variable %d (%s).",
        variable,
        details
      ),
      call. = FALSE
    )
  }
  if (length(fits) == 1L) {
    return(fits[[1L]])
  }

  loglik <- vapply(fits, margin_loglik, numeric(1))
  npars <- vapply(fits, margin_npars, numeric(1))
  valid <- is.finite(loglik) & is.finite(npars)
  if (!any(valid)) {
    stop(
      sprintf(
        "no candidate margin for variable %d has a finite log-likelihood.",
        variable
      ),
      call. = FALSE
    )
  }
  fits <- fits[valid]
  loglik <- loglik[valid]
  npars <- npars[valid]
  criterion <- switch(
    selcrit,
    loglik = -loglik,
    aic = -2 * loglik + 2 * npars,
    bic = -2 * loglik + log(length(x)) * npars
  )
  fits[[which.min(criterion)]]
}
