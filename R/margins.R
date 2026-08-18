#' Fitted marginal distribution protocol
#'
#' Fitted margins provide distribution evaluation and model metadata through a
#' small set of S3 generics. Methods for `dmargin()`, `pmargin()`, and
#' `qmargin()` dispatch on `margin`, their second argument. Metadata methods
#' dispatch on their only argument.
#'
#' A fitted margin class must implement all eight generics documented here.
#' Evaluation methods must be vectorized. `margin_type()` returns `"c"` for a
#' continuous distribution, `"d"` for an integer-supported distribution, or
#' `"zi"` for a continuous distribution with an atom at zero.
#' `margin_support()` returns the mathematical support as its lower and upper
#' bounds. The parameter count must be non-negative or `NA` when unavailable;
#' the log-likelihood may be `NA` for a fixed margin.
#'
#' The type determines left limits throughout rvinecopulib. They are `F(x)` for
#' continuous margins, `F(x - 1)` for integer-supported margins, and
#' `F(0) - P(X = 0)` at zero for zero-inflated margins. Consequently,
#' `dmargin(0, margin)` must return the atom's mass for a zero-inflated margin.
#'
#' @param x vector of evaluation points.
#' @param p vector of probabilities.
#' @param margin a fitted marginal distribution object.
#'
#' @return `dmargin()`, `pmargin()`, and `qmargin()` return numeric vectors.
#'   `margin_type()` and `margin_family_name()` return character scalars,
#'   `margin_support()` returns the two support bounds, and `margin_npars()` and
#'   `margin_loglik()` return numeric scalars.
#'
#' @name margin_protocol
NULL

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

#' @rdname margin_protocol
#' @export
margin_type <- function(margin) {
  UseMethod("margin_type")
}

#' @rdname margin_protocol
#' @export
margin_support <- function(margin) {
  UseMethod("margin_support")
}

#' @rdname margin_protocol
#' @export
margin_family_name <- function(margin) {
  UseMethod("margin_family_name")
}

#' @rdname margin_protocol
#' @export
margin_npars <- function(margin) {
  UseMethod("margin_npars")
}

#' @rdname margin_protocol
#' @export
margin_loglik <- function(margin) {
  UseMethod("margin_loglik")
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

#' @export
margin_type.default <- function(margin) {
  stop_margin_method("margin_type", margin)
}

#' @export
margin_support.default <- function(margin) {
  stop_margin_method("margin_support", margin)
}

#' @export
margin_family_name.default <- function(margin) {
  stop_margin_method("margin_family_name", margin)
}

#' @export
margin_npars.default <- function(margin) {
  stop_margin_method("margin_npars", margin)
}

#' @export
margin_loglik.default <- function(margin) {
  stop_margin_method("margin_loglik", margin)
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

#' Create a fitted or fixed custom margin
#'
#' `margin_dist()` creates a fitted or fixed margin from evaluation functions
#' and explicit metadata. The functions must each accept one vector argument
#' and return a numeric vector of the same length.
#'
#' @param d density or probability mass function.
#' @param p distribution function.
#' @param q quantile function.
#' @param family descriptive family name.
#' @param type margin type: `"c"`, `"d"`, or `"zi"`.
#' @param support mathematical support as its lower and upper bounds.
#' @param npars effective number of fitted parameters.
#' @param loglik maximized marginal log-likelihood, or `NA` for a fixed margin.
#'
#' @return A fitted margin implementing the [margin protocol][margin_protocol].
#' @export
#'
#' @examples
#' normal_margin <- margin_dist(
#'   d = function(x) dnorm(x, mean = 1, sd = 2),
#'   p = function(x) pnorm(x, mean = 1, sd = 2),
#'   q = function(p) qnorm(p, mean = 1, sd = 2),
#'   family = "norm",
#'   support = c(-Inf, Inf),
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
  support = c(-Inf, Inf),
  npars = 0,
  loglik = NA_real_
) {
  if (!is.function(d) || !is.function(p) || !is.function(q)) {
    stop("'d', 'p', and 'q' must be functions.", call. = FALSE)
  }
  validate_margin_family_name(family)
  type <- normalize_margin_types(type, "type")
  if (length(type) != 1L) {
    stop("'type' must have length one.", call. = FALSE)
  }
  validate_margin_support(support)
  validate_margin_npars(npars)
  validate_margin_loglik(loglik)
  margin <- structure(
    list(
      d = d,
      p = p,
      q = q,
      family = family,
      type = type,
      support = unname(support),
      npars = unname(npars),
      loglik = unname(loglik)
    ),
    class = c("margin_dist", "list")
  )
  validate_margin(margin)
  margin
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
margin_type.margin_dist <- function(margin) {
  margin$type
}

#' @export
margin_support.margin_dist <- function(margin) {
  margin$support
}

#' @export
margin_family_name.margin_dist <- function(margin) {
  margin$family
}

#' @export
margin_npars.margin_dist <- function(margin) {
  margin$npars
}

#' @export
margin_loglik.margin_dist <- function(margin) {
  margin$loglik
}

#' @export
logLik.margin_dist <- function(object, ...) {
  structure(object$loglik, df = object$npars, class = "logLik")
}

#' Create a fixed margin from a stats distribution
#'
#' `stats_margin()` adapts the univariate distributions documented in
#' [stats::Distributions] to the fitted-margin protocol. Distribution
#' parameters are supplied through `...`; the resulting margin is fixed and
#' therefore has zero fitted parameters and no fitted log-likelihood.
#'
#' @param family distribution suffix such as `"norm"`, `"pois"`, or `"binom"`.
#' @param ... parameters passed to the distribution's density, CDF, and quantile
#'   functions.
#'
#' @return A fixed margin implementing the [margin protocol][margin_protocol].
#' @export
#'
#' @examples
#' normal <- stats_margin("norm", mean = 1, sd = 2)
#' poisson <- stats_margin("pois", lambda = 3)
#' margin_type(normal)
#' margin_type(poisson)
stats_margin <- function(family, ...) {
  validate_margin_family_name(family)
  if (!family %in% names(stats_margin_types)) {
    stop(
      sprintf("unsupported stats distribution: \"%s\".", family),
      call. = FALSE
    )
  }
  margin <- structure(
    list(
      family = family,
      args = list(...),
      type = stats_margin_types[[family]]
    ),
    class = c("stats_margin", "list")
  )
  tryCatch(
    qmargin(c(0, 0.5, 1), margin),
    error = function(error) {
      stop(conditionMessage(error), call. = FALSE)
    }
  )
  margin$support <- unname(qmargin(c(0, 1), margin))
  validate_margin_support(margin$support)
  validate_margin(margin)
  margin
}

stats_margin_types <- c(
  beta = "c",
  cauchy = "c",
  chisq = "c",
  exp = "c",
  f = "c",
  gamma = "c",
  logis = "c",
  lnorm = "c",
  norm = "c",
  t = "c",
  unif = "c",
  weibull = "c",
  binom = "d",
  geom = "d",
  hyper = "d",
  nbinom = "d",
  pois = "d",
  signrank = "d",
  wilcox = "d"
)

#' @export
dmargin.stats_margin <- function(x, margin) {
  eval_stats_margin("d", x, "x", margin)
}

#' @export
pmargin.stats_margin <- function(x, margin) {
  eval_stats_margin("p", x, "q", margin)
}

#' @export
qmargin.stats_margin <- function(p, margin) {
  eval_stats_margin("q", p, "p", margin)
}

eval_stats_margin <- function(prefix, values, value_name, margin) {
  args <- margin$args
  args[[value_name]] <- values
  do.call(getExportedValue("stats", paste0(prefix, margin$family)), args)
}

#' @export
margin_type.stats_margin <- function(margin) {
  margin$type
}

#' @export
margin_support.stats_margin <- function(margin) {
  margin$support
}

#' @export
margin_family_name.stats_margin <- function(margin) {
  margin$family
}

#' @export
margin_npars.stats_margin <- function(margin) {
  0
}

#' @export
margin_loglik.stats_margin <- function(margin) {
  NA_real_
}

#' @export
dmargin.kde1d <- function(x, margin) {
  kde1d::dkde1d(x, margin)
}

#' @export
pmargin.kde1d <- function(x, margin) {
  kde1d::pkde1d(x, margin)
}

#' @export
qmargin.kde1d <- function(p, margin) {
  kde1d::qkde1d(p, margin)
}

#' @export
margin_type.kde1d <- function(margin) {
  normalize_margin_types(margin$type, "margin type")
}

#' @export
margin_support.kde1d <- function(margin) {
  c(
    if (is.nan(margin$xmin)) -Inf else margin$xmin,
    if (is.nan(margin$xmax)) Inf else margin$xmax
  )
}

#' @export
margin_family_name.kde1d <- function(margin) {
  "kde1d"
}

#' @export
margin_npars.kde1d <- function(margin) {
  if (is.nan(margin$edf)) NA_real_ else margin$edf
}

#' @export
margin_loglik.kde1d <- function(margin) {
  margin$loglik
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
margin_type.univariateML <- function(margin) {
  type <- attr(margin, "rvine_margin_type", exact = TRUE)
  if (is.null(type)) {
    continuous <- attr(margin, "continuous", exact = TRUE)
    type <- if (isTRUE(continuous)) {
      "c"
    } else if (isFALSE(continuous)) {
      "d"
    }
  }
  normalize_margin_types(type, "margin type")
}

#' @export
margin_support.univariateML <- function(margin) {
  support <- attr(margin, "rvine_margin_support", exact = TRUE)
  if (is.null(support)) {
    support <- qmargin(c(0, 1), margin)
  }
  unname(support)
}

#' @export
margin_family_name.univariateML <- function(margin) {
  family <- attr(margin, "rvine_margin_family", exact = TRUE)
  if (is.null(family)) {
    family <- attr(margin, "model", exact = TRUE)
  }
  if (is.null(family)) class(margin)[1L] else as.character(family)[1L]
}

#' @export
margin_npars.univariateML <- function(margin) {
  npars <- attr(stats::logLik(margin), "df")
  validate_margin_npars(npars)
  unname(npars)
}

#' @export
margin_loglik.univariateML <- function(margin) {
  as.numeric(stats::logLik(margin))
}

#' Normalize an object to the fitted-margin protocol
#'
#' `as_margin()` validates fitted margins and adapts the legacy
#' `list(distr = ...)` representation through [stats_margin()].
#'
#' @param margin a fitted margin or legacy stats distribution list.
#'
#' @return A validated fitted margin.
#' @export
as_margin <- function(margin) {
  UseMethod("as_margin")
}

#' @export
as_margin.default <- function(margin) {
  validate_margin(margin)
  margin
}

#' @export
as_margin.list <- function(margin) {
  if (has_margin_protocol(margin)) {
    validate_margin(margin)
    return(margin)
  }
  if (
    is.null(margin$distr) ||
      !is.character(margin$distr) ||
      length(margin$distr) != 1L
  ) {
    stop(
      "a legacy stats margin must be a list with one 'distr' value.",
      call. = FALSE
    )
  }
  do.call(
    stats_margin,
    c(list(family = margin$distr), margin[names(margin) != "distr"])
  )
}

has_margin_protocol <- function(margin) {
  all(vapply(
    c(
      "dmargin",
      "pmargin",
      "qmargin",
      "margin_type",
      "margin_support",
      "margin_family_name",
      "margin_npars",
      "margin_loglik"
    ),
    has_s3_method,
    logical(1),
    object = margin
  ))
}

has_s3_method <- function(generic, object) {
  any(vapply(
    class(object),
    function(class) {
      !is.null(utils::getS3method(generic, class, optional = TRUE))
    },
    logical(1)
  ))
}

validate_margin <- function(margin, x = NULL) {
  if (!has_margin_protocol(margin)) {
    stop(
      sprintf(
        "object of class <%s> does not implement the fitted-margin protocol.",
        paste(class(margin), collapse = "/")
      ),
      call. = FALSE
    )
  }
  type <- margin_type(margin)
  normalized <- normalize_margin_types(type, "margin_type() result")
  if (length(type) != 1L || !identical(type, normalized)) {
    stop(
      "margin_type() must return one canonical variable type.",
      call. = FALSE
    )
  }
  validate_margin_support(margin_support(margin))
  validate_margin_family_name(margin_family_name(margin))
  validate_margin_npars(margin_npars(margin))
  validate_margin_loglik(margin_loglik(margin))
  validate_margin_evaluation(margin, x)
  invisible(margin)
}

check_distr <- function(distr) {
  tryCatch(
    {
      as_margin(distr)
      TRUE
    },
    error = function(error) conditionMessage(error)
  )
}


stop_wo_call <- function(message) {
  stop(message, call. = FALSE)
}

validate_margin_evaluation <- function(margin, x) {
  quantiles <- qmargin(c(0.25, 0.75), margin)
  if (
    !is.numeric(quantiles) ||
      length(quantiles) != 2L ||
      anyNA(quantiles)
  ) {
    stop_wo_call("qmargin() must return one numeric value per input.")
  }
  if (is.null(x)) {
    probes <- quantiles
  } else {
    probes <- x[seq_len(min(3L, length(x)))]
    if (anyNA(x) && !anyNA(probes)) {
      probes[length(probes)] <- NA
    }
  }
  if (length(probes)) {
    density <- dmargin(probes, margin)
    proba <- pmargin(probes, margin)
    if (!is.numeric(density) || length(density) != length(probes)) {
      stop_wo_call("dmargin() must return one numeric value per input.")
    }
    if (!all(is.na(density) == is.na(probes))) {
      stop_wo_call("dmargin() must propagate missing values.")
    }
    if (any(density[!is.na(density)] < 0)) {
      stop_wo_call("dmargin() returned a negative value.")
    }
    if (!is.numeric(proba) || length(proba) != length(probes)) {
      stop_wo_call("pmargin() must return one numeric value per input.")
    }
    if (!all(is.na(proba) == is.na(probes))) {
      stop_wo_call("pmargin() must propagate missing values.")
    }
    if (any(proba[!is.na(proba)] < 0 | proba[!is.na(proba)] > 1)) {
      stop_wo_call("pmargin() returned a value outside [0, 1].")
    }
  }
  invisible(margin)
}


validate_margin_family_name <- function(name) {
  if (
    !is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)
  ) {
    stop_wo_call("a margin family name must be a non-empty character scalar.")
  }
  invisible(name)
}

validate_margin_support <- function(support) {
  if (
    !is.numeric(support) ||
      length(support) != 2L ||
      anyNA(support) ||
      support[1L] > support[2L]
  ) {
    stop_wo_call(
      "margin support must contain ordered numeric lower and upper bounds."
    )
  }
  invisible(support)
}

validate_margin_npars <- function(npars) {
  if (
    !is.numeric(npars) ||
      length(npars) != 1L ||
      !(is.finite(npars) && npars >= 0 || is.na(npars) && !is.nan(npars))
  ) {
    stop_wo_call("'npars' must be finite and non-negative or NA.")
  }
  invisible(npars)
}

validate_margin_loglik <- function(loglik) {
  if (
    !is.numeric(loglik) ||
      length(loglik) != 1L ||
      !(is.finite(loglik) || (is.na(loglik) && !is.nan(loglik)))
  ) {
    stop_wo_call("'loglik' must be finite or NA.")
  }
  invisible(loglik)
}

normalize_margin_types <- function(type, arg = "type") {
  if (!is.character(type) || anyNA(type)) {
    stop_wo_call(sprintf("'%s' must be a character vector.", arg))
  }
  normalized <- type
  normalized[type %in% c("cont", "continuous")] <- "c"
  normalized[type %in% c("disc", "discrete")] <- "d"
  normalized[type %in% c("zinf", "zero-inflated")] <- "zi"
  if (any(!normalized %in% c("c", "d", "zi"))) {
    stop_wo_call(sprintf("'%s' must contain only 'c', 'd', or 'zi'.", arg))
  }
  normalized
}

pmargin_sub <- function(x, margin) {
  type <- margin_type(margin)
  if (type == "c") {
    return(pmargin(x, margin))
  }
  if (is.ordered(x)) {
    x <- as.numeric(x) - 1L
  }
  if (type == "d") {
    return(pmargin(x - 1, margin))
  }
  probability <- pmargin(x, margin)
  at_zero <- !is.na(x) & x == 0
  probability[at_zero] <- pmax(
    probability[at_zero] - dmargin(x[at_zero], margin),
    0
  )
  probability
}

#' Margin-family fitting protocol
#'
#' Margin families describe how a candidate is fitted. Implementations declare
#' a name and their supported variable types and provide a `fit_margin()`
#' method. Every fitter receives the observations, observation weights, and the
#' requested variable type. It must either use the weights or reject them.
#'
#' @param family a margin-family specification.
#' @param x observed values.
#' @param weights observation weights, or `numeric()`.
#' @param type requested variable type.
#'
#' @return `fit_margin()` returns a fitted margin. `margin_family_types()`
#'   returns the supported variable types.
#' @name margin_family_protocol
NULL

#' @rdname margin_family_protocol
#' @export
fit_margin <- function(family, x, weights = numeric(), type = "c") {
  UseMethod("fit_margin")
}

#' @rdname margin_family_protocol
#' @export
margin_family_types <- function(family) {
  UseMethod("margin_family_types")
}

#' @export
fit_margin.default <- function(family, x, weights = numeric(), type = "c") {
  stop(
    sprintf(
      "object of class <%s> does not implement the margin-family protocol.",
      paste(class(family), collapse = "/")
    ),
    call. = FALSE
  )
}

#' @export
margin_family_types.default <- function(family) {
  stop(
    sprintf(
      "no margin_family_types() method is available for class <%s>.",
      paste(class(family), collapse = "/")
    ),
    call. = FALSE
  )
}

#' Define a custom margin family
#'
#' `margin_family()` defines a candidate family for [vine()]. Its fitting
#' function must accept `x`, `weights`, and `type`, and return an object
#' implementing the [fitted-margin protocol][margin_protocol]. Additional
#' settings can be captured in the fitting function's environment.
#'
#' @param fit function accepting `x`, `weights`, and `type`.
#' @param family_name descriptive family name.
#' @param types supported variable types.
#'
#' @return A margin-family specification.
#' @export
#'
#' @examples
#' normal_family <- margin_family(
#'   fit = function(x, weights, type) {
#'     if (length(weights)) {
#'       location <- weighted.mean(x, weights)
#'       scale <- sqrt(weighted.mean((x - location)^2, weights))
#'     } else {
#'       location <- mean(x)
#'       scale <- sd(x)
#'     }
#'     margin_dist(
#'       d = function(y) dnorm(y, location, scale),
#'       p = function(y) pnorm(y, location, scale),
#'       q = function(p) qnorm(p, location, scale),
#'       family = "normal",
#'       type = type,
#'       npars = 2,
#'       loglik = sum(dnorm(x, location, scale, log = TRUE))
#'     )
#'   },
#'   family_name = "normal"
#' )
margin_family <- function(fit, family_name = "custom", types = "c") {
  if (!is.function(fit)) {
    stop("'fit' must be a function.", call. = FALSE)
  }
  arguments <- names(formals(fit))
  if (
    !length(arguments) ||
      !"x" %in% arguments ||
      (!all(c("weights", "type") %in% arguments) && !"..." %in% arguments)
  ) {
    stop("'fit' must accept 'x', 'weights', and 'type'.", call. = FALSE)
  }
  validate_margin_family_name(family_name)
  types <- unique(normalize_margin_types(types, "types"))
  if (!length(types)) {
    stop("'types' must contain at least one variable type.", call. = FALSE)
  }
  structure(
    list(fit = fit, family_name = family_name, types = types),
    class = c("custom_margin_family", "margin_family")
  )
}

#' @export
fit_margin.custom_margin_family <- function(
  family,
  x,
  weights = numeric(),
  type = "c"
) {
  family$fit(x = x, weights = weights, type = type)
}

#' @export
margin_family_types.custom_margin_family <- function(family) {
  family$types
}

#' @export
margin_family_name.custom_margin_family <- function(margin) {
  margin$family_name
}

#' Define a kde1d margin family
#'
#' The arguments configure [kde1d::kde1d()] and travel with the family object,
#' rather than being interpreted by [vine()].
#'
#' @param xmin,xmax,mult,bw,deg arguments passed to [kde1d::kde1d()].
#'
#' @return A margin-family specification supporting all rvinecopulib variable
#'   types.
#' @export
kde1d_family <- function(xmin = NaN, xmax = NaN, mult = 1, bw = NA, deg = 2) {
  if (
    !is.numeric(xmin) ||
      length(xmin) != 1L ||
      (is.na(xmin) && !is.nan(xmin))
  ) {
    stop("'xmin' must be one number or NaN.", call. = FALSE)
  }
  if (
    !is.numeric(xmax) ||
      length(xmax) != 1L ||
      (is.na(xmax) && !is.nan(xmax))
  ) {
    stop("'xmax' must be one number or NaN.", call. = FALSE)
  }
  if (!is.nan(xmin) && !is.nan(xmax) && xmin >= xmax) {
    stop("'xmin' must be smaller than 'xmax'.", call. = FALSE)
  }
  if (
    !is.numeric(mult) || length(mult) != 1L || !is.finite(mult) || mult <= 0
  ) {
    stop("'mult' must be one positive finite number.", call. = FALSE)
  }
  if (
    length(bw) != 1L ||
      !(is.na(bw) && !is.nan(bw) || is.numeric(bw) && is.finite(bw) && bw > 0)
  ) {
    stop("'bw' must be one positive finite number or NA.", call. = FALSE)
  }
  if (
    !is.numeric(deg) ||
      length(deg) != 1L ||
      is.na(deg) ||
      !deg %in% 0:2
  ) {
    stop("'deg' must be 0, 1, or 2.", call. = FALSE)
  }
  structure(
    list(xmin = xmin, xmax = xmax, mult = mult, bw = as.numeric(bw), deg = deg),
    class = c("kde1d_margin_family", "margin_family")
  )
}

#' @export
fit_margin.kde1d_margin_family <- function(
  family,
  x,
  weights = numeric(),
  type = "c"
) {
  xmin <- family$xmin
  xmax <- family$xmax
  levels <- attr(x, "margin_levels", exact = TRUE)
  if (length(levels)) {
    if (is.nan(xmin)) {
      xmin <- 0
    }
    if (is.nan(xmax)) xmax <- length(levels) - 1
  }
  kde1d::kde1d(
    x,
    xmin = xmin,
    xmax = xmax,
    type = type,
    mult = family$mult,
    bw = family$bw,
    deg = family$deg,
    weights = weights
  )
}

#' @export
margin_family_types.kde1d_margin_family <- function(family) {
  c("c", "d", "zi")
}

#' @export
margin_family_name.kde1d_margin_family <- function(margin) {
  "kde1d"
}

#' Define a univariateML margin family
#'
#' `univariateML_family()` adapts one of univariateML's available models to the
#' margin-family protocol.
#'
#' @param family a univariateML family name.
#'
#' @return A parametric margin-family specification.
#' @export
univariateML_family <- function(family) {
  check_univariateML()
  validate_margin_family_name(family)
  metadata <- univariateML::univariateML_metadata[[paste0("ml", family)]]
  if (is.null(metadata)) {
    stop(sprintf("unknown univariateML family: \"%s\".", family), call. = FALSE)
  }
  support_type <- attr(metadata$support, "type")
  structure(
    list(
      family = family,
      types = if (identical(support_type, "Z")) "d" else "c"
    ),
    class = c("univariateML_margin_family", "margin_family")
  )
}

#' @export
fit_margin.univariateML_margin_family <- function(
  family,
  x,
  weights = numeric(),
  type = "c"
) {
  check_univariateML()
  if (length(weights)) {
    stop(
      "univariateML margin families do not support observation weights.",
      call. = FALSE
    )
  }
  fit <- getExportedValue("univariateML", paste0("ml", family$family))(x)
  attr(fit, "rvine_margin_family") <- family$family
  attr(fit, "rvine_margin_type") <- type
  attr(fit, "rvine_margin_support") <- unname(qmargin(c(0, 1), fit))
  fit
}

#' @export
margin_family_types.univariateML_margin_family <- function(family) {
  family$types
}

#' @export
margin_family_name.univariateML_margin_family <- function(margin) {
  margin$family
}

check_univariateML <- function() {
  if (!requireNamespace("univariateML", quietly = TRUE)) {
    stop(
      paste(
        "package 'univariateML' is required to fit named parametric",
        "margin families in vine()."
      ),
      call. = FALSE
    )
  }
}

has_margin_family_protocol <- function(family) {
  all(vapply(
    c("fit_margin", "margin_family_types", "margin_family_name"),
    has_s3_method,
    logical(1),
    object = family
  ))
}

validate_margin_family <- function(family) {
  if (!has_margin_family_protocol(family)) {
    stop(
      sprintf(
        "object of class <%s> does not implement the margin-family protocol.",
        paste(class(family), collapse = "/")
      ),
      call. = FALSE
    )
  }
  validate_margin_family_name(margin_family_name(family))
  types <- margin_family_types(family)
  normalized <- normalize_margin_types(types, "margin family types")
  if (!length(normalized) || !identical(types, normalized)) {
    stop(
      "margin_family_types() must return supported canonical variable types.",
      call. = FALSE
    )
  }
  invisible(family)
}

expand_margin_family_set <- function(
  family_set,
  d,
  variable_names = NULL
) {
  if (!is.list(family_set) || has_margin_family_protocol(family_set)) {
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
  if (has_margin_family_protocol(candidates)) {
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
      "margin family candidates must be names or family objects.",
      call. = FALSE
    )
  }
  if (!length(candidates)) {
    stop("margin family sets must be non-empty.", call. = FALSE)
  }
  candidates <- unlist(
    lapply(candidates, expand_margin_family_alias),
    recursive = FALSE
  )
  candidates <- lapply(candidates, as_margin_family)
  lapply(candidates, validate_margin_family)
  candidates[!duplicated(candidates)]
}

expand_margin_family_alias <- function(candidate) {
  if (!is.character(candidate) || length(candidate) != 1L) {
    return(list(candidate))
  }
  if (candidate %in% c("nonpar", "nonparametric")) {
    return(list("kde1d"))
  }
  if (candidate %in% c("all", "par", "parametric")) {
    check_univariateML()
    families <- as.list(univariateML::univariateML_models)
    if (candidate == "all") c(list("kde1d"), families) else families
  } else {
    list(candidate)
  }
}

as_margin_family <- function(family) {
  if (has_margin_family_protocol(family)) {
    return(family)
  }
  if (!is.character(family) || length(family) != 1L || is.na(family)) {
    stop(
      "margin family candidates must be names or family objects.",
      call. = FALSE
    )
  }
  if (family == "kde1d") {
    return(kde1d_family())
  }
  univariateML_family(family)
}

fit_margin_candidate <- function(x, family, type, weights) {
  fit <- fit_margin(family, x, weights = weights, type = type)
  validate_margin(fit, x)
  if (margin_type(fit) != type) {
    stop(
      sprintf(
        "fitted family \"%s\" declared type \"%s\", expected \"%s\".",
        margin_family_name(family),
        margin_type(fit),
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
  weights,
  selcrit,
  variable
) {
  compatible <- vapply(
    candidates,
    function(candidate) type %in% margin_family_types(candidate),
    logical(1)
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
      suppressWarnings(fit_margin_candidate(x, candidate, type, weights)),
      error = function(error) {
        errors <<- c(
          errors,
          sprintf(
            "%s: %s",
            margin_family_name(candidate),
            conditionMessage(error)
          )
        )
        NULL
      }
    )
  })
  fitted <- !vapply(fits, is.null, logical(1))
  fits <- fits[fitted]
  if (!length(fits)) {
    stop(
      sprintf(
        "could not fit a margin for variable %d (%s).",
        variable,
        paste(errors, collapse = "; ")
      ),
      call. = FALSE
    )
  }
  if (length(fits) == 1L) {
    return(fits[[1L]])
  }
  loglik <- vapply(fits, margin_loglik, numeric(1))
  valid <- is.finite(loglik)
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
  npars <- vapply(fits, margin_npars, numeric(1))
  if (selcrit != "loglik") {
    valid <- is.finite(npars)
    if (!any(valid)) {
      stop(
        sprintf(
          "no candidate margin for variable %d has a finite parameter count.",
          variable
        ),
        call. = FALSE
      )
    }
    fits <- fits[valid]
    loglik <- loglik[valid]
    npars <- npars[valid]
  }
  criterion <- switch(
    selcrit,
    loglik = -loglik,
    aic = -2 * loglik + 2 * npars,
    bic = -2 * loglik + log(sum(!is.na(x))) * npars
  )
  fits[[which.min(criterion)]]
}

#' Declare zero-inflated data
#'
#' Marks a numeric vector as continuous with an atom at zero. The class is
#' preserved when stored in or subset from a data frame and is detected by
#' [vine()]. Alternatively, use `var_types = "zi"`.
#'
#' @param x a numeric vector.
#'
#' @return `x` with an additional `zero_inflated` class.
#' @export
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
