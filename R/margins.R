#' Fitted marginal distribution protocol
#'
#' Fitted margins provide distribution evaluation and model metadata through a
#' small set of S3 generics. Methods for `dmargin()`, `pmargin()`, and
#' `qmargin()` dispatch on `margin`, their second argument. `margin_info()`
#' dispatches on its only argument.
#'
#' A fitted margin class must implement all four generics documented here.
#' Evaluation methods must be vectorized. `margin_info()` returns a list with
#' entries `family_name`, `type`, `support`, `npars`, and `loglik`. `type` is
#' `"c"` for a continuous distribution, `"d"` for an integer-supported
#' distribution, or `"zi"` for a continuous distribution with an atom at zero.
#' `support` contains the mathematical lower and upper bounds. `npars` must be
#' non-negative or `NA` when unavailable; `loglik` may be `NA` for a fixed
#' margin.
#'
#' The type determines left limits throughout rvinecopulib. They are `F(x)` for
#' continuous margins, `F(x - 1)` for integer-supported margins, and
#' `F(0) - P(X = 0)` at zero for zero-inflated margins. Consequently,
#' `dmargin(0, margin)` must return the atom's mass for a zero-inflated margin.
#'
#' @param x vector of evaluation points.
#' @param p vector of probabilities.
#' @param margin a fitted marginal distribution object.
#' @param object a fitted margin or margin-family specification.
#'
#' @return `dmargin()`, `pmargin()`, and `qmargin()` return numeric vectors.
#'   `margin_info()` returns a list of model information.
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
margin_info <- function(object) {
  UseMethod("margin_info")
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
#' @param npars effective number of parameters.
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
      info = list(
        family_name = family,
        type = type,
        support = unname(support),
        npars = unname(npars),
        loglik = unname(loglik)
      )
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
margin_info.margin_dist <- function(object) {
  object$info
}

#' @export
logLik.margin_dist <- function(object, ...) {
  structure(
    object$info$loglik,
    df = object$info$npars,
    class = "logLik"
  )
}

#' Create a fixed margin from a stats distribution
#'
#' `stats_margin()` adapts the univariate distributions documented in
#' [stats::Distributions] to the fitted-margin protocol. Distribution
#' parameters are supplied through `...`; the resulting margin reports the
#' distribution's parameter count but has no fitted log-likelihood.
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
#' margin_info(normal)
#' margin_info(poisson)
stats_margin <- function(family, ...) {
  validate_margin_family_name(family)
  if (!family %in% names(stats_margin_types)) {
    stop(
      sprintf("unsupported stats distribution: \"%s\".", family),
      call. = FALSE
    )
  }
  args <- list(...)
  d <- getExportedValue("stats", paste0("d", family))
  p <- getExportedValue("stats", paste0("p", family))
  q <- getExportedValue("stats", paste0("q", family))
  structure(
    margin_dist(
      d = function(x) do.call(d, c(list(x = x), args)),
      p = function(x) do.call(p, c(list(q = x), args)),
      q = function(probabilities) {
        do.call(q, c(list(p = probabilities), args))
      },
      family = family,
      type = stats_margin_types[[family]],
      support = unname(do.call(q, c(list(p = c(0, 1)), args))),
      npars = get_stats_margin_npars(family, args)
    ),
    class = c("stats_margin", "margin_dist", "list")
  )
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

get_stats_margin_npars <- function(family, args) {
  npars <- c(
    beta = 2,
    cauchy = 2,
    chisq = 2,
    exp = 1,
    f = 3,
    gamma = 2,
    logis = 2,
    lnorm = 2,
    norm = 2,
    t = 2,
    unif = 2,
    weibull = 2,
    binom = 2,
    geom = 1,
    hyper = 3,
    nbinom = 2,
    pois = 1,
    signrank = 1,
    wilcox = 2
  )[[family]]
  if (family %in% c("chisq", "t") && !"ncp" %in% names(args)) {
    npars <- npars - 1
  }
  if (family == "weibull" && !"scale" %in% names(args)) {
    npars <- npars - 1
  }
  npars
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
margin_info.kde1d <- function(object) {
  list(
    family_name = "kde1d",
    type = normalize_margin_types(object$type, "margin type"),
    support = c(
      if (is.nan(object$xmin)) -Inf else object$xmin,
      if (is.nan(object$xmax)) Inf else object$xmax
    ),
    npars = if (is.nan(object$edf)) NA_real_ else object$edf,
    loglik = object$loglik
  )
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
margin_info.univariateML <- function(object) {
  type <- attr(object, "rvine_margin_type", exact = TRUE)
  if (is.null(type)) {
    continuous <- attr(object, "continuous", exact = TRUE)
    type <- if (isTRUE(continuous)) {
      "c"
    } else if (isFALSE(continuous)) {
      "d"
    }
  }
  support <- attr(object, "rvine_margin_support", exact = TRUE)
  if (is.null(support)) {
    support <- qmargin(c(0, 1), object)
  }
  family <- attr(object, "rvine_margin_family", exact = TRUE)
  if (is.null(family)) {
    family <- attr(object, "model", exact = TRUE)
  }
  family <- if (is.null(family)) class(object)[1L] else as.character(family)[1L]
  loglik <- stats::logLik(object)
  npars <- attr(loglik, "df")
  validate_margin_npars(npars)
  list(
    family_name = family,
    type = normalize_margin_types(type, "margin type"),
    support = unname(support),
    npars = unname(npars),
    loglik = as.numeric(loglik)
  )
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
  if (has_margin_protocol(margin)) {
    validate_margin(margin)
    return(margin)
  }
  if (!inherits(margin, "list")) {
    validate_margin(margin)
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
  has_s3_methods(margin, c("dmargin", "pmargin", "qmargin", "margin_info"))
}

has_s3_methods <- function(object, generics) {
  all(vapply(
    generics,
    function(generic) {
      any(vapply(
        class(object),
        function(class) {
          !is.null(utils::getS3method(generic, class, optional = TRUE))
        },
        logical(1)
      ))
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
  info <- margin_info(margin)
  if (
    !is_margin_info(
      info,
      c("family_name", "type", "support", "npars", "loglik")
    )
  ) {
    stop(
      paste0(
        "margin_info() must return a list with 'family_name', 'type', ",
        "'support', 'npars', and 'loglik'."
      ),
      call. = FALSE
    )
  }
  normalized <- normalize_margin_types(info$type, "margin_info()$type")
  if (length(info$type) != 1L || !identical(info$type, normalized)) {
    stop(
      "margin_info()$type must be one canonical variable type.",
      call. = FALSE
    )
  }
  validate_margin_support(info$support)
  validate_margin_family_name(info$family_name)
  validate_margin_npars(info$npars)
  validate_margin_loglik(info$loglik)
  validate_margin_evaluation(margin, x)
  invisible(margin)
}

is_margin_info <- function(info, fields) {
  is.list(info) &&
    !is.null(names(info)) &&
    !any(!nzchar(names(info))) &&
    !anyDuplicated(names(info)) &&
    all(fields %in% names(info))
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
  type <- margin_info(margin)$type
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
#' Margin families describe how a candidate is fitted. Implementations provide
#' `fit_margin()` and `margin_info()` methods. Their `margin_info()` result must
#' contain `family_name` and `types`. Every fitter receives the observations,
#' observation weights, and requested variable type. It must either use the
#' weights or reject them.
#'
#' @param family a margin-family specification.
#' @param x observed values.
#' @param weights observation weights, or `numeric()`.
#' @param type requested variable type.
#'
#' @return `fit_margin()` returns a fitted margin.
#' @name margin_family_protocol
NULL

#' @rdname margin_family_protocol
#' @export
fit_margin <- function(family, x, weights = numeric(), type = "c") {
  UseMethod("fit_margin")
}

#' Define a custom margin family
#'
#' `margin_family()` defines a candidate family for [vine()]. Its fitting
#' function must accept `x`, `weights`, and `type`, and return an object
#' implementing the [fitted-margin protocol][margin_protocol]. Additional
#' named fitting arguments can be stored in `fit_args`.
#'
#' @param fit function accepting `x`, `weights`, and `type`.
#' @param family_name descriptive family name.
#' @param types supported variable types.
#' @param fit_args named list of additional arguments passed to `fit`.
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
margin_family <- function(
  fit,
  family_name = "custom",
  types = "c",
  fit_args = list()
) {
  if (!is.function(fit)) {
    stop("'fit' must be a function.", call. = FALSE)
  }
  arguments <- names(formals(fit))
  if (!all(c("x", "weights", "type") %in% arguments)) {
    stop("'fit' must accept 'x', 'weights', and 'type'.", call. = FALSE)
  }
  if (
    !is.list(fit_args) ||
      (length(fit_args) &&
        (is.null(names(fit_args)) ||
          any(!nzchar(names(fit_args))) ||
          anyDuplicated(names(fit_args)) ||
          any(names(fit_args) %in% c("x", "weights", "type"))))
  ) {
    stop(
      "'fit_args' must be uniquely named and cannot contain 'x', 'weights', or 'type'.",
      call. = FALSE
    )
  }
  if (
    length(fit_args) &&
      !"..." %in% arguments &&
      any(!names(fit_args) %in% arguments)
  ) {
    stop("every 'fit_args' entry must be accepted by 'fit'.", call. = FALSE)
  }
  validate_margin_family_name(family_name)
  types <- unique(normalize_margin_types(types, "types"))
  if (!length(types)) {
    stop("'types' must contain at least one variable type.", call. = FALSE)
  }
  structure(
    list(
      fit = fit,
      fit_args = fit_args,
      info = list(family_name = family_name, types = types)
    ),
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
  do.call(
    family$fit,
    c(list(x = x, weights = weights, type = type), family$fit_args)
  )
}

#' @export
margin_info.custom_margin_family <- function(object) {
  object$info
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
  margin_family(
    fit = fit_kde1d_margin,
    fit_args = list(
      xmin = xmin,
      xmax = xmax,
      mult = mult,
      bw = as.numeric(bw),
      deg = deg
    ),
    family_name = "kde1d",
    types = c("c", "d", "zi")
  )
}

fit_kde1d_margin <- function(
  x,
  weights,
  type,
  xmin,
  xmax,
  mult,
  bw,
  deg
) {
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
    mult = mult,
    bw = bw,
    deg = deg,
    weights = weights
  )
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
  margin_family(
    fit = fit_univariateML_margin,
    fit_args = list(family = family),
    family_name = family,
    types = if (identical(support_type, "Z")) "d" else "c"
  )
}

fit_univariateML_margin <- function(
  x,
  weights,
  type,
  family
) {
  check_univariateML()
  if (length(weights)) {
    stop(
      "univariateML margin families do not support observation weights.",
      call. = FALSE
    )
  }
  fit <- getExportedValue("univariateML", paste0("ml", family))(x)
  attr(fit, "rvine_margin_family") <- family
  attr(fit, "rvine_margin_type") <- type
  attr(fit, "rvine_margin_support") <- unname(qmargin(c(0, 1), fit))
  fit
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
  has_s3_methods(family, c("fit_margin", "margin_info"))
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
  info <- margin_info(family)
  if (!is_margin_info(info, c("family_name", "types"))) {
    stop(
      "margin_info() must return a list with 'family_name' and 'types'.",
      call. = FALSE
    )
  }
  validate_margin_family_name(info$family_name)
  normalized <- normalize_margin_types(info$types, "margin_info()$types")
  if (!length(normalized) || !identical(info$types, normalized)) {
    stop(
      "margin_info()$types must contain supported canonical variable types.",
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

validate_margin_candidate <- function(fit, x, family_name, type) {
  validate_margin(fit, x)
  fit_info <- margin_info(fit)
  if (fit_info$type != type) {
    stop(
      sprintf(
        "fitted family \"%s\" declared type \"%s\", expected \"%s\".",
        family_name,
        fit_info$type,
        type
      ),
      call. = FALSE
    )
  }
  observed <- x[!is.na(x)]
  if (
    any(observed < fit_info$support[1L]) ||
      any(observed > fit_info$support[2L])
  ) {
    stop(
      sprintf(
        "fitted family \"%s\" has support that excludes observed values.",
        family_name
      ),
      call. = FALSE
    )
  }
  invisible(fit)
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
    function(candidate) type %in% margin_info(candidate)$types,
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
  family_names <- vapply(
    candidates,
    function(candidate) margin_info(candidate)$family_name,
    character(1)
  )
  problems <- character()
  fits <- lapply(candidates, function(candidate) {
    family_name <- margin_info(candidate)$family_name
    fit <- tryCatch(
      withCallingHandlers(
        fit_margin(candidate, x, weights = weights, type = type),
        warning = function(condition) {
          problems <<- c(
            problems,
            sprintf(
              "%s warned: %s",
              family_name,
              conditionMessage(condition)
            )
          )
          invokeRestart("muffleWarning")
        }
      ),
      error = function(error) {
        problems <<- c(
          problems,
          sprintf(
            "%s failed: %s",
            family_name,
            conditionMessage(error)
          )
        )
        NULL
      }
    )
    if (is.null(fit)) {
      return(NULL)
    }
    validate_margin_candidate(fit, x, family_name, type)
    fit
  })
  fitted <- !vapply(fits, is.null, logical(1))
  fits <- fits[fitted]
  family_names <- family_names[fitted]
  if (!length(fits)) {
    stop(
      sprintf(
        "could not fit a margin for variable %d (%s).",
        variable,
        paste(unique(problems), collapse = "; ")
      ),
      call. = FALSE
    )
  }
  info <- lapply(fits, margin_info)
  loglik <- vapply(info, `[[`, numeric(1), "loglik")
  valid <- is.finite(loglik)
  problems <- c(
    problems,
    sprintf(
      "%s rejected: fitted margin has no finite log-likelihood",
      family_names[!valid]
    )
  )
  if (!any(valid)) {
    stop(
      sprintf(
        "could not fit a margin for variable %d (%s).",
        variable,
        paste(unique(problems), collapse = "; ")
      ),
      call. = FALSE
    )
  }
  fits <- fits[valid]
  info <- info[valid]
  family_names <- family_names[valid]
  loglik <- loglik[valid]
  npars <- vapply(info, `[[`, numeric(1), "npars")
  if (selcrit != "loglik") {
    valid <- is.finite(npars)
    if (any(valid)) {
      problems <- c(
        problems,
        sprintf(
          "%s rejected: fitted margin has no finite parameter count",
          family_names[!valid]
        )
      )
      fits <- fits[valid]
      family_names <- family_names[valid]
      loglik <- loglik[valid]
      npars <- npars[valid]
    } else if (length(fits) > 1L) {
      stop(
        sprintf(
          paste0(
            "could not select a margin for variable %d because no ",
            "candidate has a finite parameter count."
          ),
          variable
        ),
        call. = FALSE
      )
    } else {
      problems <- c(
        problems,
        sprintf(
          paste0(
            "%s has no finite parameter count; AIC and BIC are ",
            "unavailable"
          ),
          family_names
        )
      )
    }
  }
  criterion <- if (length(fits) == 1L) {
    1
  } else {
    switch(
      selcrit,
      loglik = -loglik,
      aic = -2 * loglik + 2 * npars,
      bic = -2 * loglik + log(sum(!is.na(x))) * npars
    )
  }
  if (length(problems)) {
    warning(
      sprintf(
        "margin selection for variable %d reported problems (%s).",
        variable,
        paste(unique(problems), collapse = "; ")
      ),
      call. = FALSE
    )
  }
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
