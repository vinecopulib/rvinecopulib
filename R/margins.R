#' Marginal distribution protocol
#'
#' Evaluate a fitted marginal distribution. Margin classes can integrate with
#' `rvinecopulib` by implementing these three S3 methods and a [logLik()] method
#' whose return value has a finite `df` attribute.
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
NULL

#' Create a custom marginal distribution
#'
#' `margin_dist()` creates a fitted or fixed margin from three minimal
#' evaluation functions. The functions must each accept one vector argument.
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
  if (!is.character(family) || length(family) != 1L || is.na(family)) {
    stop("'family' must be a single character string.", call. = FALSE)
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
  if (!is.numeric(loglik) || length(loglik) != 1L) {
    stop("'loglik' must be a single number or NA.", call. = FALSE)
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
