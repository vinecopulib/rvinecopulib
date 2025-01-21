#' Internal: Turn vector input into a matrix with two columns
#'
#' @param u input data
#' @param to_col if `u` is a vector, then `to_col = FALSE` (respectively
#' `to_col = TRUE`) transforms it into a matrix with a single row (respectively
#' single column)
#'
#'
#' @return either a matrix, or an error if u is neither a matrix, data.frame,
#' or a vector
#'
#' @noRd
if_vec_to_matrix <- function(u, to_col = FALSE) {
  if (is.null(u))
    return(NULL)
  assert_that(is.numeric(u) | is.data.frame(u))
  if (NCOL(u) == 1) {
    if (to_col) {
      u <- matrix(u, length(u), 1)
    } else {
      u <- matrix(u, 1, length(u))
    }
  }
  if (!is.matrix(u)) {
    u <- as.matrix(u)
  }

  u
}

#' Internal: Convert arguments to `bicop_dist` object.
#' @param family the family as passed in function call.
#' @param rotation the rotation as passed in function call.
#' @param parameters the parameters as passed in function call.
#' @return A `bicop_dist` object.
#' @noRd
args2bicop <- function(family, rotation, parameters, var_types = c("c", "c")) {
  if (all(inherits(family, "bicop_dist"))) {
    return(family)
  } else {
    if (missing(rotation)) {
      rotation <- 0
    }
    if (missing(parameters)) {
      parameters <- numeric(0)
    }
    assert_that(is.string(family), is.number(rotation), is.numeric(parameters))
    return(bicop_dist(family, rotation, parameters, var_types))
  }
}

process_family_set <- function(family_set, par_method) {
  family_set <- check_and_match_family_set(family_set)
  family_set <- expand_family_set(family_set)
  if (par_method == "itau") {
    if (any(!(family_set %in% family_set_itau))) {
      warning("Only families (",
        paste(family_set_itau, collapse = ", "),
        ") can be used with ", "'par_method = ", '"itau"', "'; ",
        "reducing family set.",
        call. = FALSE
      )
      family_set <- intersect(family_set, family_set_itau)
    }
  }

  family_set
}

#' Internal: Expand shortcuts in the familyset.
#' @noRd
expand_family_set <- function(family_set) {
  unique(unlist(lapply(family_set, expand_family)))
}

expand_family <- function(family) {
  switch(
    family,
    "archimedean" = family_set_archimedean,
    "elliptical" = family_set_elliptical,
    "ev" = family_set_extreme_value,
    "bbs" = family_set_bb,
    "oneparametric" = family_set_onepar,
    "twoparametric" = family_set_twopar,
    "threeparametric" = family_set_threepar,
    "parametric" = family_set_parametric,
    "nonparametric" = family_set_nonparametric,
    "itau" = family_set_itau,
    "all" = family_set_all,
    family # default is no expansion
  )
}

#' Internal: Checks whether all families are known (including partial matching).
#' @noRd
check_and_match_family_set <- function(family_set) {
  matched_fams <- family_set_all_defs[pmatch(family_set, family_set_all_defs)]
  if (any(is.na(matched_fams))) {
    stop(
      "unknown families in family_set: ",
      paste0('"', family_set[is.na(matched_fams)], '"', collapse = ", ")
    )
  }
  matched_fams
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot
# objects) - cols:   Number of columns in layout - layout: A matrix specifying
# the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE), then
# plot 1 will go in the upper left, 2 will go in the upper right, and 3 will go
# all the way across the bottom.
#
multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL) {
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots <- length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots / cols)),
      ncol = cols, nrow = ceiling(numPlots / cols)
    )
  }

  if (numPlots == 1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(
      layout =
        grid::grid.layout(
          nrow(layout),
          ncol(layout)
        )
    ))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]],
        vp = grid::viewport(
          layout.pos.row = matchidx$row,
          layout.pos.col = matchidx$col
        )
      )
    }
  }
}

# Get the depth of a list
depth <- function(this) ifelse(is.list(this), 1L + max(sapply(this, depth)), 0L)

supported_distributions <- c(
  "beta", "cauchy", "chisq", "exp", "f", "gamma",
  "logis", "lnorm", "norm", "t", "unif", "weibull"
)

#' @importFrom stats pbeta qbeta qbeta dcauchy pcauchy qcauchy dchisq pchisq
#' @importFrom stats qchisq dexp pexp qexp df pf qf dgamma pgamma qgamma
#' @importFrom stats dlnorm plnorm qlnorm dt pt qt dunif punif qunif
#' @importFrom stats dweibull pweibull qweibull
check_distr <- function(distr) {

  ## if provided with a kde1d object, then there is nothing to do
  if (inherits(distr, "kde1d")) {
    return(TRUE)
  }

  ## basic sanity checks
  if (!is.list(distr)) {
    return("a distribution should be a kde1d object or a list")
  }
  if (!any(is.element(names(distr), "distr"))) {
    return("a distribution should be a kde1d object or a list with a 'distr' element")
  }
  nn <- distr[["distr"]]
  if (!is.element(nn, supported_distributions)) {
    return("the provided name does not belong to supported distributions")
  }

  ## check that the provided parameters are consistent with the distribution
  qfun <- get(paste0("q", nn))
  par <- distr[names(distr) != "distr"]
  par$p <- 0.5
  e <- tryCatch(do.call(qfun, par), error = function(e) e)
  if (any(class(e) == "error")) {
    return(e$message)
  }

  return(TRUE)
}

get_npars_distr <- function(distr) {
  switch(distr$distr,
    beta = 2,
    cauchy = 2,
    chisq = ifelse("ncp" %in% names(distr), 2, 1),
    exp = 1,
    f = 3,
    gamma = 2,
    lnorm = 2,
    norm = 2,
    logis = 2,
    t = ifelse("ncp" %in% names(distr), 2, 1),
    unif = 2,
    weibull = ifelse("scale" %in% names(distr), 2, 1)
  )
}

#' @noRd
#' @importFrom assertthat assert_that on_failure<-
#' @importFrom assertthat is.number is.string is.flag is.scalar
in_set <- function(el, set) {
  all(el %in% set)
}


on_failure(in_set) <- function(call, env) {
  paste0(
    deparse(call$el),
    " must be one of {",
    paste0(eval(call$set, env), collapse = ", "),
    "}."
  )
}

correct_var_types <- function(var_types, data) {
  is.character(var_types) && in_set(var_types, c("c", "d"))
}

on_failure(correct_var_types) <- function(call, env) {
  paste0("var_types must be vector with elements 'c' or 'd'.")
}

#' Pseudo-Observations
#'
#' Compute the pseudo-observations for the given data matrix.
#'
#' @param x vector or matrix random variates to be converted (column wise) to
#' pseudo-observations.
#' @param ties_method similar to `ties.method` of [rank()] (only `"average"`,
#' `"first"` and `"random"` currently available).
#' @param lower_tail `logical` which, if `FALSE``, returns the
#'   pseudo-observations when applying the empirical marginal survival
#'   functions.
#' @details
#' Given `n` realizations \eqn{x_i=(x_{i1}, \ldots,x_{id})},
#' \eqn{i \in \left\lbrace 1, \ldots,n \right\rbrace }
#' of a random vector `X`, the pseudo-observations are defined via
#' \eqn{u_{ij}=r_{ij}/(n+1)} for
#' \eqn{i \in \left\lbrace 1, \ldots,n \right\rbrace}
#' and
#' \eqn{j \in \left\lbrace 1, \ldots,d \right\rbrace }, where
#' \eqn{r_{ij}} denotes the rank of \eqn{x_{ij}} among all \eqn{x_{kj}},
#' \eqn{k \in \left\lbrace 1, \ldots,n \right\rbrace }.
#'
#' The pseudo-observations can thus also be computed by component-wise applying
#' the empirical distribution functions to the data and scaling the result by
#' \eqn{n/(n+1)}. This asymptotically negligible scaling factor is used to force
#' the variates to fall inside the open unit hypercube, for example, to avoid
#' problems with density evaluation at the boundaries.
#'
#' When `lower_tail = FALSE`, then `pseudo_obs()` simply returns
#' `1 - pseudo_obs()`.
#'
#' @return
#' a vector of matrix of the same dimension as the input containing the
#' pseudo-observations.
#' @examples
#' # pseudo-observations for a vector
#' pseudo_obs(rnorm(10))
#'
#' # pseudo-observations for a matrix
#' pseudo_obs(cbind(rnorm(10), rnorm(10)))
#' @export
pseudo_obs <- function(x, ties_method = "average", lower_tail = TRUE) {
  assert_that(is.scalar(lower_tail) && is.logical(lower_tail))
  assert_that(is.character(ties_method) && is.scalar(ties_method))
  assert_that(in_set(ties_method, c("average", "first", "random")))
  assert_that(is.numeric(x) || is.matrix(x) || is.data.frame(x))

  x[] <- pseudo_obs_cpp(if_vec_to_matrix(x, TRUE), ties_method)
  if (!lower_tail) {
    x <- 1 - x
  }
  x
}

#' Corrected Empirical CDF
#'
#' The empirical CDF with tail correction, ensuring that its output is never
#' 0 or 1.
#'
#' @details The corrected empirical CDF is defined as
#' \deqn{
#' F_n(x) = \frac{1}{n + 1} \max\biggl\{1, \sum_{i = 1}^n 1(X_i \le x)\biggr\}
#' }
#'
#' @param x numeric vector of observations
#'
#' @return A function with signature `function(x)` that returns \eqn{F_n(x)}.
#'
#' @importFrom stats ecdf
#' @export
#' @examples
#' # fit ECDF on simulated data
#' x <- rnorm(100)
#' cdf <- emp_cdf(x)
#'
#' # output is bounded away from 0 and 1
#' cdf(-50)
#' cdf(50)
emp_cdf <- function(x) {
  assert_that(is.numeric(x))
  n <- length(x)
  Fn <- ecdf(x)
  function(xx) pmax(n * Fn(xx), 1) / (n + 1)
}

#' Truncates output of model data frames.
#'
#' @param x a `data.frame` whose print output should be truncated.
#' @noRd
#' @export
print.summary_df <- function(x, ..., rows = 1:10) {
  assert_that(is.numeric(rows), all(rows > 0))
  rows <- intersect(seq_len(nrow(x)), rows)
  x_print <- x[rows, ]
  cat("# A data.frame:", nrow(x), "x", ncol(x), "\n")
  print.data.frame(x_print, digits = 2, row.names = FALSE)
  if (nrow(x) > length(rows)) {
    cat("# ... with", nrow(x) - length(rows), "more rows\n")
  }
  print_varname_legend(x)
  invisible(x)
}

#' internal function
#' @noRd
print_truncation_info <- function(x) {
  if (dim(x)[2] < dim(x)[1] - 1) {
    cat(", ", dim(x)[2], "-truncated", sep = "")
  }
  cat("\n")
}

#' internal function
#' @noRd
print_fit_info <- function(x) {
  ll <- logLik(x)
  cat("nobs =", x$nobs, "  ")
  cat("logLik =", round(ll[1], 2), "  ")
  cat("npars =", round(attr(ll, "df"), 2), "  ")
  cat("AIC =", round(-2 * ll[1] + 2 * attr(ll, "df"), 2), "  ")
  cat("BIC =", round(-2 * ll[1] + log(x$nobs) * attr(ll, "df"), 2), "  ")
  cat("\n")
}

#' internal function
#' @noRd
print_varname_legend <- function(object) {
  # show names if provided
  var_names <- attr(object, "var_names")
  d <- length(var_names)
  if ((d > 0) && any(var_names != paste0("V", 1:d))) {
    linelen <- 80
    d <- length(var_names)
    cat("\n")
    cat("---\n")
    txt <- paste0(1, " <-> ", var_names[1])
    for (i in 2:(d - 1)) {
      if (nchar(txt) > linelen) {
        cat(txt, ",\n", sep = "")
        txt <- paste0(i, " <-> ", var_names[i])
      } else {
        txt <- paste0(txt, ",   ", i, " <-> ", var_names[i])
      }
    }
    if (nchar(txt) > linelen) {
      cat(txt, ",\n", sep = "")
      txt <- paste0(d, " <-> ", var_names[d])
    } else {
      txt <- paste0(txt, ",   ", d, " <-> ", var_names[d])
    }
    cat(txt, "\n")
  }
}

#' internal function : synchronize C++ random number generators with R
#' @importFrom stats runif
#' @noRd
get_seeds <- function() {
  as.numeric(sprintf("%20.0f", runif(20, 1e6, 1e7)))
}


