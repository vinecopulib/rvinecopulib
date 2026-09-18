context("Margin protocols")

scored_margin_family <- function(
  family,
  loglik,
  npars,
  fitted_type = "c",
  types = fitted_type,
  error = NULL
) {
  margin_family(
    fit = function(x, weights, type) {
      if (!is.null(error)) {
        stop(error, call. = FALSE)
      }
      margin_dist(
        d = dnorm,
        p = pnorm,
        q = qnorm,
        family = family,
        type = fitted_type,
        npars = npars,
        loglik = loglik
      )
    },
    family_name = family,
    types = types
  )
}

test_that("margin_dist validates and implements the fitted-margin protocol", {
  expect_error(margin_dist(1, pnorm, qnorm), "must be functions")
  expect_error(margin_dist(dnorm, 1, qnorm), "must be functions")
  expect_error(margin_dist(dnorm, pnorm, 1), "must be functions")
  expect_error(margin_dist(dnorm, pnorm, qnorm, family = ""), "family")
  expect_error(margin_dist(dnorm, pnorm, qnorm, type = "mixed"), "type")
  expect_error(
    margin_dist(dnorm, pnorm, qnorm, type = c("c", "d")),
    "length one"
  )
  expect_error(margin_dist(dnorm, pnorm, qnorm, support = c(1, 0)), "support")
  expect_error(margin_dist(dnorm, pnorm, qnorm, npars = -1), "npars")
  expect_error(margin_dist(dnorm, pnorm, qnorm, loglik = Inf), "loglik")
  expect_error(
    margin_dist(function(x) "bad", pnorm, qnorm),
    "dmargin"
  )

  margin <- margin_dist(
    dnorm,
    pnorm,
    qnorm,
    family = "normal",
    support = c(-Inf, Inf),
    npars = 1.5,
    loglik = -4
  )
  expect_s3_class(margin, "margin_dist")
  expect_true(has_margin_protocol(margin))
  expect_identical(as_margin(margin), margin)
  expect_equal(dmargin(0, margin), dnorm(0))
  expect_equal(pmargin(0, margin), 0.5)
  expect_equal(qmargin(0.5, margin), 0)
  expect_equal(
    margin_info(margin),
    list(
      family_name = "normal",
      type = "c",
      support = c(-Inf, Inf),
      npars = 1.5,
      loglik = -4
    )
  )
  expect_equal(as.numeric(logLik(margin)), -4)
  expect_equal(attr(logLik(margin), "df"), 1.5)
})

test_that("stats distributions are fixed fitted-margin implementations", {
  normal <- stats_margin("norm", mean = 1, sd = 2)
  poisson <- stats_margin("pois", lambda = 3)

  expect_equal(dmargin(1, normal), dnorm(1, 1, 2))
  expect_equal(pmargin(1, normal), pnorm(1, 1, 2))
  expect_equal(qmargin(0.5, normal), 1)
  expect_identical(margin_info(normal)$type, "c")
  expect_equal(margin_info(normal)$support, c(-Inf, Inf))
  expect_identical(margin_info(normal)$npars, 2)
  expect_true(is.na(margin_info(normal)$loglik))

  expect_equal(dmargin(2, poisson), dpois(2, 3))
  expect_identical(margin_info(poisson)$type, "d")
  expect_equal(margin_info(poisson)$support, c(0, Inf))
  expect_identical(margin_info(poisson)$npars, 1)
  expect_identical(margin_info(stats_margin("chisq", df = 3))$npars, 1)
  expect_identical(
    margin_info(stats_margin("chisq", df = 3, ncp = 1))$npars,
    2
  )
  expect_identical(margin_info(stats_margin("weibull", shape = 2))$npars, 1)
  expect_error(stats_margin("not-a-distribution"), "unsupported")
  expect_error(stats_margin("gamma"), "shape")
})

test_that("legacy stats lists normalize once through the stats adapter", {
  margin <- as_margin(list(distr = "beta", shape1 = 2, shape2 = 3))
  expect_s3_class(margin, "stats_margin")
  expect_identical(margin_info(margin)$family_name, "beta")
  expect_equal(dmargin(0.5, margin), dbeta(0.5, 2, 3))
  expect_error(as_margin(list(not_a_distribution = TRUE)), "legacy stats")
})

test_that("margin_family validates its fitting contract", {
  fit <- function(x, weights, type) {
    margin_dist(dnorm, pnorm, qnorm, type = type)
  }
  expect_error(margin_family(1), "fit")
  expect_error(margin_family(identity), "x.*weights.*type")
  expect_error(margin_family(fit, family_name = ""), "family")
  expect_error(margin_family(fit, types = character()), "at least one")
  expect_error(margin_family(fit, types = "mixed"), "types")
  expect_error(margin_family(fit, fit_args = list(1)), "uniquely named")
  expect_error(
    margin_family(fit, fit_args = list(type = "c")),
    "cannot contain"
  )
  expect_error(margin_family(fit, fit_args = list(extra = 1)), "accepted")

  family <- margin_family(
    fit,
    family_name = "mixed support",
    types = c("continuous", "discrete", "zero-inflated", "c")
  )
  expect_s3_class(family, "margin_family")
  expect_equal(
    margin_info(family),
    list(family_name = "mixed support", types = c("c", "d", "zi"))
  )
  expect_s3_class(fit_margin(family, 1:3, type = "c"), "margin_dist")

  configured_fit <- function(x, weights, type, location) {
    margin_dist(
      function(y) dnorm(y, location),
      function(y) pnorm(y, location),
      function(p) qnorm(p, location),
      type = type
    )
  }
  configured <- margin_family(configured_fit, fit_args = list(location = 2))
  expect_equal(qmargin(0.5, fit_margin(configured, 1:3)), 2)
})

test_that("fitted-margin validation probes at most the first three observations", {
  margin <- margin_dist(
    d = function(x) {
      if (any(x > 3, na.rm = TRUE)) {
        stop("evaluated too far", call. = FALSE)
      }
      dnorm(x)
    },
    p = function(x) {
      if (any(x > 3, na.rm = TRUE)) {
        stop("evaluated too far", call. = FALSE)
      }
      pnorm(x)
    },
    q = qnorm
  )
  expect_no_error(validate_margin(margin, c(1, 2, NA, 4, 5)))
  expect_no_error(validate_margin(margin, c(1, 2, 3, 4, NA)))

  no_na_propagation <- margin_dist(
    d = function(x) rep(1, length(x)),
    p = function(x) rep(0.5, length(x)),
    q = qnorm
  )
  expect_no_error(validate_margin(no_na_propagation, 1:5))
  expect_error(
    validate_margin(no_na_propagation, c(1, 2, 3, NA)),
    "propagate missing values"
  )
})

dmargin.frontend_test_margin <- function(x, margin) {
  dnorm(x, margin$location, margin$scale)
}

pmargin.frontend_test_margin <- function(x, margin) {
  pnorm(x, margin$location, margin$scale)
}

qmargin.frontend_test_margin <- function(p, margin) {
  qnorm(p, margin$location, margin$scale)
}

margin_info.frontend_test_margin <- function(object) {
  list(
    family_name = "frontend-normal",
    type = object$type,
    support = c(-Inf, Inf),
    npars = 2,
    loglik = object$loglik
  )
}

fit_margin.frontend_test_family <- function(
  family,
  x,
  weights = numeric(),
  type = "c"
) {
  structure(
    list(
      location = if (length(weights)) weighted.mean(x, weights) else mean(x),
      scale = 1,
      type = type,
      loglik = -sum(x^2)
    ),
    class = c("frontend_test_margin", "list")
  )
}

margin_info.frontend_test_family <- function(object) {
  list(family_name = "frontend-normal", types = c("c", "zi"))
}

register_margin_test_methods <- function() {
  for (generic in c(
    "dmargin",
    "pmargin",
    "qmargin",
    "margin_info"
  )) {
    registerS3method(
      generic,
      "frontend_test_margin",
      get(paste0(generic, ".frontend_test_margin")),
      envir = asNamespace("rvinecopulib")
    )
  }
  for (generic in c("fit_margin", "margin_info")) {
    registerS3method(
      generic,
      "frontend_test_family",
      get(paste0(generic, ".frontend_test_family")),
      envir = asNamespace("rvinecopulib")
    )
  }
}

test_that("independent S3 classes can implement both protocols", {
  register_margin_test_methods()
  family <- structure(list(), class = "frontend_test_family")
  expect_true(has_margin_family_protocol(family))
  expect_no_error(validate_margin_family(family))

  margin <- fit_margin(family, c(-1, 1), weights = c(1, 3), type = "c")
  expect_true(has_margin_protocol(margin))
  expect_identical(as_margin(margin), margin)
  expect_equal(margin$location, 0.5)
  expect_equal(
    margin_info(margin),
    list(
      family_name = "frontend-normal",
      type = "c",
      support = c(-Inf, Inf),
      npars = 2,
      loglik = -2
    )
  )
})

test_that("malformed protocol implementations fail near their boundary", {
  expect_error(
    dmargin(0, structure(list(), class = "unknown_margin")),
    "no applicable method.*dmargin"
  )
  expect_error(
    as_margin(structure(list(), class = "unknown_margin")),
    "protocol"
  )

  register_margin_test_methods()
  noncanonical <- structure(
    list(location = 0, scale = 1, type = "continuous", loglik = -1),
    class = "frontend_test_margin"
  )
  expect_error(as_margin(noncanonical), "canonical")
})

test_that("margin types determine central left-limit formulas", {
  continuous <- margin_dist(dnorm, pnorm, qnorm)
  discrete <- stats_margin("pois", lambda = 1)
  zi <- margin_dist(
    d = function(x) ifelse(x == 0, 0.3, ifelse(x > 0, 0.7 * dexp(x), 0)),
    p = function(x) ifelse(x < 0, 0, 0.3 + 0.7 * pexp(x)),
    q = function(p) ifelse(p <= 0.3, 0, qexp(pmax((p - 0.3) / 0.7, 0))),
    type = "zi",
    support = c(0, Inf)
  )

  expect_equal(
    normalize_margin_types(c("cont", "discrete", "zinf")),
    c("c", "d", "zi")
  )
  expect_error(normalize_margin_types("mixed"), "only 'c', 'd', or 'zi'")
  expect_equal(
    pmargin_sub(c(-Inf, 0, Inf, NA), continuous),
    pnorm(c(-Inf, 0, Inf, NA))
  )
  expect_equal(
    pmargin_sub(c(0, 1, Inf, NA), discrete),
    ppois(c(-1, 0, Inf, NA), 1)
  )
  expect_equal(pmargin_sub(c(-1, 0, 1, NA), zi), c(0, 0, pmargin(1, zi), NA))

  ordered_x <- ordered(c("low", "high"), levels = c("low", "high"))
  expect_equal(pmargin_sub(ordered_x, discrete), ppois(c(-1, 0), 1))
})

test_that("zero_inflated is a data-frame-safe type declaration", {
  expect_error(zero_inflated(letters), "numeric")
  x <- zero_inflated(setNames(c(0, NA, 2), c("a", "b", "c")))
  expect_s3_class(x, "zero_inflated")
  expect_equal(unclass(x), setNames(c(0, NA, 2), c("a", "b", "c")))
  data <- data.frame(x = x, y = 1:3)
  expect_s3_class(data$x, "zero_inflated")
  expect_s3_class(data[1:2, , drop = FALSE]$x, "zero_inflated")
})

test_that("kde1d family configuration owns all method-specific controls", {
  expect_error(kde1d_family(xmin = NA), "xmin")
  expect_error(kde1d_family(xmax = NA), "xmax")
  expect_error(kde1d_family(xmin = 1, xmax = 0), "smaller")
  expect_error(kde1d_family(mult = 0), "mult")
  expect_error(kde1d_family(bw = 0), "bw")
  expect_error(kde1d_family(deg = 3), "deg")

  family <- kde1d_family(xmin = 0, xmax = 5, mult = 2, bw = 0.5, deg = 1)
  expect_equal(
    margin_info(family),
    list(family_name = "kde1d", types = c("c", "d", "zi"))
  )
  fit <- fit_margin(family, runif(40, 0, 5), type = "c")
  expect_s3_class(fit, "kde1d")
  expect_equal(margin_info(fit)$support, c(0, 5))
  expect_equal(fit$mult, 2)
  expect_equal(fit$bw, 0.5)
  expect_equal(fit$deg, 1)
})

test_that("family sets normalize aliases, nesting, and variable names", {
  common <- expand_margin_family_set(c("kde1d", "kde1d"), 2)
  expect_equal(length(common), 2)
  expect_equal(lengths(common), c(1L, 1L))
  expect_s3_class(common[[1]][[1]], "custom_margin_family")

  custom <- scored_margin_family("custom", -2, 1)
  nested <- expand_margin_family_set(
    list(second = list("kde1d", custom), first = "nonparametric"),
    2,
    c("first", "second")
  )
  expect_identical(margin_info(nested[[1]][[1]])$family_name, "kde1d")
  expect_identical(margin_info(nested[[2]][[1]])$family_name, "kde1d")
  expect_s3_class(nested[[2]][[2]], "margin_family")

  variants <- normalize_margin_candidates(list(
    kde1d_family(mult = 1),
    kde1d_family(mult = 2),
    kde1d_family(mult = 2)
  ))
  expect_equal(length(variants), 2)
  expect_equal(
    vapply(variants, function(family) family$fit_args$mult, numeric(1)),
    c(1, 2)
  )

  expect_error(expand_margin_family_set(list("kde1d"), 2), "one entry")
  expect_error(expand_margin_family_set(list(), 1), "one entry")
  expect_error(expand_margin_family_set(1, 1), "candidates")
  expect_error(expand_margin_family_set(NA_character_, 1), "cannot contain NA")
  expect_error(expand_margin_family_set(character(), 1), "non-empty")
  expect_error(normalize_margin_candidates(list()), "non-empty")
  expect_error(as_margin_family(1), "names or family objects")
  expect_error(
    expand_margin_family_set(
      list(first = "kde1d", other = "kde1d"),
      2,
      c("first", "second")
    ),
    "each variable name exactly once"
  )
})

test_that("parametric aliases use the optional univariateML adapter", {
  expect_error(
    testthat::with_mocked_bindings(
      expand_margin_family_set("parametric", 1),
      check_univariateML = function() {
        stop("optional adapter unavailable", call. = FALSE)
      },
      .package = "rvinecopulib"
    ),
    "optional adapter unavailable"
  )

  skip_if_not_installed("univariateML")
  expect_error(univariateML_family("not-a-family"), "unknown")
  parametric <- expand_margin_family_set("par", 1)[[1]]
  all_families <- expand_margin_family_set("all", 1)[[1]]
  expect_setequal(
    vapply(
      parametric,
      function(family) margin_info(family)$family_name,
      character(1)
    ),
    univariateML::univariateML_models
  )
  expect_identical(margin_info(all_families[[1]])$family_name, "kde1d")

  # Objects fitted directly by univariateML do not carry rvinecopulib's cached
  # metadata, so the adapter must recover it from the original attributes.
  raw_continuous <- univariateML::mlnorm(seq(-2, 2, length.out = 30))
  raw_discrete <- univariateML::mlpois(rep(0:4, 6))
  expect_identical(margin_info(raw_continuous)$type, "c")
  expect_identical(margin_info(raw_discrete)$type, "d")
  expect_equal(margin_info(raw_continuous)$support, c(-Inf, Inf))
  expect_identical(margin_info(raw_continuous)$family_name, "Normal")
})

test_that("selection uses only generic fitted-margin metadata", {
  simple <- scored_margin_family("simple", -100, 1)
  complex <- scored_margin_family("complex", -97.5, 3)
  select <- function(criterion) {
    select_margin(
      rnorm(100),
      list(simple, complex),
      "c",
      numeric(),
      criterion,
      1
    )
  }

  expect_identical(margin_info(select("loglik"))$family_name, "complex")
  expect_identical(margin_info(select("aic"))$family_name, "complex")
  expect_identical(margin_info(select("bic"))$family_name, "simple")
})

test_that("selection handles incompatibility, candidate failures, and bad fits", {
  good <- scored_margin_family("good", -10, 1)
  failed <- scored_margin_family("failed", -9, 1, error = "optimizer exploded")
  expect_warning(
    selected <- select_margin(
      rnorm(20),
      list(failed, good),
      "c",
      numeric(),
      "aic",
      2
    ),
    "failed.*optimizer exploded"
  )
  expect_identical(
    margin_info(selected)$family_name,
    "good"
  )
  expect_error(
    select_margin(rnorm(20), list(failed), "c", numeric(), "aic", 2),
    "optimizer exploded"
  )

  wrong_type <- scored_margin_family(
    "wrong",
    -1,
    1,
    fitted_type = "d",
    types = "c"
  )
  expect_error(
    select_margin(rnorm(20), list(wrong_type), "c", numeric(), "aic", 1),
    "declared type"
  )

  malformed <- margin_family(
    function(x, weights, type) 1,
    family_name = "malformed"
  )
  expect_error(
    select_margin(rnorm(20), list(malformed), "c", numeric(), "aic", 1),
    "fitted-margin protocol"
  )
  expect_warning(
    selected <- select_margin(
      rnorm(20),
      list(malformed, good),
      "c",
      numeric(),
      "aic",
      1
    ),
    "malformed failed.*fitted-margin protocol"
  )
  expect_identical(margin_info(selected)$family_name, "good")

  missing_loglik <- scored_margin_family("missing-loglik", NA_real_, 1)
  expect_error(
    select_margin(
      rnorm(20),
      list(missing_loglik),
      "c",
      numeric(),
      "aic",
      1
    ),
    "no finite log-likelihood"
  )
  missing_npars <- scored_margin_family("missing-npars", -10, NA_real_)
  expect_warning(
    selected <- select_margin(
      rnorm(20),
      list(missing_npars),
      "c",
      numeric(),
      "aic",
      1
    ),
    "AIC and BIC are unavailable"
  )
  expect_identical(margin_info(selected)$family_name, "missing-npars")
  another_missing_npars <- scored_margin_family(
    "another-missing-npars",
    -11,
    NA_real_
  )
  expect_error(
    select_margin(
      rnorm(20),
      list(missing_npars, another_missing_npars),
      "c",
      numeric(),
      "aic",
      1
    ),
    "no candidate has a finite parameter count"
  )
  expect_identical(
    margin_info(select_margin(
      rnorm(20),
      list(missing_npars),
      "c",
      numeric(),
      "loglik",
      1
    ))$family_name,
    "missing-npars"
  )

  noisy <- margin_family(
    function(x, weights, type) {
      warning("optimizer converged at the boundary")
      fit_margin(good, x, weights, type)
    },
    family_name = "noisy"
  )
  expect_warning(
    select_margin(rnorm(20), list(noisy), "c", numeric(), "aic", 1),
    "noisy warned.*boundary"
  )

  invalid_support <- margin_family(
    function(x, weights, type) {
      margin_dist(
        dnorm,
        pnorm,
        qnorm,
        family = "invalid-support",
        support = c(0, Inf),
        npars = 1,
        loglik = -10
      )
    },
    family_name = "invalid-support"
  )
  expect_error(
    select_margin(-1:1, list(invalid_support), "c", numeric(), "aic", 1),
    "support.*excludes"
  )

  expect_error(
    select_margin(rnorm(20), list(good), "d", numeric(), "aic", 1),
    "no candidate margin"
  )
})

test_that("selection rejects non-finite fitted log-likelihoods", {
  skip_if_not_installed("univariateML")
  x <- c(0, seq(0.1, 2, length.out = 30))
  rayleigh <- univariateML_family("rayleigh")
  good <- scored_margin_family("good", -100, 1)

  expect_error(
    select_margin(x, list(rayleigh), "c", numeric(), "aic", 1),
    "could not fit.*rayleigh failed.*loglik.*finite or NA"
  )
  expect_warning(
    selected <- select_margin(
      x,
      list(rayleigh, good),
      "c",
      numeric(),
      "aic",
      1
    ),
    "rayleigh failed.*loglik.*finite or NA"
  )
  expect_identical(margin_info(selected)$family_name, "good")
})

test_that("every family receives weights and decides how to handle them", {
  x <- c(-2, -1, 0, 3)
  weights <- c(1, 2, 3, 8)
  weighted <- margin_family(
    function(x, weights, type) {
      location <- weighted.mean(x, weights)
      scale <- sqrt(weighted.mean((x - location)^2, weights))
      fit <- margin_dist(
        function(y) dnorm(y, location, scale),
        function(y) pnorm(y, location, scale),
        function(p) qnorm(p, location, scale),
        family = "weighted-normal",
        type = type,
        npars = 2,
        loglik = sum(weights * dnorm(x, location, scale, log = TRUE))
      )
      attr(fit, "fit_weights") <- weights
      fit
    },
    family_name = "weighted-normal"
  )
  selected <- select_margin(x, list(weighted), "c", weights, "aic", 1)
  expect_equal(attr(selected, "fit_weights"), weights)
  expect_equal(qmargin(0.5, selected), weighted.mean(x, weights))

  unsupported <- margin_family(
    function(x, weights, type) {
      if (length(weights)) {
        stop("weights are unsupported", call. = FALSE)
      }
      margin_dist(dnorm, pnorm, qnorm, family = "unsupported", type = type)
    },
    family_name = "unsupported"
  )
  expect_error(
    select_margin(x, list(unsupported), "c", weights, "aic", 1),
    "weights are unsupported"
  )
  x <- rnorm(30)
  weights <- seq_along(x)
  expect_warning(
    selected <- select_margin(
      x,
      list(unsupported, kde1d_family()),
      "c",
      weights,
      "aic",
      1
    ),
    "unsupported failed.*weights"
  )
  expect_s3_class(selected, "kde1d")
})

test_that("univariateML family metadata and weight rejection are explicit", {
  skip_if_not_installed("univariateML")
  normal <- univariateML_family("norm")
  poisson <- univariateML_family("pois")
  expect_identical(margin_info(normal)$types, "c")
  expect_identical(margin_info(poisson)$types, "d")
  fit <- fit_margin(normal, rnorm(30), type = "c")
  expect_s3_class(fit, "univariateML")
  expect_identical(margin_info(fit)$type, "c")
  expect_error(
    fit_margin(normal, rnorm(30), weights = rep(1, 30), type = "c"),
    "do not support observation weights"
  )
})

# Reads a saved margin in a subprocess where univariateML is not attached.
roundtrip_margin_info <- function(path) {
  script <- tempfile(fileext = ".R")
  on.exit(unlink(script), add = TRUE)
  writeLines(
    c(
      "suppressMessages(library(rvinecopulib))",
      sprintf("m <- readRDS(%s)", encodeString(path, quote = '"')),
      "i <- margin_info(m)",
      "cat(i$family_name, i$loglik, i$npars, sep = '\\n')"
    ),
    script
  )
  out <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    shQuote(script),
    stdout = TRUE,
    stderr = FALSE
  ))
  out <- out[nzchar(out)]
  list(
    family_name = out[1],
    loglik = as.numeric(out[2]),
    npars = as.numeric(out[3])
  )
}

test_that("univariateML margins survive a serialization round trip", {
  skip_if_not_installed("univariateML")
  x <- rnorm(200)
  m <- fit_margin(univariateML_family("norm"), x)
  f <- tempfile(fileext = ".rds")
  on.exit(unlink(f), add = TRUE)
  saveRDS(m, f)

  # margin_info() must not depend on univariateML being attached: its logLik()
  # method lives in that namespace, which a fresh session will not have loaded.
  info <- roundtrip_margin_info(f)
  expect_equal(info$family_name, "norm")
  expect_true(is.finite(info$loglik))
  expect_equal(info$npars, 2)
})
