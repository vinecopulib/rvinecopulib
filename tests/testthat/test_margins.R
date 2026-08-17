context("Marginal frontend")

margin_test_controls <- list(
  xmin = NaN,
  xmax = NaN,
  type = "c",
  mult = 1,
  bw = NA,
  deg = 2
)

scored_margin_family <- function(
  family,
  loglik,
  npars,
  type = "c",
  error = NULL
) {
  margin_family(
    fit = function(x) {
      if (!is.null(error)) {
        stop(error, call. = FALSE)
      }
      margin_dist(
        d = dnorm,
        p = pnorm,
        q = qnorm,
        family = family,
        type = type,
        npars = npars,
        loglik = loglik
      )
    },
    family = family,
    type = type
  )
}

test_that("margin_dist validates every constructor field", {
  expect_error(margin_dist(1, pnorm, qnorm), "must be functions")
  expect_error(margin_dist(dnorm, 1, qnorm), "must be functions")
  expect_error(margin_dist(dnorm, pnorm, 1), "must be functions")
  expect_error(margin_dist(dnorm, pnorm, qnorm, family = character()), "family")
  expect_error(margin_dist(dnorm, pnorm, qnorm, family = ""), "family")
  expect_error(
    margin_dist(dnorm, pnorm, qnorm, family = NA_character_),
    "family"
  )
  expect_error(margin_dist(dnorm, pnorm, qnorm, type = "mixed"), "type")
  expect_error(margin_dist(dnorm, pnorm, qnorm, npars = -1), "npars")
  expect_error(margin_dist(dnorm, pnorm, qnorm, npars = Inf), "npars")
  expect_error(margin_dist(dnorm, pnorm, qnorm, loglik = Inf), "loglik")
  expect_error(margin_dist(dnorm, pnorm, qnorm, loglik = NaN), "loglik")

  margin <- margin_dist(
    dnorm,
    pnorm,
    qnorm,
    family = "normal",
    npars = 1.5,
    loglik = -4
  )
  expect_s3_class(margin, "margin_dist")
  expect_equal(as.numeric(logLik(margin)), -4)
  expect_equal(attr(logLik(margin), "df"), 1.5)
  expect_equal(margin_type(margin), "c")
  expect_equal(margin_family_name(margin), "normal")
})

test_that("margin_family validates its minimal fitting contract", {
  expect_error(margin_family(1), "fit")
  expect_error(margin_family(identity, family = character()), "family")
  expect_error(margin_family(identity, family = ""), "family")
  expect_error(margin_family(identity, family = NA_character_), "family")
  expect_error(margin_family(identity, type = character()), "at least one")
  expect_error(margin_family(identity, type = "mixed"), "type")

  family <- margin_family(
    identity,
    family = "mixed support",
    type = c("continuous", "discrete", "zero-inflated", "c")
  )
  expect_s3_class(family, "margin_family")
  expect_equal(family$type, c("c", "d", "zi"))
})

test_that("zero_inflated behaves like a data-frame-safe numeric vector", {
  expect_error(zero_inflated(letters), "numeric")
  x <- zero_inflated(setNames(c(0, NA, 2), c("a", "b", "c")))
  expect_s3_class(x, "zero_inflated")
  expect_equal(unclass(x), setNames(c(0, NA, 2), c("a", "b", "c")))

  data <- data.frame(x = x, y = 1:3)
  expect_s3_class(data$x, "zero_inflated")
  expect_s3_class(data[1:2, , drop = FALSE]$x, "zero_inflated")
  expect_true(is.na(data$x[2]))
})

dmargin.frontend_test_margin <- function(x, margin) {
  dnorm(x, margin[["mean"]], margin[["sd"]])
}

pmargin.frontend_test_margin <- function(x, margin) {
  pnorm(x, margin[["mean"]], margin[["sd"]])
}

qmargin.frontend_test_margin <- function(p, margin) {
  qnorm(p, margin[["mean"]], margin[["sd"]])
}

logLik.frontend_test_margin <- function(object, ...) {
  structure(object$loglik, class = "logLik", df = object$npars)
}

register_frontend_test_margin <- function() {
  registerS3method(
    "dmargin",
    "frontend_test_margin",
    dmargin.frontend_test_margin,
    envir = asNamespace("rvinecopulib")
  )
  registerS3method(
    "pmargin",
    "frontend_test_margin",
    pmargin.frontend_test_margin,
    envir = asNamespace("rvinecopulib")
  )
  registerS3method(
    "qmargin",
    "frontend_test_margin",
    qmargin.frontend_test_margin,
    envir = asNamespace("rvinecopulib")
  )
  registerS3method(
    "logLik",
    "frontend_test_margin",
    logLik.frontend_test_margin,
    envir = asNamespace("stats")
  )
}

test_that("an independent S3 class can implement the fitted-margin protocol", {
  register_frontend_test_margin()
  margin <- structure(
    list(mean = 1, sd = 2, loglik = -10, npars = 2),
    class = "frontend_test_margin",
    type = "c",
    family = "frontend-normal"
  )
  expect_true(has_margin_protocol(margin))
  expect_true(check_distr(margin))
  expect_equal(dmargin(1, margin), dnorm(1, 1, 2))
  expect_equal(pmargin(1, margin), 0.5)
  expect_equal(qmargin(0.5, margin), 1)
  expect_equal(margin_type(margin), "c")
  expect_equal(margin_npars(margin), 2)
  expect_equal(margin_loglik(margin), -10)
  expect_equal(margin_family_name(margin), "frontend-normal")

  model <- vine_dist(
    list(margin),
    list(list(bicop_dist())),
    dvine_structure(1:2)
  )
  expect_equal(summary(model)$margins$distr, rep("frontend-normal", 2))
  expect_no_error(dvine(matrix(c(1, 1), nrow = 1), model))
})

test_that("malformed fitted-margin protocols fail informatively", {
  register_frontend_test_margin()
  expect_error(
    dmargin(0, structure(list(), class = "unknown_margin")),
    "no dmargin"
  )
  expect_error(
    pmargin(0, structure(list(), class = "unknown_margin")),
    "no pmargin"
  )
  expect_error(
    qmargin(0.5, structure(list(), class = "unknown_margin")),
    "no qmargin"
  )
  expect_error(dmargin(0, list(not_a_distribution = TRUE)), "no dmargin")

  no_type <- structure(
    list(mean = 0, sd = 1, loglik = -1, npars = 2),
    class = "frontend_test_margin"
  )
  expect_match(check_distr(no_type), "declare its type")

  bad_df <- structure(
    list(mean = 0, sd = 1, loglik = -1, npars = Inf),
    class = "frontend_test_margin",
    type = "c"
  )
  expect_match(check_distr(bad_df), "finite 'df'")
})

test_that("margin types and left limits cover aliases and boundaries", {
  continuous <- margin_dist(dnorm, pnorm, qnorm, type = "c")
  discrete <- margin_dist(
    function(x) dpois(x, 1),
    function(x) ppois(x, 1),
    function(p) qpois(p, 1),
    type = "d"
  )
  zi <- margin_dist(
    d = function(x) ifelse(x == 0, 0.3, ifelse(x > 0, 0.7 * dexp(x), 0)),
    p = function(x) ifelse(x < 0, 0, 0.3 + 0.7 * pexp(x)),
    q = function(p) ifelse(p <= 0.3, 0, qexp((p - 0.3) / 0.7)),
    type = "zi"
  )

  expect_equal(
    normalize_margin_types(c("cont", "discrete", "zinf")),
    c("c", "d", "zi")
  )
  expect_error(
    normalize_margin_types(c("c", NA_character_)),
    "character vector"
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

  legacy <- list(distr = "norm")
  expect_true(is.na(margin_loglik(legacy)))
  expect_equal(margin_type(legacy), "c")

  expect_equal(
    margin_type(structure(1, class = "other", continuous = TRUE)),
    "c"
  )
  expect_equal(
    margin_type(structure(1, class = "other", continuous = FALSE)),
    "d"
  )
  expect_equal(
    margin_family_name(structure(1, class = "other", model = "model-name")),
    "model-name"
  )
  expect_equal(
    margin_family_name(structure(1, class = "fallback-name")),
    "fallback-name"
  )
})

test_that("family sets normalize common, nested, named, and aliased inputs", {
  common <- expand_margin_family_set(c("kde1d", "kde1d"), 2)
  expect_equal(length(common), 2)
  expect_equal(lengths(common), c(1L, 1L))
  expect_identical(common[[1]][[1]], "kde1d")

  custom <- scored_margin_family("custom", -2, 1)
  nested <- expand_margin_family_set(
    list(second = list("kde1d", custom), first = "nonparametric"),
    2,
    c("first", "second")
  )
  expect_identical(nested[[1]][[1]], "kde1d")
  expect_identical(nested[[2]][[1]], "kde1d")
  expect_s3_class(nested[[2]][[2]], "margin_family")

  expect_error(expand_margin_family_set(list("kde1d"), 2), "one entry")
  expect_error(expand_margin_family_set(list(), 1), "one entry")
  expect_error(expand_margin_family_set(1, 1), "candidates")
  expect_error(expand_margin_family_set(NA_character_, 1), "cannot contain NA")
  expect_error(expand_margin_family_set(character(), 1), "non-empty")
  expect_error(expand_margin_family_set(list(list(1)), 1), "candidates")
  expect_error(normalize_margin_candidates(list()), "non-empty")
  expect_error(
    expand_margin_family_set(
      list(first = "kde1d", other = "kde1d"),
      2,
      c("first", "second")
    ),
    "each variable name exactly once"
  )
  expect_error(
    expand_margin_family_set(
      list(first = "kde1d", first = "kde1d"),
      2,
      c("first", "second")
    ),
    "each variable name exactly once"
  )
})

test_that("parametric aliases use the optional adapter explicitly", {
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
  parametric <- expand_margin_family_set("par", 1)[[1]]
  all_families <- expand_margin_family_set("all", 1)[[1]]
  expect_setequal(unlist(parametric), univariateML::univariateML_models)
  expect_identical(all_families[[1]], "kde1d")
  expect_setequal(unlist(all_families[-1]), univariateML::univariateML_models)
  expect_true(margin_candidate_supports("norm", "c"))
  expect_false(margin_candidate_supports("norm", "d"))
  expect_true(margin_candidate_supports("pois", "d"))
  expect_false(margin_candidate_supports("pois", "zi"))
  expect_error(
    margin_candidate_supports("not-a-family", "c"),
    "unknown marginal family"
  )
})

test_that("log-likelihood, AIC, and BIC selection are owned by rvinecopulib", {
  simple <- scored_margin_family("simple", -100, 1)
  complex <- scored_margin_family("complex", -97.5, 3)
  select <- function(criterion) {
    select_margin(
      rnorm(100),
      list(simple, complex),
      "c",
      margin_test_controls,
      numeric(),
      criterion,
      1
    )
  }

  expect_equal(margin_family_name(select("loglik")), "complex")
  expect_equal(margin_family_name(select("aic")), "complex")
  expect_equal(margin_family_name(select("bic")), "simple")
})

test_that("selection handles candidate failures and invalid fitted objects", {
  good <- scored_margin_family("good", -10, 1)
  failed <- scored_margin_family("failed", -9, 1, error = "optimizer exploded")
  selected <- select_margin(
    rnorm(20),
    list(failed, good),
    "c",
    margin_test_controls,
    numeric(),
    "aic",
    2
  )
  expect_equal(margin_family_name(selected), "good")

  expect_error(
    select_margin(
      rnorm(20),
      list(failed),
      "c",
      margin_test_controls,
      numeric(),
      "aic",
      2
    ),
    "optimizer exploded"
  )

  wrong_type <- scored_margin_family("wrong", -1, 1, type = "d")
  wrong_type$type <- c("c", "d")
  expect_error(
    select_margin(
      rnorm(20),
      list(wrong_type),
      "c",
      margin_test_controls,
      numeric(),
      "aic",
      1
    ),
    "declared type"
  )

  malformed <- margin_family(function(x) 1, family = "malformed")
  expect_error(
    select_margin(
      rnorm(20),
      list(malformed),
      "c",
      margin_test_controls,
      numeric(),
      "aic",
      1
    ),
    "margin protocol"
  )

  no_likelihood <- scored_margin_family("no-likelihood", NA_real_, 1)
  expect_s3_class(
    select_margin(
      rnorm(20),
      list(no_likelihood),
      "c",
      margin_test_controls,
      numeric(),
      "aic",
      1
    ),
    "margin_dist"
  )
  expect_equal(
    margin_family_name(select_margin(
      rnorm(20),
      list(no_likelihood, good),
      "c",
      margin_test_controls,
      numeric(),
      "aic",
      1
    )),
    "good"
  )
  expect_error(
    select_margin(
      rnorm(20),
      list(no_likelihood, no_likelihood),
      "c",
      margin_test_controls,
      numeric(),
      "aic",
      1
    ),
    "finite log-likelihood"
  )
  expect_error(
    select_margin(
      rnorm(20),
      list(good),
      "d",
      margin_test_controls,
      numeric(),
      "aic",
      1
    ),
    "no candidate margin"
  )
})

test_that("observation weights are passed to supporting margin candidates", {
  x <- c(-2, -1, 0, 3)
  weights <- c(1, 2, 3, 8)
  weighted <- margin_family(
    function(x, weights) {
      mu <- weighted.mean(x, weights)
      sigma <- sqrt(weighted.mean((x - mu)^2, weights))
      fit <- margin_dist(
        function(y) dnorm(y, mu, sigma),
        function(y) pnorm(y, mu, sigma),
        function(p) qnorm(p, mu, sigma),
        family = "weighted-normal",
        npars = 2,
        loglik = sum(weights * dnorm(x, mu, sigma, log = TRUE))
      )
      attr(fit, "fit_weights") <- weights
      fit
    },
    family = "weighted-normal"
  )

  selected <- select_margin(
    x,
    list(weighted),
    "c",
    margin_test_controls,
    weights,
    "aic",
    1
  )
  expect_equal(attr(selected, "fit_weights"), weights)
  expect_equal(qmargin(0.5, selected), weighted.mean(x, weights))

  data <- data.frame(first = x, second = x + 1)
  fit <- vine(
    data,
    margins_controls = list(family_set = weighted),
    copula_controls = list(family_set = "indep"),
    weights = weights
  )
  expect_equal(attr(fit$margins[[1]], "fit_weights"), weights)
  expect_equal(attr(fit$margins[[2]], "fit_weights"), weights)
})

test_that("unsupported weighted candidates fail without blocking fallbacks", {
  unsupported <- margin_family(
    function(x) margin_dist(dnorm, pnorm, qnorm, family = "unsupported"),
    family = "unsupported"
  )
  x <- rnorm(30)
  weights <- seq_along(x)

  expect_s3_class(
    select_margin(
      x,
      list(unsupported, "kde1d"),
      "c",
      margin_test_controls,
      weights,
      "aic",
      1
    ),
    "kde1d"
  )
  expect_error(
    select_margin(
      x,
      list(unsupported),
      "c",
      margin_test_controls,
      weights,
      "aic",
      1
    ),
    "unused argument.*weights"
  )

  expect_s3_class(
    select_margin(
      x,
      list(unsupported),
      "c",
      margin_test_controls,
      numeric(),
      "aic",
      1
    ),
    "margin_dist"
  )
})

test_that("the univariateML adapter rejects weights explicitly", {
  skip_if_not_installed("univariateML")
  expect_error(
    select_margin(
      rnorm(20),
      list("norm"),
      "c",
      margin_test_controls,
      rep(1, 20),
      "aic",
      1
    ),
    "univariateML.*do not support observation weights"
  )
})

test_that("default and explicit kde1d paths remain equivalent", {
  set.seed(31)
  x <- cbind(rnorm(40), rexp(40))
  controls <- list(family_set = "indep")
  implicit <- vine(x, copula_controls = controls)
  explicit <- vine(
    x,
    margins_controls = list(family_set = "kde1d"),
    copula_controls = controls
  )
  alias <- vine(
    x,
    margins_controls = list(family_set = "nonpar"),
    copula_controls = controls
  )

  expect_equal(implicit$margins, explicit$margins)
  expect_equal(implicit$margins, alias$margins)
  expect_equal(dvine(x, implicit), dvine(x, explicit))
  expect_equal(dvine(x, implicit), dvine(x, alias))
})

test_that("custom fitted vines support the complete public lifecycle", {
  family <- margin_family(
    fit = function(x) {
      mu <- mean(x)
      sigma <- sqrt(mean((x - mu)^2))
      margin_dist(
        function(y) dnorm(y, mu, sigma),
        function(y) pnorm(y, mu, sigma),
        function(p) qnorm(p, mu, sigma),
        family = "lifecycle-normal",
        npars = 2,
        loglik = sum(dnorm(x, mu, sigma, log = TRUE))
      )
    },
    family = "lifecycle-normal"
  )
  set.seed(32)
  x <- data.frame(first = rnorm(35), second = rnorm(35, 2))
  fit <- vine(
    x,
    margins_controls = list(family_set = family),
    copula_controls = list(family_set = "indep"),
    keep_data = TRUE
  )

  expect_output(print(fit), "vine distribution fit")
  info <- summary(fit)$margins
  expect_equal(info$family, rep("lifecycle-normal", 2))
  expect_equal(info$npars, rep(2, 2))
  expect_true(all(is.finite(info$loglik)))
  expect_true(is.finite(as.numeric(logLik(fit))))
  expect_equal(predict(fit, x, what = "pdf"), fitted(fit, what = "pdf"))
  set.seed(33)
  cdf <- predict(fit, x[1:3, ], what = "cdf", n_mc = 100)
  expect_true(all(cdf >= 0 & cdf <= 1))

  simulated <- rvine(5, fit)
  expect_identical(colnames(simulated), names(x))
  conditioned <- rvine(
    5,
    fit,
    x_cond = 1.25,
    conditioning_set = "second"
  )
  expect_equal(conditioned[, "second"], rep(1.25, 5))

  path <- tempfile(fileext = ".rds")
  saveRDS(fit, path)
  restored <- readRDS(path)
  unlink(path)
  expect_equal(dvine(x, restored), dvine(x, fit))
  set.seed(34)
  original_simulation <- rvine(4, fit)
  set.seed(34)
  expect_equal(rvine(4, restored), original_simulation)
})

test_that("univariateML covers continuous, discrete, and persistence workflows", {
  skip_if_not_installed("univariateML")
  set.seed(35)
  x <- data.frame(
    positive = rlnorm(80, 0, 0.4),
    count = rpois(80, 2)
  )
  fit <- vine(
    x,
    margins_controls = list(
      family_set = list(c("lnorm", "gamma"), c("pois", "geom")),
      selcrit = "bic"
    ),
    var_types = c("c", "d"),
    copula_controls = list(family_set = "indep"),
    keep_data = TRUE
  )

  expect_true(all(vapply(fit$margins, inherits, logical(1), "univariateML")))
  expect_equal(
    dmargin(x$positive, fit$margins[[1]]),
    univariateML::dml(x$positive, fit$margins[[1]])
  )
  expect_equal(
    pmargin(x$count, fit$margins[[2]]),
    univariateML::pml(x$count, fit$margins[[2]])
  )
  expect_equal(
    qmargin(c(0.1, 0.9), fit$margins[[1]]),
    univariateML::qml(c(0.1, 0.9), fit$margins[[1]])
  )
  expect_no_error(summary(fit))
  expect_no_error(predict(fit, x[1:3, ], what = "pdf"))

  path <- tempfile(fileext = ".rds")
  saveRDS(fit, path)
  restored <- readRDS(path)
  unlink(path)
  expect_equal(dvine(x, restored), dvine(x, fit))
})
