# fixes problems with change in all.equal() behavior in R 4.1.x
expect_eql <- function(...) expect_equal(..., check.environment = FALSE)
expect_equiv <- function(...) expect_equivalent(..., check.environment = FALSE)

context("Fitting 'vine' models")

set.seed(5)
u <- sapply(1:5, function(i) rnorm(30))
fit <- vine(u, copula_controls = list(family_set = "nonpar"), keep_data = TRUE)

test_that("returns proper 'vine' object", {
  expect_s3_class(fit, "vine")
  expect_s3_class(fit, "vine_dist")
  expect_identical(
    names(fit),
    c(
      "margins",
      "margins_controls",
      "copula",
      "copula_controls",
      "npars",
      "loglik",
      "data",
      "weights",
      "nobs",
      "names",
      "var_levels"
    )
  )
  expect_identical(fit$margins_controls$cores, 1)
})

test_that("S3 generics work", {
  expect_eql(predict(fit, u), fitted(fit))
  expect_eql(
    predict(fit, u, what = "cdf"),
    fitted(fit, what = "cdf"),
    tolerance = 0.01
  )
  expect_error(predict(fit, u, what = "hfunc1"))
  expect_s3_class(logLik(fit), "logLik")
  expect_length(attr(logLik(fit), "df"), 1)
  expect_length(predict(fit, u[1, ], what = "pdf"), 1)
})

test_that("print/summary generics work", {
  expect_output(print(fit))
  s <- summary(fit)
  expect_is(s$margins, c("summary_df", "data.frame"))
  expect_is(s$copula, c("summary_df", "data.frame"))
})

test_that("truncation works", {
  fit_truncated <- truncate_model(fit, trunc_lvl = 1)
  expect_silent(dvine(u, fit_truncated))
  expect_silent(rvine(50, fit_truncated))

  fit_truncated <- vine(
    u,
    copula_controls = list(
      family_set = "nonpar",
      trunc_lvl = 1
    )
  )
  expect_silent(dvine(u, fit_truncated))
  expect_silent(rvine(50, fit_truncated))
})

test_that("conditioning-aware selection is passed to the copula fit", {
  u_named <- as.data.frame(u[, 1:4])
  fit_conditioned <- vine(
    u_named,
    copula_controls = list(
      family_set = "indep",
      conditioning_set = c("V2", "V4")
    )
  )

  expect_identical(
    fit_conditioned$copula$controls$conditioning_set,
    c(2L, 4L)
  )
  expect_equal(
    sort(tail(fit_conditioned$copula$structure$order, 2)),
    c(2L, 4L)
  )

  fit_truncated <- vine(
    u_named,
    copula_controls = list(
      family_set = "indep",
      conditioning_set = c("V2", "V4"),
      trunc_lvl = 1
    )
  )
  expect_length(fit_truncated$copula$pair_copulas, 1L)
  expect_equal(
    sort(tail(fit_truncated$copula$structure$order, 2)),
    c(2L, 4L)
  )
})

test_that("margins_controls works", {
  expect_identical(eval(formals(vine)$margins_controls), list())
  fit_defaults <- vine(
    u[, 1:2],
    margins_controls = list(cores = 1),
    copula_controls = list(family_set = "indep")
  )
  expect_identical(fit_defaults$margins_controls$selcrit, "aic")
  expect_identical(fit_defaults$margins_controls$cores, 1)
  expect_true(all(vapply(
    fit_defaults$margins,
    function(margin) margin_info(margin)$family_name == "kde1d",
    logical(1)
  )))

  fit_mult <- vine(
    u,
    margins_controls = list(family_set = kde1d_family(mult = 2))
  )
  expect_eql(
    sapply(fit_mult$margins, "[[", "mult"),
    rep(2, ncol(u))
  )

  fit_xmin <- vine(
    abs(u),
    margins_controls = list(
      family_set = kde1d_family(xmin = 0, deg = 1)
    )
  )
  expect_eql(sapply(fit_xmin$margins, "[[", "xmin"), rep(0, 5))
  expect_eql(sapply(fit_xmin$margins, "[[", "deg"), rep(1, 5))

  # BEGIN legacy margins_controls compatibility
  legacy_margin_controls_state$warning_issued <- FALSE
  expect_warning(
    fit_legacy <- vine(u, margins_controls = list(mult = 2)),
    "deprecated"
  )
  expect_eql(sapply(fit_legacy$margins, "[[", "mult"), rep(2, ncol(u)))
  expect_no_warning(
    fit_legacy_xmin <- vine(abs(u), margins_controls = list(xmin = 0))
  )
  expect_eql(
    sapply(fit_legacy_xmin$margins, "[[", "mult"),
    rep(log1p(ncol(u)), ncol(u))
  )

  fit_legacy_type <- vine(
    cbind(rnorm(30), rpois(30, 1)),
    margins_controls = list(type = c("c", "d")),
    copula_controls = list(family_set = "indep")
  )
  expect_equal(
    vapply(
      fit_legacy_type$margins,
      function(margin) margin_info(margin)$type,
      character(1)
    ),
    c("c", "d")
  )
  # END legacy margins_controls compatibility
})

test_that("partial copula_controls do not retain data by default", {
  fit_default <- vine(
    u[, 1:2],
    copula_controls = list(family_set = "indep")
  )
  expect_false(fit_default$copula_controls$keep_data)
  expect_null(fit_default$copula$data)

  fit_keep <- vine(
    u[, 1:2],
    copula_controls = list(family_set = "indep", keep_data = TRUE)
  )
  expect_true(fit_keep$copula_controls$keep_data)
  expect_equal(dim(fit_keep$copula$data), c(nrow(u), 2))
})

test_that("weights work", {
  w <- rexp(nrow(u))
  expect_warning(
    fit_weights <- vine(
      u,
      copula_controls = list(family_set = "nonpar"),
      weights = w,
      keep_data = TRUE
    ),
    "AIC and BIC are unavailable"
  )
  expect_eql(fit_weights$weights, w)
  expect_false(identical(fit$margins[[1]], fit_weights$margins[[1]]))
})

test_that("margin fitting can use multiple cores", {
  margin_data <- prep_for_margins(as.data.frame(u))
  family_set <- expand_margin_family_set("kde1d", ncol(u), colnames(u))
  serial <- fit_vine_margins(
    margin_data,
    family_set,
    rep("c", ncol(u)),
    numeric(),
    "aic",
    1
  )
  multicore <- fit_vine_margins(
    margin_data,
    family_set,
    rep("c", ncol(u)),
    numeric(),
    "aic",
    2
  )
  expect_equal(multicore, serial)

  failing <- margin_family(
    function(x, weights, type) stop("margin fit failed", call. = FALSE),
    family_name = "failing"
  )
  expect_error(
    fit_vine_margins(
      margin_data[1:2],
      list(list(failing), family_set[[2]]),
      c("c", "c"),
      numeric(),
      "aic",
      2
    ),
    "margin fit failed"
  )
  expect_error(
    collect_vine_margin_results(list(NULL)),
    "variable 1.*cores = 1"
  )

  noisy <- margin_family(
    function(x, weights, type) {
      warning("margin fit warning")
      margin_dist(
        dnorm,
        pnorm,
        qnorm,
        family = "noisy",
        npars = 1,
        loglik = -10
      )
    },
    family_name = "noisy"
  )
  expect_warning(
    fits <- fit_vine_margins(
      margin_data[1:2],
      rep(list(list(noisy)), 2),
      c("c", "c"),
      numeric(),
      "aic",
      2
    ),
    "margin fit warning"
  )
  expect_length(fits, 2)

  fit_override <- vine(
    u,
    margins_controls = list(cores = 1),
    copula_controls = list(family_set = "indep"),
    cores = 2
  )
  expect_identical(fit_override$margins_controls$cores, 1)
  expect_identical(fit_override$copula_controls$cores, 2)

  fit_parallel_override <- vine(
    u[, 1:2],
    margins_controls = list(cores = 2),
    copula_controls = list(family_set = "indep"),
    cores = 1
  )
  expect_identical(fit_parallel_override$margins_controls$cores, 2)
  expect_identical(fit_parallel_override$copula_controls$cores, 1)
  # as_count() replaces the hand-rolled check and additionally bounds the value
  # by .Machine$integer.max; keep the broader set of invalid inputs
  for (invalid_cores in list(c(2, 3), NA_real_, Inf, 0, -1, 1.5)) {
    expect_error(
      vine(u, margins_controls = list(cores = invalid_cores)),
      "finite positive whole number"
    )
  }
})

test_that("custom tree criteria are available through vine", {
  n_calls <- 0L
  custom_criterion <- function(data, weights) {
    n_calls <<- n_calls + 1L
    abs(cor(data[, 1], data[, 2]))
  }
  fit_custom <- vine(
    as.data.frame(u[, 1:3]),
    copula_controls = list(
      family_set = "indep",
      tree_crit = custom_criterion
    )
  )

  expect_s3_class(fit_custom, "vine")
  expect_gt(n_calls, 0L)
  expect_identical(
    fit_custom$copula_controls$tree_crit,
    custom_criterion
  )
})

test_that("d = 1 works", {
  vc <- vine(runif(20))
  expect_eql(dim(summary(vc)$margins)[1], 1)
  expect_eql(dim(summary(vc)$copula)[1], 0)
})

test_that("discrete variables work", {
  x <- data.frame(
    x1 = rep(1:4, length.out = 50),
    x2 = qnorm((seq_len(50) - 0.5) / 50),
    x3 = floor(((seq_len(50) * 17) %% 53) / 53 * 4)
  )

  expect_no_error(
    fit <- vine(
      x,
      var_types = c("d", "c", "c"),
      margins_controls = list(
        family_set = list(
          kde1d_family(xmin = 1, xmax = 4),
          kde1d_family(),
          kde1d_family()
        )
      )
    )
  )
  x2 <- x
  x2[[1]] <- as.ordered(x2[[1]])
  expect_no_error(fit2 <- vine(x2))

  expect_equal(fit$margins[[1]]$type, "discrete")
  expect_equal(fit$copula$var_types, c("d", "c", "c"))
  expect_equiv(dvine(x, fit), dvine(x2, fit2))
  set.seed(42)
  p <- pvine(x, fit)
  set.seed(42)
  expect_equiv(p, pvine(x2, fit2))
  expect_equal(colnames(rvine(20, fit)), c("x1", "x2", "x3"))

  expect_warning(
    fit <- vine(x, var_types = c("d", "c", "zi")),
    "AIC and BIC are unavailable"
  )
  expect_equal(fit$copula$var_types, c("d", "c", "d"))
  expect_no_error(dvine(x, fit))
  expect_no_error(pvine(x, fit))
  expect_equal(colnames(rvine(20, fit)), c("x1", "x2", "x3"))
})

test_that("variable types can be declared by class or argument", {
  n <- 50
  x <- data.frame(
    continuous = rnorm(n),
    zero = zero_inflated(c(rep(0, 10), rexp(n - 10)))
  )
  expect_s3_class(x$zero, "zero_inflated")
  expect_s3_class(x[1:5, , drop = FALSE]$zero, "zero_inflated")

  expect_warning(
    fit <- vine(x, copula_controls = list(family_set = "indep")),
    "AIC and BIC are unavailable"
  )
  expect_equal(
    vapply(
      fit$margins,
      function(margin) margin_info(margin)$type,
      character(1)
    ),
    c("c", "zi")
  )
  expect_equal(fit$copula$var_types, c("c", "d"))

  y <- cbind(rnorm(n), rpois(n, 2))
  fit <- vine(
    y,
    var_types = c("c", "d"),
    copula_controls = list(family_set = "indep")
  )
  expect_equal(
    vapply(
      fit$margins,
      function(margin) margin_info(margin)$type,
      character(1)
    ),
    c("c", "d")
  )
  expect_error(
    vine(cbind(rnorm(n), runif(n)), var_types = c("c", "d")),
    "integer-valued"
  )
})

test_that("unordered factor expansion preserves missing rows", {
  x <- data.frame(
    category = factor(rep(c("a", "b", "c"), 10)),
    continuous = rnorm(30)
  )
  x$category[c(2, 11)] <- NA

  expanded <- expand_factors(x)
  factor_columns <- setdiff(names(expanded), "continuous")
  expect_equal(nrow(expanded), nrow(x))
  expect_true(all(vapply(expanded[factor_columns], is.ordered, logical(1))))
  expect_true(all(is.na(expanded[c(2, 11), factor_columns])))
  expect_equal(expanded$continuous, x$continuous)

  expect_no_error(
    fit <- vine(x, copula_controls = list(family_set = "indep"))
  )
  expect_equal(fit$nobs, nrow(x))
})

test_that("conflicting variable type declarations are rejected", {
  x <- data.frame(a = ordered(rep(1:3, 10)), b = rnorm(30))
  expect_error(vine(x, var_types = c("c", "c")), "disagree")

  z <- data.frame(a = zero_inflated(rexp(30)), b = rnorm(30))
  expect_error(vine(z, var_types = c("d", "c")), "disagree")

  # BEGIN legacy margins_controls compatibility
  suppressWarnings(
    expect_error(
      vine(
        cbind(rnorm(30), rpois(30, 1)),
        var_types = c("c", "d"),
        margins_controls = list(type = c("c", "c"))
      ),
      "disagree"
    )
  )
  # END legacy margins_controls compatibility
  expect_error(
    vine(matrix(rnorm(90), ncol = 3), var_types = c("c", "d")),
    "length one or 3"
  )
})

test_that("rvinecopulib selects among custom margin families", {
  make_location_family <- function(offset, family) {
    margin_family(
      fit = function(x, weights, type) {
        location <- mean(x) + offset
        margin_dist(
          d = function(y) dnorm(y, location, 1),
          p = function(y) pnorm(y, location, 1),
          q = function(p) qnorm(p, location, 1),
          family = family,
          type = type,
          npars = 1,
          loglik = sum(dnorm(x, location, 1, log = TRUE))
        )
      },
      family_name = family
    )
  }
  good <- make_location_family(0, "good")
  bad <- make_location_family(10, "bad")
  x <- cbind(rnorm(40), rnorm(40, 2))

  fit <- vine(
    x,
    margins_controls = list(
      family_set = list(list(good, bad), list(good, bad)),
      selcrit = "bic"
    ),
    copula_controls = list(family_set = "indep")
  )

  expect_equal(
    vapply(
      fit$margins,
      function(margin) margin_info(margin)$family_name,
      character(1)
    ),
    c("good", "good")
  )
  expect_true(is.finite(as.numeric(logLik(fit))))
  expect_equal(length(fit$margins_controls$family_set), 2)
})

test_that("custom families can fit zero-inflated margins", {
  zi_exponential <- margin_family(
    fit = function(x, weights, type) {
      atom <- mean(x == 0)
      rate <- 1 / mean(x[x > 0])
      margin_dist(
        d = function(y) {
          ifelse(y == 0, atom, (1 - atom) * dexp(y, rate))
        },
        p = function(y) {
          ifelse(y < 0, 0, atom + (1 - atom) * pexp(y, rate))
        },
        q = function(p) {
          ifelse(p <= atom, 0, qexp((p - atom) / (1 - atom), rate))
        },
        family = "zi_exponential",
        type = type,
        npars = 2,
        loglik = sum(ifelse(
          x == 0,
          log(atom),
          log(1 - atom) + dexp(x, rate, log = TRUE)
        ))
      )
    },
    family_name = "zi_exponential",
    types = "zi"
  )
  x <- data.frame(
    zero = zero_inflated(c(rep(0, 20), rexp(60))),
    continuous = rnorm(80)
  )
  fit <- vine(
    x,
    margins_controls = list(
      family_set = list(zi_exponential, "kde1d")
    ),
    copula_controls = list(family_set = "indep")
  )

  expect_equal(margin_info(fit$margins[[1]])$type, "zi")
  expect_equal(pmargin_sub(rep(0, 3), fit$margins[[1]]), rep(0, 3))
  expect_no_error(dvine(x, fit))
})

test_that("univariateML families implement the margin protocol", {
  skip_if_not_installed("univariateML")
  set.seed(12)
  x <- cbind(continuous = rnorm(80), count = rpois(80, 3))
  fit <- vine(
    x,
    margins_controls = list(
      family_set = list(c("norm", "cauchy"), c("pois", "geom")),
      selcrit = "aic"
    ),
    var_types = c("c", "d"),
    copula_controls = list(family_set = "indep")
  )

  expect_true(all(vapply(fit$margins, inherits, logical(1), "univariateML")))
  expect_equal(
    dmargin(x[, 1], fit$margins[[1]]),
    univariateML::dml(x[, 1], fit$margins[[1]])
  )
  expect_equal(
    pmargin(x[, 2], fit$margins[[2]]),
    univariateML::pml(x[, 2], fit$margins[[2]])
  )
  expect_equal(
    pmargin_sub(x[, 2], fit$margins[[2]]),
    univariateML::pml(x[, 2] - 1, fit$margins[[2]])
  )
  expect_no_error(rvine(5, fit))

  set.seed(1)
  ordered_data <- data.frame(
    category = ordered(
      rbinom(200, 4, 0.4),
      levels = 0:4,
      labels = c("very low", "low", "middle", "high", "very high")
    ),
    continuous = rnorm(200)
  )
  ordered_fit <- vine(
    ordered_data,
    margins_controls = list(family_set = list("binom", "norm")),
    copula_controls = list(family_set = "indep")
  )
  expect_no_error(dvine(ordered_data, ordered_fit))
  simulated <- rvine(5, ordered_fit)
  expect_true(is.ordered(simulated$category))
  expect_equal(levels(simulated$category), levels(ordered_data$category))
})

test_that("margin candidate controls are validated", {
  skip_if_not_installed("univariateML")
  x <- cbind(rnorm(30), rpois(30, 2))

  expect_error(
    vine(
      x,
      margins_controls = list(family_set = list("norm", "norm")),
      var_types = c("c", "d")
    ),
    "no candidate margin"
  )
  expect_error(
    vine(
      x[, 1, drop = FALSE],
      margins_controls = list(family_set = "norm", selcrit = "invalid")
    ),
    "must be 'loglik', 'aic', or 'bic'"
  )
  expect_error(
    vine(
      x[, 1, drop = FALSE],
      margins_controls = list(family_set = "norm"),
      weights = rep(1, nrow(x))
    ),
    "univariateML.*do not support observation weights"
  )
})
