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
})

test_that("S3 generics work", {
  expect_eql(predict(fit, u), fitted(fit))
  expect_eql(
    predict(fit, u, what = "cdf"),
    fitted(fit, what = "cdf"),
    tolerance = 0.01
  )
  expect_error(predict(fit, u, what = "hfunc1"))
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
})

test_that("margins_controls works", {
  default_controls <- eval(formals(vine)$margins_controls)
  expect_identical(default_controls$type, "c")

  fit_mult <- vine(u, margins_controls = list(mult = 2))
  expect_eql(
    sapply(fit_mult$margins, "[[", "mult"),
    rep(2, ncol(u))
  )

  fit_xmin <- vine(
    abs(u),
    margins_controls = list(xmin = 0, deg = 1, mult = 1:5)
  )
  expect_eql(sapply(fit_xmin$margins, "[[", "xmin"), rep(0, 5))
  expect_eql(sapply(fit_xmin$margins, "[[", "deg"), rep(1, 5))
})

test_that("weights work", {
  w <- rexp(nrow(u))
  fit_weights <- vine(
    u,
    copula_controls = list(family_set = "nonpar"),
    weights = w,
    keep_data = TRUE
  )
  expect_eql(fit_weights$weights, w)
  expect_false(identical(fit$margins[[1]], fit_weights$margins[[1]]))
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
      margin = list(
        type = c("d", "c", "c"),
        xmin = c(1, NaN, NaN),
        xmax = c(4, NaN, NaN)
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

  expect_no_error(fit <- vine(x, margin = list(type = c("d", "c", "zi"))))
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

  fit <- vine(x, copula_controls = list(family_set = "indep"))
  expect_equal(fit$margins_controls$type, c("c", "zi"))
  expect_equal(fit$copula$var_types, c("c", "d"))

  y <- cbind(rnorm(n), rpois(n, 2))
  fit <- vine(
    y,
    var_types = c("c", "d"),
    copula_controls = list(family_set = "indep")
  )
  expect_equal(fit$margins_controls$type, c("c", "d"))
  expect_error(
    vine(cbind(rnorm(n), runif(n)), var_types = c("c", "d")),
    "integer-valued"
  )
})

test_that("conflicting variable type declarations are rejected", {
  x <- data.frame(a = ordered(rep(1:3, 10)), b = rnorm(30))
  expect_error(vine(x, var_types = c("c", "c")), "disagree")

  z <- data.frame(a = zero_inflated(rexp(30)), b = rnorm(30))
  expect_error(vine(z, var_types = c("d", "c")), "disagree")

  expect_error(
    vine(
      cbind(rnorm(30), rpois(30, 1)),
      var_types = c("c", "d"),
      margins_controls = list(type = c("c", "c"))
    ),
    "disagree"
  )
  expect_error(
    vine(matrix(rnorm(90), ncol = 3), var_types = c("c", "d")),
    "length one or 3"
  )
})
