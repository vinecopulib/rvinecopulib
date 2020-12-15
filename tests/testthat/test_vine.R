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
      "margins", "margins_controls", "copula", "copula_controls",
      "npars", "loglik", "data", "weights", "nobs", "names", "var_levels"
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

  fit_truncated <- vine(u, copula_controls = list(
    family_set = "nonpar",
    trunc_lvl = 1
  ))
  expect_silent(dvine(u, fit_truncated))
  expect_silent(rvine(50, fit_truncated))
})

test_that("margins_controls works", {
  fit_mult <- vine(u, margins_controls = list(mult = 2))
  expect_eql(
    sapply(fit_mult$margins, "[[", "bw"),
    2 / log(6) * sapply(fit$margins, "[[", "bw")
  )

  fit_xmin <- vine(abs(u), margins_controls = list(xmin = 0, deg = 1, mult = 1:5))
  expect_eql(sapply(fit_xmin$margins, "[[", "xmin"), rep(0, 5))
  expect_eql(sapply(fit_xmin$margins, "[[", "deg"), rep(1, 5))
})

test_that("weights work", {
  w <- rexp(nrow(u))
  fit_weights <- vine(u, copula_controls = list(family_set = "nonpar"),
                     weights = w, keep_data = TRUE)
  expect_eql(fit_weights$weights, w)
  expect_false(identical(fit$margins[[1]], fit_weights$margins[[1]]))
})

test_that("d = 1 works", {
  vc <- vine(runif(20))
  expect_eql(dim(summary(vc)$margins)[1], 1)
  expect_eql(dim(summary(vc)$copula)[1], 0)
})
