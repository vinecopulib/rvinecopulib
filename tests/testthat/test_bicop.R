# fixes problems with change in all.equal() behavior in R 4.1.x
expect_eql <- function(...) expect_equal(..., check.environment = FALSE)
expect_equiv <- function(...) expect_equivalent(..., check.environment = FALSE)

context("Fitting 'bicop' models")

set.seed(0)
dist <- bicop_dist("gumbel", 90, 3)
u <- rbicop(20, dist)
fit <- bicop(u, keep_data = TRUE)

test_that("returns proper 'bicop' object", {
  expect_s3_class(fit, "bicop")
  expect_s3_class(fit, "bicop_dist")
  expect_identical(
    names(fit),
    c(
      "family",
      "rotation",
      "parameters",
      "var_types",
      "npars",
      "loglik",
      "data",
      "controls",
      "nobs"
    )
  )

  fit <- bicop(u, family = "tll", keep_data = FALSE)
  expect_s3_class(fit, "bicop")
  expect_s3_class(fit, "bicop_dist")
  expect_identical(
    names(fit),
    c(
      "family",
      "rotation",
      "parameters",
      "var_types",
      "npars",
      "loglik",
      "controls",
      "nobs"
    )
  )

  colnames(u) <- paste(1:2)
  expect_identical(
    names(bicop(u, family = "indep")),
    c(
      "family",
      "rotation",
      "parameters",
      "var_types",
      "npars",
      "loglik",
      "names",
      "controls",
      "nobs"
    )
  )
})

test_that("TLL fitting is equivariant to swapping arguments", {
  set.seed(42)
  x <- matrix(rnorm(400), 200, 2)
  x[, 2] <- 0.8 * x[, 1] + sqrt(1 - 0.8^2) * x[, 2]
  u <- pseudo_obs(x)

  fit_12 <- bicop(u, family_set = "tll")
  fit_21 <- bicop(u[, 2:1], family_set = "tll")

  expect_equal(
    dbicop(u, fit_12),
    dbicop(u[, 2:1], fit_21),
    tolerance = 1e-10
  )
  expect_equal(fit_12$loglik, fit_21$loglik, tolerance = 1e-10)
})

test_that("family sets (w/ partial matching)", {
  bicop(u, family = "arch")
  bicop(u, family = "nonp")
  bicop(u, family = "elli")
  bicop(u, family = "onep")
  bicop(u, family = "two")
  bicop(u, family = "bbs")
  bicop(u, family = "itau")
  expect_warning(bicop(u, par_method = "itau"))
  expect_error(bicop(u, family = "asdf"))
})

test_that("as.bicop works", {
  expect_error(as.bicop(list(stupid_argument = 10)))
  object <- list(family = "t", rotation = 0, parameters = c(0.5, 5), npars = 2)
  expect_s3_class(as.bicop(object), "bicop_dist")
  object$var_types <- c("d", "d")
  expect_eql(unlist(as.bicop(object)), unlist(object))

  vec_object <- list(
    family = "gaussian",
    rotation = 0,
    parameters = matrix(seq(-0.5, 0.5, length.out = 10), ncol = 1),
    npars = 10
  )
  expect_error(
    as.bicop(vec_object),
    "not supported by 'bicop_dist' objects"
  )
  unchecked <- as.bicop(vec_object, check = FALSE)
  expect_s3_class(unchecked, "bicop_dist")
  expect_equal(unchecked$parameters, vec_object$parameters)
})


test_that("S3 generics work", {
  expect_eql(predict(fit, u, what = "pdf"), fitted(fit, what = "pdf"))
  expect_eql(predict(fit, u, what = "cdf"), fitted(fit, what = "cdf"))
  expect_eql(predict(fit, u, what = "hfunc1"), fitted(fit, what = "hfunc1"))
  expect_eql(predict(fit, u, what = "hfunc2"), fitted(fit, what = "hfunc2"))
  expect_eql(predict(fit, u, what = "hinv1"), fitted(fit, what = "hinv1"))
  expect_eql(predict(fit, u, what = "hinv2"), fitted(fit, what = "hinv2"))
  u <- as.data.frame(u)
  expect_s3_class(logLik(fit), "logLik")
  expect_equal(
    as.numeric(logLik(fit)),
    sum(log(predict(fit, u, what = "pdf")))
  )
  expect_output(print(fit))
  summary_output <- capture.output(summary(fit))
  expect_length(summary_output, 3)
  expect_match(summary_output[1], "Bivariate copula fit")
  expect_match(summary_output[2], "^Dependence: tau =")
  expect_match(summary_output[2], "; beta =")
  expect_match(summary_output[3], "^Fit: n =")
  expect_match(summary_output[3], "; logLik = .+; df = .+; AIC = .+; BIC =")
  expect_eql(tail_dep(fit), tail_dep(as.bicop(fit)))
  expect_eql(blomqvist_beta(fit), blomqvist_beta(as.bicop(fit)))
  expect_output(print(bicop(u, family = "nonp")))
})
