context("Covariance and confidence intervals")

set.seed(20260826)

gaussian_vine_fit <- function(n = 600, d = 3, rho = 0.5) {
  S <- matrix(rho, d, d)
  diag(S) <- 1
  u <- pseudo_obs(matrix(rnorm(n * d), n, d) %*% chol(S))
  vinecop(
    u,
    family_set = "gaussian",
    structure = dvine_structure(seq_len(d)),
    keep_data = TRUE
  )
}

test_that("vcov is symmetric, positive definite and correctly named", {
  fit <- gaussian_vine_fit()
  for (step_wise in c(TRUE, FALSE)) {
    V <- vcov(fit, step_wise = step_wise)
    expect_equal(dim(V), c(3L, 3L))
    expect_equal(V, t(V))
    expect_true(all(eigen(V, only.values = TRUE)$values > 0))
    expect_equal(rownames(V), colnames(V))
    expect_true(all(grepl("^T[0-9]", rownames(V))))
  }
})

test_that("standard errors shrink at the parametric rate", {
  small <- vcov(gaussian_vine_fit(n = 500))
  large <- vcov(gaussian_vine_fit(n = 2000))
  ratio <- sqrt(diag(small)) / sqrt(diag(large))
  # quadrupling n should roughly halve the standard errors
  expect_true(all(ratio > 1.5 & ratio < 2.6))
})

test_that("confint brackets the estimates and widens with the level", {
  fit <- gaussian_vine_fit()
  est <- unlist(lapply(fit$pair_copulas, function(t) {
    lapply(t, function(p) p$parameters)
  }))
  ci90 <- confint(fit, level = 0.90)
  ci99 <- confint(fit, level = 0.99)
  expect_true(all(ci90[, 1] < est & est < ci90[, 2]))
  expect_true(all(ci99[, 1] < ci90[, 1]))
  expect_true(all(ci99[, 2] > ci90[, 2]))
  expect_equal(colnames(ci90), c("5 %", "95 %"))
})

test_that("confint covers the truth for a Gaussian vine", {
  # equicorrelated normal: tree-one parameters are rho, tree two is the
  # partial correlation rho / (1 + rho)
  rho <- 0.5
  fit <- gaussian_vine_fit(n = 4000, rho = rho)
  ci <- confint(fit)
  truth <- c(rho, rho, rho / (1 + rho))
  expect_true(all(ci[, 1] <= truth & truth <= ci[, 2]))
})

test_that("parm selects a subset", {
  fit <- gaussian_vine_fit()
  ci <- confint(fit)
  expect_equal(nrow(confint(fit, parm = 1)), 1L)
  expect_equal(rownames(confint(fit, parm = rownames(ci)[2])), rownames(ci)[2])
})

test_that("vcov works for bivariate copulas", {
  u <- pseudo_obs(
    matrix(rnorm(500 * 2), 500, 2) %*% chol(matrix(c(1, .6, .6, 1), 2))
  )
  bf <- bicop(u, family_set = "gaussian", keep_data = TRUE)
  V <- vcov(bf)
  expect_equal(dim(V), c(1L, 1L))
  expect_true(V[1, 1] > 0)
  expect_equal(nrow(confint(bf)), 1L)
})

test_that("parameter labels and parameter-free models handle edge cases", {
  bb1 <- bicop_dist("bb1", parameters = c(2, 1))
  expect_equal(par_labels(bb1), c("bb1.par1", "bb1.par2"))

  u <- matrix(runif(200), ncol = 2)
  independent <- vinecop(u, family_set = "indep", keep_data = TRUE)
  expect_equal(vcov(independent), matrix(numeric(), 0, 0))
})

test_that("there is no method for vine distributions", {
  # intervals that ignore marginal estimation are anticonservative, and
  # fit_margin() is a user-supplied black box, so no method is provided
  n <- 400
  dat <- data.frame(a = rnorm(n), b = rnorm(n))
  fit <- vine(
    dat,
    margins_controls = list(family_set = "kde1d"),
    copula_controls = list(family_set = "parametric"),
    keep_data = TRUE
  )
  expect_null(getS3method("vcov", "vine_dist", optional = TRUE))
  expect_null(getS3method("confint", "vine_dist", optional = TRUE))

  # the copula component can still be examined directly, with the caveat that
  # it conditions on the fitted margins
  expect_true(is.matrix(vcov(fit$copula, newdata = pseudo_obs(dat))))
})

test_that("vcov refuses nonparametric models and reports missing data", {
  u <- pseudo_obs(matrix(rnorm(400 * 3), 400, 3))
  tll <- vinecop(u, family_set = "tll", keep_data = TRUE)
  expect_error(vcov(tll), "nonparametric")

  nodata <- vinecop(u, family_set = "gaussian")
  expect_error(vcov(nodata), "keep_data")
  expect_silent(vcov(nodata, newdata = u))
})

test_that("newdata overrides the stored data", {
  fit <- gaussian_vine_fit(n = 2000)
  full <- sqrt(diag(vcov(fit)))
  half <- sqrt(diag(vcov(fit, newdata = fit$data[1:500, ])))
  expect_true(all(half > full))
})


test_that("A is a Jacobian under step_wise and a Hessian otherwise", {
  fit <- gaussian_vine_fit(n = 800, d = 4)
  Hs <- hessian(fit$data, fit, step_wise = TRUE)
  Hj <- hessian(fit$data, fit, step_wise = FALSE)
  # the step-wise matrix is block triangular, not symmetric, so it must not be
  # symmetrized before inversion
  expect_gt(max(abs(Hs - t(Hs))), 1e-3)
  expect_equal(Hj, t(Hj), tolerance = 1e-8)
})


test_that("sandwich standard errors track the sampling distribution", {
  # a short Monte Carlo check that the bread is oriented correctly: an
  # incorrectly symmetrized bread inflates the tree-one standard errors
  skip_on_cran()
  set.seed(99)
  d <- 3
  n <- 400
  R <- 150
  st <- dvine_structure(1:d)
  truth <- c(0.8, 0.7, 0.4)
  m0 <- vinecop_dist(
    list(
      list(
        bicop_dist("gaussian", 0, truth[1]),
        bicop_dist("gaussian", 0, truth[2])
      ),
      list(bicop_dist("gaussian", 0, truth[3]))
    ),
    st
  )
  est <- se <- matrix(NA, R, 3)
  for (r in seq_len(R)) {
    u <- rvinecop(n, m0)
    f <- vinecop(u, family_set = "gaussian", structure = st, keep_data = TRUE)
    est[r, ] <- unlist(lapply(f$pair_copulas, function(t) {
      lapply(t, function(p) p$parameters)
    }))
    se[r, ] <- sqrt(diag(vcov(f)))
  }
  ratio <- colMeans(se) / apply(est, 2, sd)
  expect_true(all(abs(ratio - 1) < 0.12))
})
