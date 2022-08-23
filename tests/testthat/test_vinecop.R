# fixes problems with change in all.equal() behavior in R 4.1.x
expect_eql <- function(...) expect_equal(..., check.environment = FALSE)
expect_equiv <- function(...) expect_equivalent(..., check.environment = FALSE)

context("Fitting 'vinecop' models")

set.seed(5)
u <- sapply(1:5, function(i) runif(30))
fit <- vinecop(u, family = "nonpar")
fit_with_data <- vinecop(u, family = "nonpar", keep_data = TRUE)

test_that("returns proper 'vinecop' object", {
  expect_s3_class(fit, "vinecop")
  expect_s3_class(fit, "vinecop_dist")
  expect_identical(
    names(fit),
    c("pair_copulas", "structure", "var_types",
      "npars", "loglik", "threshold", "controls", "nobs")
  )
  expect_identical(
    names(fit_with_data),
    c("pair_copulas", "structure", "var_types",
      "npars", "loglik", "threshold", "data", "controls", "nobs")
  )

  colnames(u) <- paste(seq_len(ncol(u)))
  expect_identical(
    names(vinecop(u, family = "indep")),
    c("pair_copulas", "structure", "var_types",
      "npars", "loglik", "threshold", "names", "controls", "nobs")
  )
})

test_that("works with structure", {
  u <- sapply(1:2, function(i) runif(30))
  expect_silent(fit <- vinecop(u, structure = matrix(c(1:2, 1:0), 2, 2)))
})

if (Sys.info()["sysname"] != "SunOS") {
  test_that("runs in parallel", {
    expect_silent(fit <- vinecop(u, cores = 2))
  })
}

test_that("S3 generics work", {
  expect_eql(predict(fit, u), fitted(fit_with_data))
  expect_eql(
    predict(fit, u, what = "cdf"),
    fitted(fit_with_data, what = "cdf"),
    tolerance = 0.01
  )
  expect_error(predict(fit, u, what = "hfunc1"))
  fit$data <- NULL
  expect_error(fitted(fit))
  expect_length(attr(logLik(fit), "df"), 1)
})

test_that("print/summary generics work", {
  expect_output(print(fit))
  expect_s3_class(s <- summary(fit), "summary_df")
  expect_is(s, "data.frame")
})

test_that("truncation works", {
  fit_truncated <- truncate_model(fit, trunc_lvl = 1)
  expect_silent(dvinecop(u, fit_truncated))
  expect_silent(rvinecop(50, fit_truncated))

  fit_truncated <- vinecop(u, family = "par", trunc_lvl = 1)
  expect_silent(dvinecop(u, fit_truncated))
  expect_silent(rvinecop(50, fit_truncated))
})

test_that("partial selection works", {
  fit_partial <- vinecop(u[, sample(1:5)],
                         structure = truncate_model(fit$structure, 1),
                         trunc_lvl = 3)
  expect_eql(unname(dim(fit_partial)[2]), 3)

  m_old <- as_rvine_matrix(fit$structure)
  m_new <- as_rvine_matrix(fit_partial$structure)
  tree1_old_edges <- c(paste(diag(m_old[5:2, ]), m_old[1, -5]),
                       paste(m_old[1, -5], diag(m_old[5:2, ])))
  expect_true(all(paste(diag(m_new[5:2, ]), m_new[1, -5]) %in% tree1_old_edges))
})

test_that("d = 1 works", {
  vc <- vinecop(runif(20), structure = rvine_structure(1))
  vc2 <- vinecop(runif(20))
  expect_identical(vc, vc2)

  expect_eql(AIC(vc), 0)
  expect_eql(mBICV(vc), 0)
  expect_eql(dim(summary(vc))[1], 0)
})
