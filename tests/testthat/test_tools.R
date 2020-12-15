# fixes problems with change in all.equal() behavior in R 4.1.x
expect_eql <- function(...) expect_equal(..., check.environment = FALSE)
expect_equiv <- function(...) expect_equivalent(..., check.environment = FALSE)

context("Exported statistical tools")

set.seed(0)

test_that("'pseudo_obs' works", {
  x <- rnorm(10)
  expect_eql(sort(pseudo_obs(x)), (1:10) / 11)
  expect_eql(sort(pseudo_obs(x, lower_tail = FALSE)), (1:10) / 11)
  expect_eql(pseudo_obs(x), 1 - pseudo_obs(x, lower_tail = FALSE))

  mat <- cbind(rnorm(10), rnorm(10))
  expect_eql(sort(pseudo_obs(mat)[, 1]), (1:10) / 11)
  expect_eql(sort(pseudo_obs(mat)[, 2]), (1:10) / 11)

  df <- data.frame(a = rnorm(10), b = rnorm(10))
  expect_eql(sort(pseudo_obs(df)[, 1]), (1:10) / 11)
  expect_eql(sort(pseudo_obs(df)[, 2]), (1:10) / 11)

  expect_is(pseudo_obs(x), "numeric")
  expect_is(pseudo_obs(mat), "matrix")
  expect_is(pseudo_obs(df), "data.frame")
  expect_eql(names(pseudo_obs(df)), c("a", "b"))
  expect_error(pseudo_obs("something"))
  expect_error(pseudo_obs(x, "something"))
})
