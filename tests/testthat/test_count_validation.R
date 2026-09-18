# fixes problems with change in all.equal() behavior in R 4.1.x
expect_eql <- function(...) expect_equal(..., check.environment = FALSE)

context("Count validation")

test_that("native counts are finite positive whole numbers", {
  invalid <- c(-1, 0, 1.5, Inf, NA_real_, .Machine$integer.max + 1)
  for (value in invalid) {
    expect_error(as_count(value, "value"), "finite positive whole number")
  }
  expect_identical(as_count(1, "value"), 1)
  expect_identical(
    as_count(.Machine$integer.max, "value"),
    2147483647
  )
})

test_that("simulation and integration counts are checked before native calls", {
  vc <- vinecop_dist(list(list(bicop_dist())), dvine_structure(1:2))
  u <- matrix(0.5, 1, 2)

  expect_error(rbicop(Inf, "indep"), "`n` must be")
  expect_error(rvinecop(-1, vc), "`n` must be")
  expect_error(pvinecop(u, vc, n_mc = -1), "`n_mc` must be")
  expect_error(rvine_structure_sim(1.5), "`d` must be")
})

test_that("thread counts are checked before native calls", {
  cop <- bicop_dist()
  vc <- vinecop_dist(list(list(cop)), dvine_structure(1:2))
  u <- matrix(0.5, 1, 2)

  expect_error(dbicop(u, cop, cores = Inf), "`cores` must be")
  expect_error(dvinecop(u, vc, cores = 1.5), "`cores` must be")
  expect_error(rosenblatt(u, vc, cores = -1), "`cores` must be")
})

test_that("structure indices cannot overflow unsigned native types", {
  expect_error(rvine_structure(Inf), "positive integers")
  expect_error(
    rvine_matrix(matrix(Inf, 1, 1)),
    "finite non-negative integers"
  )
})
