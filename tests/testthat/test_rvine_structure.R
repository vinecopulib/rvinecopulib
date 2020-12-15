# fixes problems with change in all.equal() behavior in R 4.1.x
expect_eql <- function(...) expect_equal(..., check.environment = FALSE)
expect_equiv <- function(...) expect_equivalent(..., check.environment = FALSE)

context("Class 'rvine_structure'")

mat <- matrix(c(4, 3, 2, 1, 4, 3, 2, 0, 4, 3, 0, 0, 4, 0, 0, 0), 4, 4)
mylist <- list(order = 1:4, struct_array = list(c(4, 4, 4), c(3, 3), 2))

test_that("constructor and as/is generics work", {
  expect_silent(rvs <- as_rvine_structure(mylist)) ## calls the constructor
  expect_silent(rvm <- as_rvine_matrix(mylist)) ## true coercion
  expect_silent(as_rvine_structure(mat)) ## calls the constructors
  expect_silent(as_rvine_matrix(mat)) ## true coercion
  expect_silent(as_rvine_matrix(as_rvine_matrix(mat)))
  expect_true(is.rvine_structure(rvs))
  expect_true(is.rvine_matrix(rvm))
  mat[1, 1] <- 0
  expect_silent(as_rvine_matrix(mat, validate = FALSE))
})

test_that("print/dim generics work", {
  rvs <- as_rvine_structure(mylist)
  rvm <- as_rvine_matrix(mat)
  expect_output(print(rvs))
  expect_output(print(rvm))
  expect_equiv(dim(rvs), c(4, 3))
  expect_equiv(dim(rvm), c(4, 3))
})

test_that("C- and D-vine structures work", {
  cv <- cvine_structure(1:6, trunc_lvl = 4)
  expect_equiv(dim(cv), c(6, 4))
  expect_equiv(cv$order, 1:6)
  expect_identical(cv, cvine_structure(6, 4))

  dv <- dvine_structure(1:6, trunc_lvl = 4)
  expect_equiv(dim(dv), c(6, 4))
  expect_equiv(dv$order, 1:6)
  expect_identical(dv, dvine_structure(6, 4))
})

test_that("plot functions work", {
  struct <- as_rvine_structure(mat)
  mat <- rvine_matrix(mat)
  vc <- vinecop_dist(list(replicate(3, bicop_dist(), simplify = FALSE)), struct)
  p <- plot(vc)
  ps <- plot(struct)
  pm <- plot(mat)
  p$plot_env <- NULL
  ps$plot_env <- NULL
  pm$plot_env <- NULL
  expect_equiv(p, ps)
  expect_equiv(p, pm)
})

test_that("d = 1 works", {
  struct <- rvine_structure(1)
  expect_length(struct$struct_array, 0)
  expect_length(struct$order, 1)

  mat <- as_rvine_matrix(struct)
  expect_eql(unname(dim(mat)), c(1, 0))

  expect_output(print(struct))
  expect_error(plot(struct))

  expect_silent(dvine_structure(1))
  expect_silent(cvine_structure(1))
})

