context("Class 'rvine_structure'")

mat <- matrix(c(4, 3, 2, 1, 4, 3, 2, 0, 4, 3, 0, 0, 4, 0, 0, 0), 4, 4)
mylist <- list(order = 1:4, struct_array = list(c(4, 4, 4), c(3, 3), 2))

test_that("constructor and as/is generics work", {
  expect_silent(rvs <- as_rvine_structure(mylist)) ## calls the constructor
  expect_silent(rvm <- as_rvine_matrix(mylist)) ## true coercion
  expect_silent(as_rvine_structure(mat)) ## calls the constructors
  expect_silent(as_rvine_matrix(mat)) ## true coercion
  expect_true(is.rvine_structure(rvs))
  expect_true(is.rvine_matrix(rvm))
})

test_that("print/dim generics work", {
  rvs <- as_rvine_structure(mylist)
  rvm <- as_rvine_matrix(mat)
  expect_output(print(rvs))
  expect_output(print(rvm))
  expect_equivalent(dim(rvs), c(4, 3))
  expect_equivalent(dim(rvm), c(4, 3))
})

test_that("C- and D-vine structures work", {
  cv <- cvine_structure(1:6, trunc_lvl = 4)
  expect_equivalent(dim(cv), c(6, 4))
  expect_equivalent(cv$order, 1:6)
  expect_identical(cv, cvine_structure(6, 4))

  dv <- dvine_structure(1:6, trunc_lvl = 4)
  expect_equivalent(dim(dv), c(6, 4))
  expect_equivalent(dv$order, 1:6)
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
  expect_equivalent(p, ps)
  expect_equivalent(p, pm)
})

