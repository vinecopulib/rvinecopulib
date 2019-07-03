context("Truncation generic")

mat <- matrix(c(3, 2, 1, 3, 2, 0, 3, 0, 0), 3, 3)

test_that("works with rvine_structure", {
  struct <- truncate_model(as_rvine_structure(mat), 1)

  expect_equal(unname(dim(struct)), c(3, 1))
  expect_equal(struct$order, 1:3)
  expect_length(struct$struct_array[[1]], 1)
  expect_warning(expect_identical(struct, truncate_model(struct, 2)))
  expect_error(truncate_model(struct, 7))
})

test_that("works with rvine_matrix", {
  mat <- truncate_model(as_rvine_matrix(mat), 1)

  expect_equal(unname(dim(mat)), c(3, 1))
  expect_equal(mat[1, ], c(3, 3, 3))
  expect_equal(mat[2, ], c(0, 2, 0))
  expect_equal(mat[3, ], c(1, 0, 0))
})

test_that("works with vinecop objects", {
  u <- replicate(4, runif(20))
  vc <- vinecop(u, family = "gauss")
  vc_trunc <- truncate_model(vc, 1)

  expect_equal(vc_trunc$npars, 3)
  expect_equal(
    vc_trunc$loglik,
    sum(sapply(vc$pair_copulas[[1]], function(x) x$loglik))
  )
  expect_equal(unname(dim(vc_trunc)), c(4, 1))
  expect_equal(vc_trunc$pair_copulas, vc$pair_copulas[1])
  expect_length(vc_trunc$structure$struct_array[[1]], 1)
  expect_warning(expect_identical(vc_trunc, truncate_model(vc_trunc, 2)))
})

test_that("works with vine objects", {
  x <- replicate(3, runif(20))
  vd <- vine(x, copula_controls = list(family = "gauss"))
  vd_trunc <- truncate_model(vd, 1)

  expect_lt(vd_trunc$npars, vd$npars)
  expect_lt(vd_trunc$loglik, vd$loglik)
  expect_equal(unname(dim(vd_trunc)), c(3, 1))
  expect_equal(vd_trunc$copula$pair_copulas, vd$copula$pair_copulas[1])
  expect_length(vd_trunc$copula$structure$struct_array[[1]], 1)
  expect_warning(expect_identical(vd_trunc, truncate_model(vd_trunc, 2)))
})
