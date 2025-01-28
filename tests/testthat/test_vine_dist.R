# fixes problems with change in all.equal() behavior in R 4.1.x
expect_eql <- function(...) expect_equal(..., check.environment = FALSE)
expect_equiv <- function(...) expect_equivalent(..., check.environment = FALSE)

context("Class 'vine_dist'")

set.seed(0)
bicop <- bicop_dist("bb1", 90, c(3, 2))
pcs <- list(list(bicop, bicop), list(bicop))
mat <- matrix(c(1, 2, 3, 1, 2, 0, 1, 0, 0), 3, 3)
vc <- vine_dist(list(list(distr = "norm")), pcs, mat)

test_that("constructor creates proper `vine_dist` object", {
  expect_s3_class(vc, "vine_dist")
  expect_identical(names(vc), c("margins", "copula", "npars", "loglik"))
})


test_that("d/p/r- functions work", {
  u <- rvine(50, vc)
  expect_false(any(rvine(50, vc, qrng = FALSE) == rvine(50, vc, qrng = FALSE)))
  set.seed(1)
  u <- rvine(50, vc, qrng = TRUE)
  set.seed(1)
  expect_true(all(u == rvine(50, vc, qrng = TRUE)))
  expect_gte(min(dvine(u, vc)), 0)
  expect_gte(min(pvine(u, vc, 100)), 0)
  expect_lte(max(pvine(u, vc, 100)), 1)
})

test_that("constructor catches wrong input", {
  # missing margin name
  expect_error(vine_dist(list(stupid = "norm"), pcs, mat))

  # unused margin argument
  expect_error(vine_dist(list(list(distr = "gamma", stupid = 42)), pcs, mat))

  # same with multiple margins specified
  margs <- list(list(distr = "chisq", df = 1),
                list(distr = "exp", scale = 1),  # does not have scale parameter
                list(distr = "lnorm", meanlog = 0, sdlog = 1))
  expect_error(vine_dist(margs, pcs, mat))

  # missing margin argument
  expect_error(vine_dist(list(list(distr = "beta")), pcs, mat)) # shape1 missing

  # incorrect number of argins specified
  margs <- list(distr = "f", df1 = 1, df2 = 1)
  margs <- list(margs, margs, margs, list(distr = "unif", min = 0, max = 1))
  expect_error(vine_dist(margs, cop, mat))
  expect_error(vine_dist(margs[1:2], cop, mat))

  # check npar calculation
  margs <- list(
    list(distr = "weibull", shape = 1, scale = 1),
    list(distr = "t", df = 4, ncp = 1),
    list(distr = "cauchy", location = 0, scale = 1)
  )
  expect_equiv(vine_dist(margs, pcs, mat)$npars, 6 + 6)
})

test_that("print/summary/dim generics work", {
  expect_output(print(vc))

  s <- summary(vc)
  expect_is(s$margins, "data.frame")
  expect_is(s$copula, "data.frame")
  expect_eql(nrow(s$margins), 3)
  expect_eql(ncol(s$margins), 2)
  expect_eql(nrow(s$copula), 3)
  expect_eql(ncol(s$copula), 10)

  expect_equiv(dim(vc)[1], 3)
  expect_equiv(dim(vc)[2], 2)
})

test_that("getters work", {

  # test get_matrix
  expect_equiv(as_rvine_matrix(mat), get_matrix(vc))
  expect_error(get_matrix(12))

  # test get_pair_copulas
  expect_silent(pcc <- get_pair_copula(vc, 1, 1))
  expect_eql(bicop, bicop_dist(pcc$family, pcc$rotation, pcc$parameters))
  expect_error(get_pair_copula(12, 1, 1))
  expect_error(get_pair_copula(vc, 1:2, 1))
  expect_error(get_pair_copula(vc, 1, 1:2))
  expect_error(get_pair_copula(vc, 0, 1))
  expect_error(get_pair_copula(vc, 1, 0))
  expect_error(get_pair_copula(vc, 12, 1))
  expect_error(get_pair_copula(vc, 1, 12))

  # test get_all_pair_copulas
  expect_equiv(pcs, get_all_pair_copulas(vc))
  expect_equiv(pcs[1:2], get_all_pair_copulas(vc, 1:2))
  expect_error(get_all_pair_copulas(12))
  expect_error(get_all_pair_copulas(vc, 0))
  expect_error(get_all_pair_copulas(vc, 12))

  # test other getters
  expect_equiv(get_parameters(vc, 1, 1), coef(pcs[[1]][[1]]))
  expect_equiv(
    get_all_parameters(vc),
    lapply(pcs, function(tree) lapply(tree, coef))
  )
  expect_equiv(get_ktau(vc, 1, 1), par_to_ktau(bicop))
  expect_equiv(
    get_all_ktaus(vc),
    lapply(pcs, function(tree)
      lapply(tree, function(pc) par_to_ktau(pc)))
  )
  expect_equiv(get_family(vc, 1, 1), "bb1")
  expect_equiv(
    get_all_families(vc),
    lapply(pcs, function(tree)
      lapply(tree, function(pc) pc$family))
  )

  # test printed output of getters
  expect_output(print(get_all_pair_copulas(vc)))
  expect_output(print(get_all_pair_copulas(vc, 1)))
})
