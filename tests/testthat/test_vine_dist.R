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

test_that("custom fitted margins implement the minimal protocol", {
  margin <- margin_dist(
    dnorm,
    pnorm,
    qnorm,
    family = "normal",
    npars = 2,
    loglik = -10
  )
  custom <- vine_dist(list(margin), pcs, mat)
  legacy <- vine_dist(list(list(distr = "norm")), pcs, mat)

  expect_identical(custom$margins, rep(list(margin), 3))
  expect_equal(dmargin(c(-1, 0, 1), margin), dnorm(c(-1, 0, 1)))
  expect_equal(pmargin(c(-1, 0, 1), margin), pnorm(c(-1, 0, 1)))
  expect_equal(qmargin(c(0.1, 0.5, 0.9), margin), qnorm(c(0.1, 0.5, 0.9)))
  expect_equal(dvine(matrix(0, 1, 3), custom), dvine(matrix(0, 1, 3), legacy))
  expect_equal(summary(custom)$margins$distr, rep("normal", 3))
  expect_equal(custom$npars, legacy$npars)

  path <- tempfile(fileext = ".rds")
  saveRDS(custom, path)
  restored <- readRDS(path)
  unlink(path)
  expect_equal(dvine(matrix(0, 1, 3), restored), dvine(matrix(0, 1, 3), custom))
})

test_that("margin_dist validates its minimal contract", {
  expect_error(margin_dist(1, pnorm, qnorm), "must be functions")
  expect_error(margin_dist(dnorm, pnorm, qnorm, type = "mixed"), "'type'")
  expect_error(margin_dist(dnorm, pnorm, qnorm, npars = -1), "'npars'")
  expect_equal(attr(logLik(margin_dist(dnorm, pnorm, qnorm)), "df"), 0)
})

test_that("left-limit probabilities follow the declared margin type", {
  continuous <- margin_dist(dnorm, pnorm, qnorm)
  discrete <- margin_dist(
    function(x) dpois(x, 1),
    function(x) ppois(x, 1),
    function(p) qpois(p, 1),
    type = "d"
  )
  zi <- margin_dist(
    d = function(x) {
      ifelse(x == 0, 0.25, ifelse(x > 0, 0.75 * dexp(x), 0))
    },
    p = function(x) {
      ifelse(x < 0, 0, ifelse(x == 0, 0.25, 0.25 + 0.75 * pexp(x)))
    },
    q = function(p) {
      ifelse(p <= 0.25, 0, qexp((p - 0.25) / 0.75))
    },
    type = "zi"
  )

  expect_equal(pmargin_sub(c(-1, 0, 1), continuous), pnorm(c(-1, 0, 1)))
  expect_equal(pmargin_sub(c(0, 3), discrete), ppois(c(-1, 2), 1))
  expect_equal(
    pmargin_sub(c(-1, 0, 1), zi),
    c(0, 0, pmargin(1, zi))
  )
})

test_that("legacy fixed margins retain all supported distributions", {
  margins <- list(
    list(distr = "beta", shape1 = 2, shape2 = 3),
    list(distr = "cauchy", location = 0, scale = 2),
    list(distr = "chisq", df = 4),
    list(distr = "exp", rate = 2),
    list(distr = "f", df1 = 4, df2 = 8),
    list(distr = "gamma", shape = 2, rate = 3),
    list(distr = "logis", location = 1, scale = 2),
    list(distr = "lnorm", meanlog = 0, sdlog = 0.5),
    list(distr = "norm", mean = 1, sd = 2),
    list(distr = "t", df = 5),
    list(distr = "unif", min = -1, max = 2),
    list(distr = "weibull", shape = 2, scale = 3)
  )

  for (margin in margins) {
    fixed <- as_margin(margin)
    q <- qmargin(0.37, fixed)
    expect_equal(pmargin(q, fixed), 0.37, tolerance = 1e-8)
    expect_true(is.finite(dmargin(q, fixed)))
  }
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

test_that("rvine performs conditional simulation on the data scale", {
  indep <- bicop_dist()
  vc_cond <- vine_dist(
    rep(list(list(distr = "norm")), 3),
    list(list(indep, indep), list(indep)),
    dvine_structure(1:3)
  )
  vc_cond$names <- vc_cond$copula$names <- c("a", "b", "c")

  x <- rvine(
    10,
    vc_cond,
    x_cond = c(1.25, -0.5),
    conditioning_set = c("c", "b")
  )
  expect_identical(dim(x), c(10L, 3L))
  expect_identical(colnames(x), c("a", "b", "c"))
  expect_equal(x[, "c"], rep(1.25, 10))
  expect_equal(x[, "b"], rep(-0.5, 10))

  x_obs <- cbind(seq(-1, 1, length.out = 5), rep(0.25, 5))
  x <- rvine(5, vc_cond, x_cond = x_obs)
  expect_equal(x[, "b"], x_obs[, 1])
  expect_equal(x[, "c"], x_obs[, 2])

  expect_error(
    rvine(5, vc_cond, conditioning_set = "a"),
    "requires 'x_cond'"
  )
  expect_error(rvine(1.5, vc_cond), "not a count")
  expect_error(rvine(5, vc_cond, qrng = 1), "not a flag")
  expect_error(rvine(5, vc_cond, cores = 0), "not greater than 0")
  expect_error(
    rvine(5, vc_cond, x_cond = list(0.5), conditioning_set = "a"),
    "must be a vector, matrix, or data frame"
  )
  expect_error(
    rvine(5, vc_cond, x_cond = matrix(0.5, 2, 1), conditioning_set = "a"),
    "must have one or 'n' rows"
  )
  expect_error(
    rvine(5, vc_cond, x_cond = c(0.5, 0.5), conditioning_set = "a"),
    "one column per conditioning variable"
  )
  expect_error(
    rvine(5, vc_cond, x_cond = matrix(0.5, 1, 3)),
    "between 1 and d - 1 columns"
  )
  expect_error(
    rvine(5, vc_cond, x_cond = matrix(numeric(), 5, 0)),
    "between 1 and d - 1 columns"
  )
})

test_that("constructor catches wrong input", {
  # missing margin name
  expect_error(vine_dist(list(stupid = "norm"), pcs, mat))

  # unused margin argument
  expect_error(vine_dist(list(list(distr = "gamma", stupid = 42)), pcs, mat))

  # same with multiple margins specified
  margs <- list(
    list(distr = "chisq", df = 1),
    list(distr = "exp", scale = 1), # does not have scale parameter
    list(distr = "lnorm", meanlog = 0, sdlog = 1)
  )
  expect_error(vine_dist(margs, pcs, mat))

  # missing margin argument
  expect_error(vine_dist(list(list(distr = "beta")), pcs, mat)) # shape1 missing

  # incorrect number of argins specified
  margs <- list(distr = "f", df1 = 1, df2 = 1)
  margs <- list(margs, margs, margs, list(distr = "unif", min = 0, max = 1))
  expect_error(vine_dist(margs, cop, mat))
  expect_error(vine_dist(margs[1:2], cop, mat))

  # fixed stats margins retain their distribution parameter counts
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
    lapply(pcs, function(tree) lapply(tree, function(pc) par_to_ktau(pc)))
  )
  expect_equiv(get_family(vc, 1, 1), "bb1")
  expect_equiv(
    get_all_families(vc),
    lapply(pcs, function(tree) lapply(tree, function(pc) pc$family))
  )

  # test printed output of getters
  expect_output(print(get_all_pair_copulas(vc)))
  expect_output(print(get_all_pair_copulas(vc, 1)))
})
