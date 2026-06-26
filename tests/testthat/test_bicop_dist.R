# fixes problems with change in all.equal() behavior in R 4.1.x
expect_eql <- function(...) expect_equal(..., check.environment = FALSE)
expect_equiv <- function(...) expect_equivalent(..., check.environment = FALSE)

context("Class 'bicop_dist'")

test_that("constructor creates proper bicop_dist object", {
  dist <- bicop_dist("gumbel", 90, 3)
  expect_s3_class(dist, "bicop_dist")
  expect_identical(
    names(dist),
    c("family", "rotation", "parameters", "var_types", "npars")
  )
})


test_that("checks for family/rotation/parameters consistency", {
  expect_error(bicop_dist("asdf", 90, 3))
  expect_error(bicop_dist("frank", 3, 3))
  expect_error(bicop_dist("frank", 90, -3))
  expect_error(bicop_dist("frank", 90, 1:2))
})

test_that("partial matching for family names", {
  expect_eql(bicop_dist("ind")$family, "indep")
  expect_eql(bicop_dist("gauss")$family, "gaussian")
  expect_error(bicop_dist("g"))
})

test_that("d/p/r/h functions work", {
  set.seed(0)
  dist <- bicop_dist("bb1", 270, c(1, 2))
  u <- rbicop(50, "bb1", 270, c(1, 2))
  expect_false(any(
    rbicop(50, dist, qrng = FALSE) == rbicop(50, dist, qrng = FALSE)
  ))

  set.seed(1)
  u1 <- rbicop(50, dist, qrng = TRUE)
  set.seed(1)
  expect_true(all(u1 == rbicop(50, dist, qrng = TRUE)))
  expect_gte(min(dbicop(c(0.1, 0.2), dist)), 0)
  expect_gte(min(dbicop(u, dist)), 0)
  expect_gte(min(pbicop(u, dist)), 0)
  expect_lte(max(pbicop(u, dist)), 1)
  expect_lte(max(hbicop(c(0.1, 0.2), 1, dist)), 1)
  expect_lte(max(hbicop(u, 2, dist)), 1)
  expect_lte(max(hbicop(c(0.1, 0.2), 1, dist, inverse = TRUE)), 1)
  expect_lte(max(hbicop(u, 2, dist, inverse = TRUE)), 1)
})

test_that("plot functions work", {
  dist <- bicop_dist("gumbel", 90, 3)
  # we could check some values in the plot objects
  expect_silent(p <- plot(dist))
  expect_silent(p <- contour(dist))
  expect_silent(p <- plot(dist, margins = "norm"))
  expect_silent(p <- contour(dist, margins = "unif"))
  expect_silent(p <- plot(dist, margins = "exp"))
  expect_silent(p <- contour(dist, margins = "flexp"))
})

test_that("parameter <-> tau conversion works", {
  dist <- bicop_dist("joe", 90, 3)

  expect_eql(coef(dist), dist$parameters)

  # one-parameter family
  tau <- par_to_ktau(dist)
  expect_identical(tau, par_to_ktau("joe", 90, 3))

  par <- ktau_to_par(dist, tau)
  expect_identical(par, ktau_to_par("joe", tau))

  expect_eql(3, par[1])

  # two-parameter
  tau <- par_to_ktau("bb1", 0, c(1, 2))
  expect_error(ktau_to_par("bb1", 0.5))

  # rotationless
  ktau_to_par("frank", 0.5)
  ktau_to_par("frank", -0.5)
})

test_that("print method produces output", {
  dist <- bicop_dist("indep")
  expect_output(print(dist))
  expect_output(summary(dist))
})

test_that("getters work", {
  dist <- bicop_dist("gumbel", 90, 3)

  # test get_pair_copula
  expect_identical(dist, get_pair_copula(dist))
  expect_warning(get_pair_copula(dist, 1))
  expect_warning(get_pair_copula(dist, NA, 1))

  # test other getters
  expect_equiv(get_parameters(dist), coef(dist))
  expect_equiv(get_ktau(dist), par_to_ktau(dist))
  expect_equiv(get_family(dist), "gumbel")
})

test_that("works with TLL family", {
  par <- matrix(1, 30, 30)
  expect_eql(bicop_dist("tll", 0, par)$parameters, par)
  expect_error(bicop_dist("tll", 0, par[-1, ]))
  expect_warning(bicop_dist("tll", 0, 2 * par))
})

test_that("vectorized parameters work for dbicop/pbicop/hbicop", {
  set.seed(7)
  n <- 60
  u <- matrix(runif(2 * n), ncol = 2)
  pars <- cbind(seq(-0.8, 0.8, length.out = n))

  d_vec <- dbicop(u, "gaussian", 0, pars)
  p_vec <- pbicop(u, "gaussian", 0, pars)
  h_vec <- hbicop(u, 1, "gaussian", 0, pars)
  hi_vec <- hbicop(u, 2, "gaussian", 0, pars, inverse = TRUE)

  d_row <- vapply(
    seq_len(n),
    function(i) dbicop(u[i, ], "gaussian", 0, pars[i, , drop = FALSE]),
    numeric(1)
  )
  p_row <- vapply(
    seq_len(n),
    function(i) pbicop(u[i, ], "gaussian", 0, pars[i, , drop = FALSE]),
    numeric(1)
  )
  h_row <- vapply(
    seq_len(n),
    function(i) hbicop(u[i, ], 1, "gaussian", 0, pars[i, , drop = FALSE]),
    numeric(1)
  )
  hi_row <- vapply(
    seq_len(n),
    function(i) {
      hbicop(
        u[i, ],
        2,
        "gaussian",
        0,
        pars[i, , drop = FALSE],
        inverse = TRUE
      )
    },
    numeric(1)
  )

  expect_equal(d_vec, d_row)
  expect_equal(p_vec, p_row)
  expect_equal(h_vec, h_row)
  expect_equal(hi_vec, hi_row)
})

test_that("vectorized parameters are rejected by bicop_dist objects", {
  set.seed(8)
  n <- 40
  pars <- cbind(seq(1.1, 4, length.out = n), seq(1.05, 1.9, length.out = n))
  expect_error(
    bicop_dist("bb1", 0, pars),
    "not supported by 'bicop_dist' objects"
  )
})

test_that("vectorized parameters work for rotated two-parameter families", {
  set.seed(10)
  n <- 24
  u <- matrix(runif(2 * n), ncol = 2)
  pars <- cbind(seq(1.1, 3.5, length.out = n), seq(1.05, 1.8, length.out = n))

  d_vec <- dbicop(u, "bb1", 270, pars)
  p_vec <- pbicop(u, "bb1", 270, pars)
  h_vec <- hbicop(u, 2, "bb1", 270, pars)
  hi_vec <- hbicop(u, 1, "bb1", 270, pars, inverse = TRUE)

  d_row <- vapply(
    seq_len(n),
    function(i) dbicop(u[i, ], "bb1", 270, pars[i, , drop = FALSE]),
    numeric(1)
  )
  p_row <- vapply(
    seq_len(n),
    function(i) pbicop(u[i, ], "bb1", 270, pars[i, , drop = FALSE]),
    numeric(1)
  )
  h_row <- vapply(
    seq_len(n),
    function(i) hbicop(u[i, ], 2, "bb1", 270, pars[i, , drop = FALSE]),
    numeric(1)
  )
  hi_row <- vapply(
    seq_len(n),
    function(i) {
      hbicop(
        u[i, ],
        1,
        "bb1",
        270,
        pars[i, , drop = FALSE],
        inverse = TRUE
      )
    },
    numeric(1)
  )

  expect_equal(d_vec, d_row)
  expect_equal(p_vec, p_row)
  expect_equal(h_vec, h_row)
  expect_equal(hi_vec, hi_row)
})

test_that("vectorized parameter sanity checks are enforced", {
  set.seed(9)
  u <- matrix(runif(20), ncol = 2)
  pars_bad_n <- matrix(seq(-0.5, 0.5, length.out = nrow(u) - 1), ncol = 1)
  pars_sim <- matrix(c(-0.3, 0.1), ncol = 1)

  expect_error(
    dbicop(u, "gaussian", 0, pars_bad_n),
    "parameters\\.rows\\(\\) must equal u\\.rows\\(\\)"
  )
  expect_error(
    hbicop(u, 1, "gaussian", 0, pars_bad_n),
    "parameters\\.rows\\(\\) must equal u\\.rows\\(\\)"
  )
  expect_error(
    rbicop(10, "gaussian", 0, pars_sim),
    "not simulation"
  )
})
