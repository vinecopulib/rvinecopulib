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

test_that("dependence measures are computed", {
  dist <- bicop_dist("clayton", 0, 2)
  td <- tail_dep(dist)

  expect_identical(
    dimnames(td),
    list(
      variable1 = c("lower", "upper"),
      variable2 = c("lower", "upper")
    )
  )
  expect_eql(td["lower", "lower"], 2^(-1 / 2))
  expect_eql(td["upper", "upper"], 0)
  expect_eql(blomqvist_beta(dist), 4 / sqrt(7) - 1)
})

test_that("print method produces output", {
  dist <- bicop_dist("indep")
  expect_output(print(dist))
  expect_output(
    summary(dist),
    "Dependence: tau = 0.00; beta = 0.00; tail dependence: none"
  )
  expect_output(
    summary(bicop_dist("clayton", 0, 2)),
    "tail dependence: LL = 0.71"
  )
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

test_that("vectorized parameters work for every one-parameter family (#320)", {
  # Regression guard for the Eigen self-aliasing bug in bicop_wrap(): the
  # earlier Gaussian-only test could not catch it because Gaussian's lower
  # bound (-1) swallows the subnormal garbage, whereas clayton/gumbel/joe
  # (positive lower bounds) throw "parameters exceed lower bound". Loop over
  # all one-parameter families, generating valid per-row parameters generically
  # from a grid of Kendall's taus via ktau_to_par().
  set.seed(11)
  n <- 40
  u <- matrix(runif(2 * n), ncol = 2)
  taus <- seq(0.1, 0.7, length.out = n)

  for (fam in family_set_onepar) {
    pars <- matrix(
      vapply(taus, function(tt) as.numeric(ktau_to_par(fam, tt)), numeric(1)),
      ncol = 1
    )

    d_vec <- dbicop(u, fam, 0, pars)
    p_vec <- pbicop(u, fam, 0, pars)
    h1_vec <- hbicop(u, 1, fam, 0, pars)
    h2_vec <- hbicop(u, 2, fam, 0, pars)
    hi1_vec <- hbicop(u, 1, fam, 0, pars, inverse = TRUE)
    hi2_vec <- hbicop(u, 2, fam, 0, pars, inverse = TRUE)

    by_row <- function(f) vapply(seq_len(n), f, numeric(1))
    d_row <- by_row(function(i) dbicop(u[i, ], fam, 0, pars[i, , drop = FALSE]))
    p_row <- by_row(function(i) pbicop(u[i, ], fam, 0, pars[i, , drop = FALSE]))
    h1_row <- by_row(function(i) {
      hbicop(u[i, ], 1, fam, 0, pars[i, , drop = FALSE])
    })
    h2_row <- by_row(function(i) {
      hbicop(u[i, ], 2, fam, 0, pars[i, , drop = FALSE])
    })
    hi1_row <- by_row(function(i) {
      hbicop(u[i, ], 1, fam, 0, pars[i, , drop = FALSE], inverse = TRUE)
    })
    hi2_row <- by_row(function(i) {
      hbicop(u[i, ], 2, fam, 0, pars[i, , drop = FALSE], inverse = TRUE)
    })

    expect_eql(d_vec, d_row, info = fam)
    expect_eql(p_vec, p_row, info = fam)
    expect_eql(h1_vec, h1_row, info = fam)
    expect_eql(h2_vec, h2_row, info = fam)
    expect_eql(hi1_vec, hi1_row, info = fam)
    expect_eql(hi2_vec, hi2_row, info = fam)
    expect_true(all(is.finite(d_vec)), info = fam)
  }
})

test_that("reported Clayton vectorized-parameter crash is fixed (#320)", {
  set.seed(123)
  n <- 10
  u <- matrix(runif(n * 2), ncol = 2)
  par <- matrix(runif(n, 2, 3), ncol = 1)
  res <- dbicop(u = u, family = "clayton", parameters = par)
  expect_length(res, n)
  expect_true(all(is.finite(res)))
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

  expect_error(
    dbicop(u, "gaussian", 0, pars_bad_n),
    "parameters\\.rows\\(\\) must equal u\\.rows\\(\\)"
  )
  expect_error(
    hbicop(u, 1, "gaussian", 0, pars_bad_n),
    "parameters\\.rows\\(\\) must equal u\\.rows\\(\\)"
  )
})

test_that("rbicop supports vectorized parameters", {
  n <- 40
  pars <- matrix(seq(-0.8, 0.8, length.out = n), ncol = 1)

  set.seed(123)
  uniforms <- rbicop(n, "indep")
  expected <- uniforms
  expected[, 2] <- hbicop(
    uniforms,
    1,
    "gaussian",
    0,
    pars,
    inverse = TRUE
  )

  set.seed(123)
  simulated <- rbicop(n, "gaussian", 0, pars)
  expect_equal(simulated, expected)

  set.seed(123)
  expect_equal(rbicop(seq_len(n), "gaussian", 0, pars), expected)

  expect_error(
    rbicop(n - 1L, "gaussian", 0, pars),
    "one row per row of u"
  )
})

test_that("vectorized rbicop is synchronized with the R seed", {
  n <- 25
  pars <- matrix(seq(-0.8, 0.8, length.out = n), ncol = 1)

  for (use_qrng in c(FALSE, TRUE)) {
    set.seed(314)
    u1 <- rbicop(n, "gaussian", 0, pars, qrng = use_qrng)
    state_after <- .Random.seed

    set.seed(314)
    u2 <- rbicop(n, "gaussian", 0, pars, qrng = use_qrng)
    expect_identical(u2, u1)
    expect_identical(.Random.seed, state_after)

    set.seed(314)
    runif(20)
    expect_identical(.Random.seed, state_after)

    set.seed(315)
    u3 <- rbicop(n, "gaussian", 0, pars, qrng = use_qrng)
    expect_false(identical(u3, u1))
  }
})

test_that("rbicop vectorization handles rotations and multiple parameters", {
  n <- 30
  pars <- cbind(
    seq(1.1, 3.5, length.out = n),
    seq(1.05, 1.8, length.out = n)
  )

  set.seed(456)
  uniforms <- rbicop(n, "indep", qrng = TRUE)
  expected <- uniforms
  expected[, 2] <- hbicop(
    uniforms,
    1,
    "bb1",
    270,
    pars,
    inverse = TRUE
  )

  set.seed(456)
  expect_equal(rbicop(n, "bb1", 270, pars, qrng = TRUE), expected)
})

test_that("rbicop validates vectorized simulation inputs", {
  pars <- matrix(c(-0.3, 0.1), ncol = 1)

  expect_error(rbicop(1.5, "gaussian", 0, 0.2), "not a count")
  expect_error(rbicop(3, "gaussian", 0, pars), "one row per row of u")
  expect_error(
    rbicop(2, "gaussian", 0, matrix(c(0.1, NA_real_), ncol = 1)),
    "must not contain NaN or Inf"
  )
  expect_error(
    rbicop(2, "gaussian", 0, matrix(c(0.1, -1.1), ncol = 1)),
    "out of bounds"
  )
})

test_that("scores and hessian work for bivariate models", {
  set.seed(13)
  u_deriv <- matrix(runif(80, 0.1, 0.9), ncol = 2)
  cop <- bicop_dist("gaussian", parameters = 0.4)

  s <- scores(u_deriv, cop)
  h <- hessian(u_deriv, cop)
  expect_equal(dim(s), c(nrow(u_deriv), 1L))
  expect_equal(dim(h), c(1L, 1L))
  expect_equal(scores(u_deriv, cop, cores = 2), s)
  expect_equal(hessian(u_deriv, cop, cores = 2), h)

  eps <- 1e-5
  log_density <- function(par) log(dbicop(u_deriv, "gaussian", 0, par))
  s_fd <- (log_density(0.4 + eps) - log_density(0.4 - eps)) / (2 * eps)
  h_fd <- mean(
    (log_density(0.4 + eps) - 2 * log_density(0.4) + log_density(0.4 - eps)) /
      eps^2
  )
  expect_equal(drop(s), s_fd, tolerance = 1e-7)
  expect_equal(drop(h), h_fd, tolerance = 1e-5)
})

test_that("bivariate derivatives support per-observation parameters", {
  set.seed(14)
  n <- 30
  u_deriv <- matrix(runif(2 * n, 0.1, 0.9), ncol = 2)
  pars <- matrix(seq(-0.7, 0.7, length.out = n), ncol = 1)
  cop <- bicop_dist("gaussian", parameters = 0)

  s <- scores(u_deriv, cop, parameters = pars)
  h <- hessian(u_deriv, cop, parameters = pars)
  expect_equal(scores(u_deriv, cop, parameters = drop(pars)), s)
  expect_equal(hessian(u_deriv, cop, parameters = drop(pars)), h)
  s_row <- vapply(
    seq_len(n),
    function(i) {
      scores(
        u_deriv[i, ],
        bicop_dist("gaussian", parameters = pars[i, 1])
      )
    },
    numeric(1)
  )
  h_row <- vapply(
    seq_len(n),
    function(i) {
      hessian(
        u_deriv[i, ],
        bicop_dist("gaussian", parameters = pars[i, 1])
      )
    },
    numeric(1)
  )

  expect_equal(drop(s), s_row)
  expect_equal(drop(h), mean(h_row))
})

test_that("bivariate derivative safeguards are enforced", {
  u_deriv <- matrix(c(0.2, 0.3, 0.7, 0.8), ncol = 2)
  cop <- bicop_dist("gaussian", parameters = 0.4)

  expect_error(scores(u_deriv, cop, cores = 0), "not greater than 0")
  expect_error(hessian(u_deriv, cop, cores = 0), "not greater than 0")
  expect_error(scores(u_deriv, cop, parameters = "bad"), "not a numeric")
  expect_error(
    scores(u_deriv, cop, parameters = matrix(0.2, 1, 1)),
    "one row per row of u"
  )
  expect_error(
    hessian(u_deriv, cop, parameters = matrix(c(0.2, NA), 2, 1)),
    "must not contain NaN or Inf"
  )
  expect_error(
    scores(u_deriv, cop, parameters = matrix(c(0.2, 1.1), 2, 1)),
    "out of bounds"
  )
  expect_error(
    scores(
      cbind(u_deriv, u_deriv[, 1] - 0.01),
      bicop_dist(
        "gaussian",
        parameters = 0.4,
        var_types = c("d", "c")
      )
    ),
    "continuous"
  )
  expect_error(
    hessian(u_deriv, bicop_dist("tll")),
    "not implemented for the TLL"
  )
})

test_that("dbicop exposes selected density derivatives", {
  u <- matrix(c(0.21, 0.37, 0.44, 0.62, 0.73, 0.28), ncol = 2)
  cop <- bicop_dist("gaussian", parameters = 0.35)
  eps <- 1e-6

  u_plus <- u_minus <- u
  u_plus[, 1] <- u_plus[, 1] + eps
  u_minus[, 1] <- u_minus[, 1] - eps
  du1_fd <- (dbicop(u_plus, cop) - dbicop(u_minus, cop)) / (2 * eps)

  expect_equal(dbicop(u, cop, deriv = "u1"), du1_fd, tolerance = 1e-6)
  expect_equal(
    dbicop(u, cop, deriv = "par"),
    dbicop(u, cop, deriv = "par1")
  )
  expect_equal(
    dbicop(u, cop, deriv = c("u1", "par1")),
    dbicop(u, cop, deriv = c("par1", "u1"))
  )
  expect_equal(dbicop(u, cop, log = TRUE), log(dbicop(u, cop)))
  expect_equal(
    dbicop(u, cop, log = TRUE, deriv = "par1"),
    drop(scores(u, cop))
  )
  expect_equal(
    mean(dbicop(u, cop, log = TRUE, deriv = c("par1", "par1"))),
    drop(hessian(u, cop))
  )
})

test_that("hbicop exposes selected h-function derivatives", {
  u <- matrix(c(0.18, 0.31, 0.47, 0.66, 0.79, 0.23), ncol = 2)
  cop <- bicop_dist("bb1", rotation = 270, parameters = c(1.2, 1.4))
  eps <- 1e-6

  cop_plus <- bicop_dist("bb1", 270, c(1.2 + eps, 1.4))
  cop_minus <- bicop_dist("bb1", 270, c(1.2 - eps, 1.4))
  dh1_dpar1_fd <-
    (hbicop(u, 1, cop_plus) - hbicop(u, 1, cop_minus)) / (2 * eps)

  expect_equal(
    hbicop(u, 1, cop, deriv = "par1"),
    dh1_dpar1_fd,
    tolerance = 1e-6
  )
  expect_equal(hbicop(u, 1, cop, deriv = "u2"), dbicop(u, cop))
  expect_equal(hbicop(u, 2, cop, deriv = "u1"), dbicop(u, cop))
  expect_equal(
    hbicop(u, 1, cop, deriv = c("par1", "u2")),
    dbicop(u, cop, deriv = "par1")
  )
  expect_equal(
    hbicop(u, 2, cop, deriv = c("u2", "par2")),
    hbicop(u, 2, cop, deriv = c("par2", "u2"))
  )
})

test_that("function derivatives support observation-specific parameters", {
  u <- matrix(c(0.2, 0.35, 0.5, 0.65, 0.8, 0.25), ncol = 2)
  parameters <- matrix(c(-0.5, 0.1, 0.6), ncol = 1)
  evaluators <- list(
    function(x, par, cores = 1) {
      dbicop(x, "gaussian", parameters = par, deriv = "par1", cores = cores)
    },
    function(x, par, cores = 1) {
      dbicop(
        x,
        "gaussian",
        parameters = par,
        deriv = c("u1", "par1"),
        cores = cores
      )
    },
    function(x, par, cores = 1) {
      dbicop(
        x,
        "gaussian",
        parameters = par,
        log = TRUE,
        deriv = "par1",
        cores = cores
      )
    },
    function(x, par, cores = 1) {
      dbicop(
        x,
        "gaussian",
        parameters = par,
        log = TRUE,
        deriv = c("par1", "par1"),
        cores = cores
      )
    },
    function(x, par, cores = 1) {
      hbicop(
        x,
        1,
        "gaussian",
        parameters = par,
        deriv = "par1",
        cores = cores
      )
    },
    function(x, par, cores = 1) {
      hbicop(
        x,
        1,
        "gaussian",
        parameters = par,
        deriv = c("u1", "par1"),
        cores = cores
      )
    },
    function(x, par, cores = 1) {
      hbicop(
        x,
        2,
        "gaussian",
        parameters = par,
        deriv = "par1",
        cores = cores
      )
    },
    function(x, par, cores = 1) {
      hbicop(
        x,
        2,
        "gaussian",
        parameters = par,
        deriv = c("u2", "par1"),
        cores = cores
      )
    }
  )

  for (evaluate in evaluators) {
    expected <- vapply(seq_len(nrow(u)), function(i) {
      evaluate(u[i, ], parameters[i, 1])
    }, numeric(1))
    expect_equal(evaluate(u, parameters, cores = 2), expected)
  }
  expect_equal(
    dbicop(u, "gaussian", parameters = parameters, cores = 2),
    dbicop(u, "gaussian", parameters = parameters)
  )
  expect_equal(
    hbicop(u, 1, "gaussian", parameters = parameters, inverse = TRUE, cores = 2),
    hbicop(u, 1, "gaussian", parameters = parameters, inverse = TRUE)
  )
})

test_that("function derivative safeguards are enforced", {
  u <- matrix(c(0.2, 0.3, 0.7, 0.8), ncol = 2)
  cop <- bicop_dist("gaussian", parameters = 0.4)

  expect_error(dbicop(u, cop, deriv = 1), "character vector")
  expect_error(dbicop(u, cop, deriv = character()), "length one or two")
  expect_error(
    dbicop(u, cop, deriv = c("u1", "u2", "par1")),
    "length one or two"
  )
  expect_error(dbicop(u, cop, deriv = NA_character_), "character vector")
  expect_error(dbicop(u, cop, deriv = "u3"), "components of 'deriv'")
  expect_error(dbicop(u, cop, deriv = "par0"), "components of 'deriv'")
  expect_error(dbicop(u, cop, deriv = "par2"), "parameter")
  expect_error(dbicop(u, cop, deriv = "u1", cores = 0), "not greater than 0")
  expect_error(
    hbicop(u, 1, cop, inverse = TRUE, deriv = "u1"),
    "inverse h-functions"
  )
  expect_error(dbicop(u, bicop_dist("tll"), deriv = "u1"), "not implemented")
  expect_error(
    dbicop(
      cbind(u, u[, 1] - 0.01),
      "gaussian",
      parameters = 0.4,
      var_types = c("d", "c"),
      deriv = "u1"
    ),
    "continuous variable types"
  )
  expect_error(
    hbicop(
      u,
      1,
      "gaussian",
      parameters = matrix(c(0.2, 0.3, 0.4), ncol = 1),
      deriv = "par1"
    ),
    "one row per row of u"
  )
  expect_error(
    rvinecopulib:::bicop_deriv_cpp(u, cop, "pdf", "u1", 3, 1),
    "order must be one or two"
  )
  expect_error(
    rvinecopulib:::bicop_deriv_cpp(u, cop, "unknown", "u1", 1, 1),
    "unknown derivative target"
  )
})
