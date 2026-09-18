# fixes problems with change in all.equal() behavior in R 4.1.x
expect_eql <- function(...) expect_equal(..., check.environment = FALSE)
expect_equiv <- function(...) expect_equivalent(..., check.environment = FALSE)

expect_pdf_full_triangular_vectors <- function(x, vc, n) {
  arrays <- c("pdf_edges", "hfunc1", "hfunc2", "hfunc1_sub", "hfunc2_sub")
  sub_arrays <- c("hfunc1_sub", "hfunc2_sub")
  expect_type(x, "list")
  expect_type(x$pdf, "double")
  expect_null(dim(x$pdf))
  expect_length(x$pdf, n)
  for (array in arrays) {
    expect_type(x[[array]], "list")
    expected_length <- if (array %in% sub_arrays && all(vc$var_types == "c")) {
      0L
    } else {
      dim(vc)["trunc_lvl"]
    }
    expect_length(x[[array]], expected_length)
    for (tree in seq_along(x[[array]])) {
      expect_type(x[[array]][[tree]], "list")
      expect_length(x[[array]][[tree]], dim(vc)[1] - tree)
      for (edge in seq_along(x[[array]][[tree]])) {
        expect_type(x[[array]][[tree]][[edge]], "double")
        expect_null(dim(x[[array]][[tree]][[edge]]))
        expect_true(length(x[[array]][[tree]][[edge]]) %in% c(0L, n))
      }
    }
  }
}

context("Discrete variables")

set.seed(5)

test_that("bicop_dist works", {
  cop <- bicop_dist("gum", 90, 4, c("d", "c"))
  expect_identical(cop$var_types, c("d", "c"))

  u <- rbicop(10, "gum", 90, 4)
  u <- cbind(u, u)
  u[, 1] <- ceiling(u[, 1] * 3) / 3
  u[, 3] <- floor(u[, 3] * 3) / 3

  # only check for errors
  dbicop(u, cop)
  pbicop(u, cop)
  hbicop(u, 1, cop)
  hbicop(u, 2, cop, TRUE)
  rbicop(10, cop)

  dbicop(u, "gum", 90, 4, c("d", "c"))
  pbicop(u, "gum", 90, 4, c("d", "c"))
  hbicop(u, 2, "gum", 90, 4, FALSE, c("d", "c"))
  hbicop(u, 1, "gum", 90, 4, TRUE, c("d", "c"))
})


test_that("bicop works", {
  u <- rbicop(1000, "gum", 90, 4)
  u <- cbind(u, u)
  u[, 1] <- ceiling(u[, 1] * 3) / 3
  u[, 3] <- floor(u[, 3] * 3) / 3
  cop <- bicop(
    u,
    family = "gumbel",
    var_types = c("d", "c"),
    presel = FALSE,
    keep_data = TRUE
  )
  expect_eql(cop$family, "gumbel")
  expect_eql(cop$rotation, 90)
  expect_eql(cop$parameters[1], 4, tol = 0.5)
  expect_identical(cop$var_types, c("d", "c"))
})

test_that("discrete TLL fitting uses interval left limits", {
  set.seed(42)
  z1 <- rnorm(120)
  z2 <- 0.7 * z1 + sqrt(1 - 0.7^2) * rnorm(120)
  x1 <- qpois(pnorm(z1), 2)
  x2 <- qpois(pnorm(z2), 3)
  u <- cbind(
    ppois(x1, 2),
    ppois(x2, 3),
    ppois(x1 - 1, 2),
    ppois(x2 - 1, 3)
  )
  u_narrow <- u
  u_narrow[, 3:4] <- 0.5 * (u_narrow[, 1:2] + u_narrow[, 3:4])

  fit <- bicop(u, family_set = "tll", var_types = c("d", "d"))
  fit_narrow <- bicop(
    u_narrow,
    family_set = "tll",
    var_types = c("d", "d")
  )

  expect_gt(max(abs(fit$parameters - fit_narrow$parameters)), 1e-3)
})

test_that("discrete density collapses narrow intervals at their midpoint", {
  cop_dd <- bicop_dist(
    "gaussian",
    parameters = 0.7,
    var_types = c("d", "d")
  )
  cop_cd <- bicop_dist(
    "gaussian",
    parameters = 0.7,
    var_types = c("c", "d")
  )
  cop_dc <- bicop_dist(
    "gaussian",
    parameters = 0.7,
    var_types = c("d", "c")
  )

  narrow_first <- matrix(c(0.400001, 0.7, 0.4, 0.2), nrow = 1)
  expected_first <- matrix(c(0.4000005, 0.7, 0.2), nrow = 1)
  expect_equal(
    dbicop(narrow_first, cop_dd),
    dbicop(expected_first, cop_cd),
    tolerance = 1e-12
  )

  narrow_second <- matrix(c(0.8, 0.300001, 0.1, 0.3), nrow = 1)
  expected_second <- matrix(c(0.8, 0.3000005, 0.1), nrow = 1)
  expect_equal(
    dbicop(narrow_second, cop_dd),
    dbicop(expected_second, cop_dc),
    tolerance = 1e-12
  )
})

# -----------------------------------------------------------------

set.seed(0)
bicop <- bicop_dist("bb1", 90, c(3, 2))
pcs <- list(list(bicop, bicop), list(bicop))
mat <- matrix(c(1, 2, 3, 1, 2, 0, 1, 0, 0), 3, 3)
var_types <- c("d", "d", "c")

test_that("vinecop_dist works", {
  cop <- vinecop_dist(pcs, mat, var_types)
  expect_eql(cop$var_types, var_types)

  u <- rvinecop(20, cop)
  u <- cbind(u, u)
  u[, 1] <- ceiling(u[, 1] * 3) / 3
  u[, 3] <- floor(u[, 3] * 3) / 3
  u[, 2] <- ceiling(u[, 2] * 10) / 10
  u[, 3] <- floor(u[, 4] * 10) / 10

  # only check for errors
  dvinecop(u, cop)
  pdf_full <- dvinecop(u, cop, keep_all = TRUE)
  expect_pdf_full_triangular_vectors(pdf_full, cop, nrow(u))
  expect_eql(pdf_full$pdf, dvinecop(u, cop))
  pvinecop(u, cop)
})

test_that("vinecop works", {
  u <- replicate(3, runif(20))
  u <- cbind(u, u)
  u[, 1] <- ceiling(u[, 1] * 3) / 3
  u[, 2] <- ceiling(u[, 2] * 10) / 10
  u[, 4] <- floor(u[, 4] * 3) / 3
  u[, 5] <- floor(u[, 5] * 10) / 10
  colnames(u) <- paste("x", 1:3) |> rep(2)
  cop <- vinecop(u, var_types = var_types, family = "tll")
  summary(cop)
  expect_identical(cop$var_types, var_types)
})


# -----------------------------------------------------------------

test_that("vine works", {
  n <- 20
  x1 <- rnorm(n)
  x2 <- ordered(sample(5, n, TRUE), 1:5)
  x3 <- x1 + as.numeric(x2) + rnorm(n, sd = 0.5)
  x <- data.frame(x1, x2, x3)

  fit <- vine(x)
  sim <- rvine(n * 10, fit)
  expect_true(is.data.frame(sim))
  expect_eql(sort(unique(sim[, 2])), ordered(1:5))

  summary(fit)
  expect_identical(fit$copula$var_types, c("c", "d", "c"))
  # only check for errors
  dvine(x, fit)
  pvine(x, fit)
  rvine(10, fit)
})

test_that("rvine computes discrete conditioning limits internally", {
  n <- 30
  x <- data.frame(
    x1 = rnorm(n),
    x2 = ordered(sample(5, n, TRUE), 1:5),
    x3 = rnorm(n)
  )
  fit <- vine(
    x,
    copula_controls = list(
      family_set = "indep",
      conditioning_set = c("x3", "x2")
    )
  )
  x_cond <- data.frame(
    x3 = 1.25,
    x2 = ordered(3, levels = 1:5)
  )

  sim <- rvine(
    20,
    fit,
    x_cond = x_cond,
    conditioning_set = c("x3", "x2")
  )
  expect_equal(sim$x3, rep(1.25, 20))
  expect_identical(sim$x2, ordered(rep(3, 20), levels = 1:5))

  sim_factor <- rvine(
    5,
    fit,
    x_cond = ordered(3, levels = 1:5),
    conditioning_set = "x2"
  )
  expect_identical(sim_factor$x2, ordered(rep(3, 5), levels = 1:5))
})
