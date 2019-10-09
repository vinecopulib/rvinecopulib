context("Discrete variables")

set.seed(5)

test_that("bicop_dist works", {
  cop <- bicop_dist("gum", 90, 4, c("d", "c"))
  expect_identical(cop$var_types, c("d", "c"))

  u <- rbicop(1000, "gum", 90, 4)
  u_sub <- u
  u[, 1] <- ceiling(u[, 1] * 3) / 3
  u_sub[, 1] <- floor(u_sub[, 1] * 3) / 3

  # only check for errors
  dbicop(u, cop, u_sub = u_sub)
  pbicop(u, cop, u_sub = u_sub)
  hbicop(u, 1, cop, u_sub = u_sub)
  hbicop(u, 2, cop, TRUE, u_sub = u_sub)
  rbicop(10, cop)

  dbicop(u, "gum", 90, 4, c("d", "c"), u_sub)
  pbicop(u, "gum", 90, 4, c("d", "c"), u_sub)
  hbicop(u, 2, "gum", 90, 4, FALSE, c("d", "c"), u_sub)
  hbicop(u, 1, "gum", 90, 4, TRUE, c("d", "c"), u_sub)
})


test_that("bicop works", {
  u <- rbicop(1000, "gum", 90, 4)
  u_sub <- u
  u[, 1] <- ceiling(u[, 1] * 3) / 3
  u_sub[, 1] <- floor(u_sub[, 1] * 3) / 3
  cop <- bicop(u, family = "gumbel", data_sub = u_sub, var_types = c("d", "c"),
               presel = FALSE, keep_data = TRUE)
  expect_equal(cop$family, "gumbel")
  expect_equal(cop$rotation, 90)
  expect_equal(cop$parameters[1], 4, tol = 0.5)
  expect_identical(cop$var_types, c("d", "c"))
  expect_identical(cop$data_sub, u_sub)#
})

# -----------------------------------------------------------------

set.seed(0)
bicop <- bicop_dist("bb1", 90, c(3, 2))
pcs <- list(list(bicop, bicop), list(bicop))
mat <- matrix(c(1, 2, 3, 1, 2, 0, 1, 0, 0), 3, 3)
var_types <- c("d", "d", "c")

test_that("vinecop_dist works", {
  cop <- vinecop_dist(pcs, mat, var_types)
  expect_equal(cop$var_types, var_types)

  u <- rvinecop(20, cop)
  u_sub <- u
  u[, 1] <- ceiling(u[, 1] * 3) / 3
  u_sub[, 1] <- floor(u_sub[, 1] * 3) / 3
  u[, 2] <- ceiling(u[, 2] * 10) / 10
  u_sub[, 2] <- floor(u_sub[, 2] * 10) / 10

  # only check for errors
  dvinecop(u, cop, u_sub = u_sub)
  pvinecop(u, cop, u_sub = u_sub)
})

test_that("vinecop works", {
  u <- replicate(3, runif(200))
  u_sub <- u
  u[, 1] <- ceiling(u[, 1] * 3) / 3
  u_sub[, 1] <- floor(u_sub[, 1] * 3) / 3
  u[, 2] <- ceiling(u[, 2] * 10) / 10
  u_sub[, 2] <- floor(u_sub[, 2] * 10) / 10
  cop <- vinecop(u, data_sub = u_sub, var_types = var_types, family = "tll")
  summary(cop)
  expect_identical(cop$var_types, var_types)
})
