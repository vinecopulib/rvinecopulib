context("Discrete variables")

test_that("bicop works", {
  u <- rbicop(1000, "gum", 90, 4)
  u_sub <- u
  u[, 1] <- ceiling(u[, 1] * 3) / 3
  u_sub[, 1] <- floor(u_sub[, 1] * 3) / 3
  cop <- bicop(u, family = "gumbel", data_sub = u_sub,
               var_types = c("d", "c"), presel = FALSE)
  expect_equal(cop$family, "gumbel")
  expect_equal(cop$rotation, 90)
  expect_equal(cop$parameters[1], 4, tol = 0.5)
  expect_identical(cop$var_types, c("d", "c"))

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
