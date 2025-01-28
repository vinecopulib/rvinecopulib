# fixes problems with change in all.equal() behavior in R 4.1.x
expect_eql <- function(...) expect_equal(..., check.environment = FALSE)
expect_equiv <- function(...) expect_equivalent(..., check.environment = FALSE)

context("(Inverse) Rosenblatt transform")

pc <- bicop_dist("bb1", 90, c(3, 2))
pcs <- list(list(pc, pc), list(pc))
mat <- matrix(c(1, 2, 3, 1, 2, 0, 1, 0, 0), 3, 3)
vc <- vinecop_dist(pcs, mat)
vd <- vine_dist(list(list(distr = "norm")), pcs, mat)

test_that("rosenblatt works with bivariate copulas", {
  u <- rbicop(20, pc)
  expect_eql(inverse_rosenblatt(rosenblatt(u, pc), pc), u)
  pc <- bicop(u, family = "clay")
  expect_eql(inverse_rosenblatt(rosenblatt(u, pc), pc), u)
})

test_that("rosenblatt works with vine copulas", {
  u <- rvinecop(20, vc)
  expect_eql(inverse_rosenblatt(rosenblatt(u, vc), vc), u)
  vc <- vinecop(u, structure = mat, family = "clay")
  expect_eql(inverse_rosenblatt(rosenblatt(u, vc), vc), u)
})

test_that("discrete rosenblatt works with vine copulas", {
  u <- rvinecop(2000, vc)
  uu <- cbind(u, u[, 2])

  thresh <- 0.05
  uu[u[, 2] < thresh, 2] <- 1e-10
  uu[u[, 2] < thresh, 4] <- thresh

  vc_c <- vc <- vinecop(uu, var_types = c("c", "d", "c"),
                        structure = mat, family = "clay")
  vc_c$var_types = rep("c", 3)

  v <- inverse_rosenblatt(rosenblatt(uu, vc), vc_c)
  expect_eql(v, u, tol = 2 * thresh)

  # other format
  uu <- cbind(u, u)
  uu[u[, 2] < thresh, 2] <- 1e-10
  uu[u[, 2] < thresh, 5] <- thresh
  v <- inverse_rosenblatt(rosenblatt(uu, vc), vc_c)
  expect_eql(v, u, tol = 2 * thresh)
})

test_that("rosenblatt works with vine distribution", {
  u <- rvine(20, vd)
  expect_eql(inverse_rosenblatt(rosenblatt(u, vd), vd), u)
  vd <- vine(u, copula_controls = list(structure = mat, family = "clay"))
  expect_equiv(inverse_rosenblatt(rosenblatt(u, vd), vd), u)

})

test_that("discrete rosenblatt works with vine distributions", {
  x <- data.frame(
    x1 = as.ordered(sample(1:4, 50, replace = TRUE)),
    x2 = rnorm(50),
    x3 = rbinom(50, 3, 0.5)
  )
  vd <- vine(x, margins_controls = list(type = c("d", "c", "zi")),
             copula_controls = list(family_set = "gauss"))
  u <- rvine(50, vd)
  expect_true(mean(u[, 3] == 0) > 0.02)
  expect_eql(colnames(rosenblatt(x, vd)), colnames(x))
})
