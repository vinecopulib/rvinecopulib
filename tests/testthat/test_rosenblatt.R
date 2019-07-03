context("Inverse Rosenblatt transform")

pc <- bicop_dist("bb1", 90, c(3, 2))
pcs <- list(list(pc, pc), list(pc))
mat <- matrix(c(1, 2, 3, 1, 2, 0, 1, 0, 0), 3, 3)
vc <- vinecop_dist(pcs, mat)
vd <- vine_dist(list(distr = "norm"), pcs, mat)

test_that("rosenblatt works with bivariate copulas", {
  u <- rbicop(20, pc)
  expect_equal(inverse_rosenblatt(rosenblatt(u, pc), pc), u)
  pc <- bicop(u, family = "clay")
  expect_equal(inverse_rosenblatt(rosenblatt(u, pc), pc), u)
})

test_that("rosenblatt works with vine copulas", {
  u <- rvinecop(20, vc)
  expect_equal(inverse_rosenblatt(rosenblatt(u, vc), vc), u)
  vc <- vinecop(u, structure = mat, family = "clay")
  expect_equal(inverse_rosenblatt(rosenblatt(u, vc), vc), u)
})

test_that("rosenblatt works with vine distribution", {
  u <- rvine(20, vd)
  expect_equal(inverse_rosenblatt(rosenblatt(u, vd), vd), u)
  vd <- vine(u, copula_controls = list(structure = mat, family = "clay"))
  expect_equal(inverse_rosenblatt(rosenblatt(u, vd), vd), u)
})
