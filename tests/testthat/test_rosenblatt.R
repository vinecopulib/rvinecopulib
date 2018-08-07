context("Inverse Rosenblatt transform")

pc <- bicop_dist("bb1", 90, c(3, 2))
pcs <- list(list(pc, pc), list(pc))
mat <- matrix(c(1, 2, 3, 1, 2, 0, 1, 0, 0), 3, 3) 
vc <- vinecop_dist(pcs, mat)
vd <- vine_dist(list(distr = "norm"), pcs, mat)

test_that("inverse rosenblatt works with bivariate copulas", {
    u <- replicate(2, runif(50))
    expect_equal(
        inverse_rosenblatt(u, pc), 
        cbind(u[, 1], hbicop(u, 1, pc, inverse = TRUE))
    )
    pc <- bicop(u, family = "clay")
    expect_equal(
        inverse_rosenblatt(u, pc), 
        cbind(u[, 1], hbicop(u, 1, pc, inverse = TRUE))
    )
})

test_that("inverse rosenblatt works with vine copulas", {
    u <- replicate(3, runif(50))
    expect_equal(inverse_rosenblatt(u, vc)[, 1], u[, 1])
    vc <- vinecop(u, structure = mat, family = "clay")
    expect_equal(inverse_rosenblatt(u, vc)[, 1], u[, 1])
})

test_that("inverse rosenblatt works with vine distribution", {
    u <- replicate(3, runif(50))
    expect_equal(inverse_rosenblatt(u, vd)[, 1], qnorm(u[, 1]))
    vd <- vine(u, copula_controls = list(structure = mat, family = "clay"))
    expect_equal(
        inverse_rosenblatt(u, vd)[, 1], 
        kde1d::qkde1d(u[, 1], vd$margins[[1]])
    )
})
