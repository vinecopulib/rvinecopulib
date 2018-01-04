context("Class 'bicop_dist'")

test_that("constructor creates proper bicop_dist object", {
    dist <- bicop_dist("gumbel", 90, 3)
    expect_s3_class(dist, "bicop_dist")
    expect_identical(names(dist), c("family", "rotation", "parameters", "npars"))
})


test_that("checks for family/rotation/parameters consistency", {
    expect_error(bicop_dist("asdf", 90, 3))
    expect_error(bicop_dist("frank", 3, 3))
    expect_error(bicop_dist("frank", 90, -3))
    expect_error(bicop_dist("frank", 90, 1:2))
})

test_that("partial matching for family names", {
    expect_equal(bicop_dist("ind")$family, "indep")
    expect_equal(bicop_dist("gauss")$family, "gaussian")
    expect_error(bicop_dist("g"))
})

test_that("d/p/r/h functions work", {
    set.seed(0)
    dist <- bicop_dist("bb1", 270, c(1, 2))
    u <- rbicop(50, "bb1", 270, c(1, 2))
    u <- rbicop(50, dist, U = u)
    expect_error(rbicop(50, dist, U = u[, -1]))
    expect_error(rbicop(50, dist, U = u[-1, ]))
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
    
    # one-parameter family
    tau <- par_to_tau(dist)
    expect_identical(tau, par_to_tau("joe", 90, 3))
    
    par <- tau_to_par(dist, tau)
    expect_identical(par, tau_to_par("joe", tau))
    
    expect_equal(3, par[1])
    
    # two-parameter
    tau <- par_to_tau("bb1", 0, c(1, 2))
    expect_error(tau_to_par("bb1", 0.5))
})

test_that("print method produces output", {
    dist <- bicop_dist("indep")
    expect_output(print(dist))
    expect_output(summary(dist))
})