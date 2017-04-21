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
    bicop_dist("ind")
    bicop_dist("gauss") 
})

test_that("d/r/hx/hix functions work", {
    dist <- bicop_dist("bb1", 270, c(1, 2))
    u <- rbicop(500, "bb1", 270, c(1, 2))
    expect_gte(min(dbicop(c(0.1, 0.2), dist)), 0)
    expect_gte(min(dbicop(u, dist)), 0)
    expect_lte(max(h1bicop(c(0.1, 0.2), dist)), 1)
    expect_lte(max(h2bicop(u, dist)), 1)
    expect_lte(max(hi1bicop(c(0.1, 0.2), dist)), 1)
    expect_lte(max(hi2bicop(u, dist)), 1)
})