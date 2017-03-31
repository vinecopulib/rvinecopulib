context("Class 'bicop_dist'")

test_that("constructor creates proper bicop_dist object", {
    dist <- bicop_dist("gumbel", 90, 3)
    expect_s3_class(dist, "bicop_dist")
    expect_identical(names(dist), c("family", "rotation", "parameters"))
})


test_that("checks for family/rotation/parameters consistency", {
    expect_error(bicop_dist("asdf", 90, 3))
    expect_error(bicop_dist("frank", 3, 3))
    expect_error(bicop_dist("frank", 90, -3))
    expect_error(bicop_dist("frank", 90, 1:2))
})
