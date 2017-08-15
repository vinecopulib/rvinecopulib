context("Class 'vinecop_dist'")

bicop <- bicop_dist("bb1", 90, c(3, 2))
pcs <- list(list(bicop, bicop), list(bicop))
mat <- matrix(c(1, 2, 3, 1, 2, 0, 1, 0, 0), 3, 3)
vc <- vinecop_dist(pcs, mat)


test_that("constructor creates proper vinecop_dist object", {
    expect_s3_class(vc, "vinecop_dist")
    expect_identical(names(vc), c("pair_copulas", "matrix"))
})


test_that("d/p/r- functions work", {
    u <- rvinecop(50, vc)
    expect_gte(min(dvinecop(u, vc)), 0)
    expect_gte(min(pvinecop(u, vc, 100)), 0)
    expect_lte(max(pvinecop(u, vc, 100)), 1)
})

test_that("constructor catches wrong R-vine matrix", {
    mat[3, 3] <- 5
    expect_error(vinecop_dist(pcs, mat))
})

test_that("print/summary generics work", {
    expect_output(print(vc))
    expect_output(s <- summary(vc))
    expect_is(s, "data.frame")
    expect_equal(nrow(s), 3)
    expect_equal(ncol(s), 7)
})