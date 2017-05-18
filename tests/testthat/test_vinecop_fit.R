context("Fitting vine copula models")

u <- sapply(1:5, function(i) runif(30))

test_that("returns proper 'bicop_fit' object", {
    fit <- vinecop_select(u, "nonpar")
    expect_s3_class(fit, "vinecop_fit")
    expect_s3_class(fit, "vinecop_dist")
    expect_identical(names(fit),  c("pair_copulas", "matrix", "data", "controls"))
    
    fit <- vinecop_select(u, "nonpar", keep_data = FALSE)
    expect_identical(names(fit),  c("pair_copulas", "matrix", "controls"))
})
