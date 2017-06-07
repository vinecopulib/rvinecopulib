context("Fitting vine copula models")

set.seed(5)
u <- sapply(1:5, function(i) runif(30))
fit <- vinecop_select(u, "nonpar")
fit_no_data <- vinecop_select(u, "nonpar", keep_data = FALSE)

test_that("returns proper 'bicop' object", {
    expect_s3_class(fit, "vinecop")
    expect_s3_class(fit, "vinecop_dist")
    expect_identical(names(fit),  c("pair_copulas", "matrix", "data", "controls"))
    expect_identical(names(fit_no_data),  c("pair_copulas", "matrix", "controls"))
})

test_that("S3 generics work", {
    expect_equal(predict(fit, u), fitted(fit))
    expect_equal(
        predict(fit, u, what = "cdf"), 
        fitted(fit, what = "cdf"),
        tolerance = 0.01
    )
    expect_error(predict(fit, u, what = "hfunc1"))
    expect_length(attr(logLik(fit), "df"), 1)
    expect_error(logLik(fit_no_data))
})