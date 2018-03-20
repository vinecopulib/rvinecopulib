context("Fitting `vine` models")

set.seed(5)
u <- sapply(1:5, function(i) rnorm(30))
fit <- vine(u, copula_controls = list(family_set = "nonpar"))
fit_no_data <- vine(u, copula_controls = list(family_set = "nonpar", 
                                              keep_data = FALSE))

test_that("returns proper 'vine' object", {
    expect_s3_class(fit, "vine")
    expect_s3_class(fit, "vine_dist")
    expect_identical(
        names(fit),  
        c("margins", "margins_controls", "copula", 
          "copula_controls", "data", "nobs")
    )
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

test_that("print/summary generics work", {
    expect_output(print(fit))
    expect_output(s <- summary(fit))
    expect_is(s, "data.frame")
    expect_length(attr(s, "info"), 5)
})