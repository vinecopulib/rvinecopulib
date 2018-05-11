context("Fitting 'vine' models")

set.seed(5)
u <- sapply(1:5, function(i) rnorm(30))
fit <- vine(u, copula_controls = list(family_set = "nonpar"))

test_that("returns proper 'vine' object", {
    expect_s3_class(fit, "vine")
    expect_s3_class(fit, "vine_dist")
    expect_identical(
        names(fit),  
        c("margins", "margins_controls", "copula", 
          "copula_controls", "npars", "loglik", "data", "nobs")
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
})

test_that("print/summary generics work", {
    expect_output(print(fit))
    expect_s3_class(s <- summary(fit), "vinecop_dist_summary")
    expect_is(s, "data.frame")
})

test_that("discrete data work", {
    n <- 1e2
    x1 <- rnorm(n)
    x2 <- ordered(sample(5, n, TRUE), 1:5)
    x3 <- x1 + as.numeric(x2) + rnorm(n, sd = 0.5)
    
    my_data <- data.frame(x1, x2, x3)
    fit <- vine(my_data)
    sim <- rvine(n * 10, fit)
    expect_equal(sort(unique(sim[,2])), 1:5)
})

