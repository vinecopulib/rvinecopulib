context("Fitting 'vine' models")

set.seed(5)
u <- sapply(1:5, function(i) rnorm(30))
fit <- vine(u, copula_controls = list(family_set = "nonpar"), keep_data = TRUE)

test_that("returns proper 'vine' object", {
    expect_s3_class(fit, "vine")
    expect_s3_class(fit, "vine_dist")
    expect_identical(
        names(fit),  
        c("margins", "margins_controls", "copula", 
          "copula_controls", "npars", "loglik", "data", "nobs", "names")
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
    s <- summary(fit)
    expect_is(s$margins, c("summary_df", "data.frame"))
    expect_is(s$copula, c("summary_df", "data.frame"))
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

test_that("truncation works", {

    fit_truncated <- truncate_model(fit, trunc_lvl = 1)
    expect_silent(dvine(u, fit_truncated))
    expect_silent(rvine(50, fit_truncated))
    
    fit_truncated <- vine(u, copula_controls = list(family_set = "nonpar", 
                                                    trunc_lvl = 1))
    expect_silent(dvine(u, fit_truncated))
    expect_silent(rvine(50, fit_truncated))
})

test_that("margins_controls works", {
    fit_mult <- vine(u, margins_controls = list(mult = 2))
    expect_equal(sapply(fit_mult$margins, "[[", "bw"), 
                 2 / log(6) * sapply(fit$margins, "[[", "bw"))
    
    fit_xmin <- vine(abs(u), margins_controls = list(xmin = 0))
    expect_equal(sapply(fit_xmin$margins, "[[", "xmin"), rep(0, 5))
})
