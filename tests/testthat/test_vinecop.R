context("Fitting 'vinecop' models")

set.seed(5)
u <- sapply(1:5, function(i) runif(30))
fit <- vinecop(u, "nonpar")
fit_with_data <- vinecop(u, "nonpar", keep_data = TRUE)

test_that("returns proper 'vinecop' object", {
    expect_s3_class(fit, "vinecop")
    expect_s3_class(fit, "vinecop_dist")
    expect_identical(
        names(fit),  
        c("pair_copulas", "structure", "npars", "loglik", "threshold", "controls", "nobs")
    )
    expect_identical(
        names(fit_with_data), 
        c("pair_copulas", "structure", "npars", "loglik", "threshold", "data", "controls", "nobs")
    )
    
    colnames(u) <- paste(seq_len(ncol(u)))
    expect_identical(
        names(vinecop(u, "indep")), 
        c("pair_copulas", "structure", "npars", "loglik", "threshold", "names", "controls", "nobs")
    )
})

test_that("works with structure", {
    u <- sapply(1:2, function(i) runif(30))
    expect_silent(fit <- vinecop(u, structure = matrix(c(1:2, 1:0), 2, 2)))
})

if (Sys.info()['sysname'] != "SunOS") {
    test_that("runs in parallel", {
        expect_silent(fit <- vinecop(u, cores = 2))
    })   
}

test_that("S3 generics work", {
    expect_equal(predict(fit, u), fitted(fit_with_data))
    expect_equal(
        predict(fit, u, what = "cdf"), 
        fitted(fit_with_data, what = "cdf"),
        tolerance = 0.01
    )
    expect_error(predict(fit, u, what = "hfunc1"))
    expect_length(attr(logLik(fit), "df"), 1)
})

test_that("print/summary generics work", {
    expect_output(print(fit))
    expect_s3_class(s <- summary(fit), "summary_df")
    expect_is(s, "data.frame")
})

test_that("truncation works", {
    
    fit_truncated <- truncate_model(fit, trunc_lvl = 1)
    expect_silent(dvinecop(u, fit_truncated))
    expect_silent(rvinecop(50, fit_truncated))
    
    fit_truncated <- vinecop(u, "par", trunc_lvl = 1)
    expect_silent(dvinecop(u, fit_truncated))
    expect_silent(rvinecop(50, fit_truncated))
})
