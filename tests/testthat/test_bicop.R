context("Fitting 'bicop' models")

dist <- bicop_dist("gumbel", 90, 3)
u <- rbicop(20, dist)

test_that("returns proper 'bicop' object", {
    fit <- bicop(u, "clayton")
    expect_s3_class(fit, "bicop")
    expect_s3_class(fit, "bicop_dist")
    expect_identical(
        names(fit), 
        c("family", "rotation", "parameters", "npars", "data", "controls")
    )
    
    fit <- bicop(u, "tll0", keep_data = FALSE)
    expect_s3_class(fit, "bicop")
    expect_s3_class(fit, "bicop_dist")
    expect_identical(
        names(fit), 
        c("family", "rotation", "parameters", "npars", "controls")
    )
})

test_that("partial matching for family set names", {
    bicop(u, "arch")
    bicop(u, "nonp")
    expect_error(bicop(u, "asdf"))
})
