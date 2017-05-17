context("Fitting 'bicop_fit' models")

dist <- bicop_dist("gumbel", 90, 3)
u <- rbicop(20, dist)

test_that("returns proper 'bicop_fit' object", {
    fit <- bicop_fit(u, "clayton")
    expect_s3_class(fit, "bicop_fit")
    expect_s3_class(fit, "bicop_dist")
    expect_identical(
        names(fit), 
        c("family", "rotation", "parameters", "npars", "data", "controls")
    )
    
    fit <- bicop_fit(u, "tll0", keep_data = FALSE)
    expect_s3_class(fit, "bicop_fit")
    expect_s3_class(fit, "bicop_dist")
    expect_identical(
        names(fit), 
        c("family", "rotation", "parameters", "npars", "controls")
    )
})

test_that("partial matching for family set names", {
    bicop_fit(u, "arch")
    bicop_fit(u, "nonp")
    expect_error(bicop_fit(u, "asdf"))
})
