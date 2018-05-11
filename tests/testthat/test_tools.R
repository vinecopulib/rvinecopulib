context("Exported statistical tools")

set.seed(0)

test_that("'pseudo_obs' works", {
    
    x <- rnorm(10)
    expect_equal(sort(pseudo_obs(x)), (1:10)/11)
    expect_equal(sort(pseudo_obs(x, lower_tail = FALSE)), (1:10)/11)
    expect_equal(pseudo_obs(x), 1 - pseudo_obs(x, lower_tail = FALSE))
    
    u <- pseudo_obs(cbind(rnorm(10), rnorm(10)))
    expect_equal(sort(u[, 1]), (1:10)/11)
    expect_equal(sort(u[, 2]), (1:10)/11)
    
    expect_error(pseudo_obs("something"))
    expect_error(pseudo_obs(x, "something"))
})
