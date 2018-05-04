context("Class 'vine_dist'")

set.seed(0)
bicop <- bicop_dist("bb1", 90, c(3, 2))
pcs <- list(list(bicop, bicop), list(bicop))
mat <- matrix(c(1, 2, 3, 1, 2, 0, 1, 0, 0), 3, 3)
vc <- vine_dist(list(name = "norm"), pcs, mat)

test_that("constructor creates proper `vine_dist` object", {
    expect_s3_class(vc, "vine_dist")
    expect_identical(names(vc), c("margins", "copula", "npars", "loglik"))
})


test_that("d/p/r- functions work", {
    u <- rvine(50, vc)
    u <- rvine(50, vc, pnorm(u))
    expect_gte(min(dvine(u, vc)), 0)
    expect_gte(min(pvine(u, vc, 100)), 0)
    expect_lte(max(pvine(u, vc, 100)), 1)
})


test_that("constructor catches wrong input", {
    # missing margin name
    expect_error(vine_dist(list(stupid = "norm"), cop))
    
    # unused margin argument
    expect_error(vine_dist(list(name = "norm", stupid = 42), cop))
    
    # missing margin argument
    expect_error(vine_dist(list(name = "beta", scale1 = 1), cop))
    
    # length of margins vector do not correspond to cop
    expect_error(vine_dist(list(list(name = "norm"), 
                                list(name = "gamma", shape = 1)), cop))

})

test_that("print/summary generics work", {
    expect_output(print(vc))
    expect_output(s <- summary(vc))
    expect_is(s, "data.frame")
    expect_equal(nrow(s), 3)
    expect_equal(ncol(s), 7)
})

test_that("getters work", {
    
    # test get_all_pair_copulas
    expect_identical(pcs, get_all_pair_copulas(vc))
    
    # test get_all_pair_copulas
    pcc <- get_pair_copula(vc, 1, 1)
    expect_equal(bicop, bicop_dist(pcc$family,pcc$rotation, pcc$parameters))
    
    # test get_matrix
    expect_identical(mat, get_matrix(vc))
    
    # test truncate 
    expect_identical(vc$copula$pair_copulas[1:1], truncate_model(vc, 1)[[2]][[1]])
})
