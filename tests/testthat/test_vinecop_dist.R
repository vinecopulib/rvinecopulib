context("Class 'vinecop_dist'")

bicop <- bicop_dist("bb1", 90, c(3, 2))
pcs <- list(list(bicop, bicop), list(bicop))
mat <- matrix(c(1, 2, 3, 1, 2, 0, 1, 0, 0), 3, 3)
vc <- vinecop_dist(pcs, mat)


test_that("constructor creates proper vinecop_dist object", {
    expect_s3_class(vc, "vinecop_dist")
    expect_identical(names(vc), c("pair_copulas", "matrix", "npars"))
})


test_that("d/p/r- functions work", {
    u <- rvinecop(50, vc)
    u <- rvinecop(50, vc, u)
    expect_gte(min(dvinecop(u, vc)), 0)
    expect_gte(min(pvinecop(u, vc, 100)), 0)
    expect_lte(max(pvinecop(u, vc, 100)), 1)
})

test_that("constructor catches wrong R-vine matrix", {
    mat[3, 3] <- 5
    expect_error(vinecop_dist(pcs, mat))
})

test_that("print/summary generics work", {
    expect_output(print(vc))
    expect_output(s <- summary(vc))
    expect_is(s, "data.frame")
    expect_equal(nrow(s), 3)
    expect_equal(ncol(s), 7)
})

test_that("plot functions work", {
    pcs <- lapply(1:4, function(j) # pair-copulas in tree j
        lapply(runif(5-j), function(cor) bicop_dist("gaussian", 0, cor)))
    mat <- matrix(c(1, 2, 3, 4, 5, 
                    1, 2, 3, 4, 0, 
                    1, 2, 3, 0, 0, 
                    1, 2, 0, 0, 0,
                    1, 0, 0, 0, 0), 
                  5, 5)
    vc <- vinecop_dist(pcs, mat)
    
    # we could check some values in the plot objects
    expect_silent(p <- plot(vc, edge_labels = "family", var_names = "legend"))
    expect_silent(p <- plot(vc, edge_labels = "tau", var_names = "use"))
    expect_silent(p <- plot(vc, edge_labels = "pair"))
    expect_silent(p <- plot(vc, edge_labels = "family_tau"))
    expect_error(p <- plot(vc, edge_labels = "no"))
    expect_error(p <- plot(vc, var_names = "isaidno"))
    expect_error(p <- plot(vc, tree = 10))
    expect_silent(p <- plot(vc, "ALL"))
    expect_silent(p <- contour(vc, xlim = c(0.2, 0.8), ylim = c(0.2, 0.8)))
    expect_silent(p <- contour(vc, margins = "unif"))
    expect_error(p <- contour(vc, margins = "nonono"))
    expect_error(p <- contour(vc, var_names = "comeon"))
})