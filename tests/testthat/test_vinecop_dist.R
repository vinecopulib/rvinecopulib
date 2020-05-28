context("Class 'vinecop_dist'")

set.seed(0)
bicop <- bicop_dist("bb1", 90, c(3, 2))
pcs <- list(list(bicop, bicop), list(bicop))
mat <- matrix(c(1, 2, 3, 1, 2, 0, 1, 0, 0), 3, 3)
vc <- vinecop_dist(pcs, mat)


test_that("constructor creates proper vinecop_dist object", {
  expect_s3_class(vc, "vinecop_dist")
  expect_identical(
    names(vc),
    c("pair_copulas", "structure", "var_types", "npars", "loglik")
  )
})


test_that("d/p/r- functions work", {
  u <- rvinecop(50, vc)
  expect_false(any(rvinecop(50, vc, qrng = FALSE) ==
    rvinecop(50, vc, qrng = FALSE)))
  set.seed(1)
  u <- rvinecop(50, vc, qrng = TRUE)
  set.seed(1)
  expect_true(all(u == rvinecop(50, vc, qrng = TRUE)))
  expect_gte(min(dvinecop(u, vc)), 0)
  expect_gte(min(pvinecop(u, vc, 100)), 0)
  expect_lte(max(pvinecop(u, vc, 100)), 1)
})


test_that("constructor catches wrong input", {
  # wrong number of pcs
  pcs2 <- pcs
  pcs2[[1]][[2]] <- NULL
  expect_error(vinecop_dist(pcs[-1], mat))

  # not all pcs are of class 'bicop_dist'
  pcs2[[1]][[2]] <- list(this = "stupid")
  expect_error(vinecop_dist(pcs2, mat))

  # wrong R-vine matrix
  mat[3, 3] <- 5
  expect_error(vinecop_dist(pcs, mat))
})

test_that("works with truncated vines", {
  # takes and returns truncated pair_copulas list
  trunc_vine <- vinecop_dist(pcs[-2], mat)
  expect_length(trunc_vine$pair_copulas, 1)

  # summary table is truncated too
  expect_s3_class(summary(vinecop_dist(pcs[-2], mat)), "summary_df")
  expect_silent(smr <- summary(vinecop_dist(pcs[-2], mat)))
  expect_equal(nrow(smr), 2)
})

test_that("print/summary/dim generics work", {
  expect_output(print(vc))
  expect_s3_class(summary(vc), "summary_df")
  expect_silent(s <- summary(vc))
  expect_is(s, "data.frame")
  expect_equal(nrow(s), 3)
  expect_equal(ncol(s), 10)
  expect_equivalent(dim(vc), c(3, 2))
})

test_that("plot functions work", {
  pcs <- lapply(1:4, function(j) # pair-copulas in tree j
    lapply(runif(5 - j), function(cor) bicop_dist("gaussian", 0, cor)))
  mat <- matrix(
    c(
      1, 2, 3, 4, 5,
      1, 2, 3, 4, 0,
      1, 2, 3, 0, 0,
      1, 2, 0, 0, 0,
      1, 0, 0, 0, 0
    ),
    5, 5
  )
  vc <- vinecop_dist(pcs, mat)

  # we could check some values in the plot objects
  expect_silent(p <- plot(vc, edge_labels = "family", var_names = "legend"))
  expect_silent(p <- plot(vc, edge_labels = "tau", var_names = "use"))
  expect_silent(p <- plot(vc, edge_labels = "pair"))
  expect_silent(p <- plot(vc, edge_labels = "family_tau"))
  expect_silent(p <- plot(vc, var_names = "hide"))
  expect_error(p <- plot(vc, edge_labels = "no"))
  expect_error(p <- plot(vc, var_names = "isaidno"))
  expect_error(p <- plot(vc, tree = 10))
  expect_silent(p <- plot(vc, "ALL"))
  expect_silent(p <- contour(vc, xlim = c(0.2, 0.8), ylim = c(0.2, 0.8)))
  expect_silent(p <- contour(vc, margins = "unif"))
  expect_error(p <- contour(vc, margins = "nonono"))
  expect_error(p <- contour(vc, var_names = "comeon"))

  # contour for truncated vines
  vc$pair_copulas[[4]] <- NULL
  expect_silent(p <- contour(vc, margins = "unif"))
})

test_that("getters work", {

  # test get_structure
  expect_silent(pcc <- get_structure(vc))
  expect_error(get_structure(12))

  # test get_matrix
  expect_equivalent(as_rvine_matrix(mat), get_matrix(vc))
  expect_error(get_matrix(12))

  # test get_pair_copulas
  expect_silent(pcc <- get_pair_copula(vc, 1, 1))
  expect_equal(bicop, bicop_dist(pcc$family, pcc$rotation, pcc$parameters))
  expect_error(get_pair_copula(12, 1, 1))
  expect_error(get_pair_copula(vc, 1:2, 1))
  expect_error(get_pair_copula(vc, 1, 1:2))
  expect_error(get_pair_copula(vc, 0, 1))
  expect_error(get_pair_copula(vc, 1, 0))
  expect_error(get_pair_copula(vc, 12, 1))
  expect_error(get_pair_copula(vc, 1, 12))

  # test get_all_pair_copulas
  expect_equivalent(pcs, get_all_pair_copulas(vc))
  expect_equivalent(pcs[1:2], get_all_pair_copulas(vc, 1:2))
  expect_error(get_all_pair_copulas(12))
  expect_error(get_all_pair_copulas(vc, 0))
  expect_error(get_all_pair_copulas(vc, 12))

  # test other getters
  expect_equivalent(get_parameters(vc, 1, 1), coef(pcs[[1]][[1]]))
  expect_equivalent(
    get_all_parameters(vc),
    lapply(pcs, function(tree) lapply(tree, coef))
  )
  expect_equivalent(get_ktau(vc, 1, 1), par_to_ktau(bicop))
  expect_equivalent(
    get_all_ktaus(vc),
    lapply(pcs, function(tree)
      lapply(tree, function(pc) par_to_ktau(pc)))
  )
  expect_equivalent(get_family(vc, 1, 1), "bb1")
  expect_equivalent(
    get_all_families(vc),
    lapply(pcs, function(tree)
      lapply(tree, function(pc) pc$family))
  )

  # test printed output of getters
  expect_output(print(get_all_pair_copulas(vc)))
  expect_output(print(get_all_pair_copulas(vc, 1)))
})


test_that("d = 1 works", {
  vc <- vinecop_dist(list(), rvine_structure(1))
  u <- runif(5)
  expect_equal(dim(summary(vc))[1], 0)

  expect_equal(dvinecop(u, vc), rep(1, 5))
  expect_equal(pvinecop(u, vc), u, tol = 1e-2)
  expect_equal(c(rosenblatt(u, vc)), u)
  expect_equal(c(inverse_rosenblatt(u, vc)), u)
  expect_silent(rvinecop(10, vc))
  expect_error(plot(vc))
  expect_warning(contour(vc))
})
