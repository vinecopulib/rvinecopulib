# fixes problems with change in all.equal() behavior in R 4.1.x
expect_eql <- function(...) expect_equal(..., check.environment = FALSE)
expect_equiv <- function(...) expect_equivalent(..., check.environment = FALSE)

expect_pdf_full_triangular_vectors <- function(x, vc, n) {
  arrays <- c("pdf_edges", "hfunc1", "hfunc2", "hfunc1_sub", "hfunc2_sub")
  sub_arrays <- c("hfunc1_sub", "hfunc2_sub")
  expect_type(x, "list")
  expect_type(x$pdf, "double")
  expect_null(dim(x$pdf))
  expect_length(x$pdf, n)
  for (array in arrays) {
    expect_type(x[[array]], "list")
    expected_length <- if (
      array %in% sub_arrays && all(vc$var_types == "c")
    ) {
      0L
    } else {
      dim(vc)["trunc_lvl"]
    }
    expect_length(x[[array]], expected_length)
    for (tree in seq_along(x[[array]])) {
      expect_type(x[[array]][[tree]], "list")
      expect_length(x[[array]][[tree]], dim(vc)[1] - tree)
      for (edge in seq_along(x[[array]][[tree]])) {
        expect_type(x[[array]][[tree]][[edge]], "double")
        expect_null(dim(x[[array]][[tree]][[edge]]))
        expect_true(length(x[[array]][[tree]][[edge]]) %in% c(0L, n))
      }
    }
  }
}

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
  expect_false(any(
    rvinecop(50, vc, qrng = FALSE) == rvinecop(50, vc, qrng = FALSE)
  ))
  set.seed(1)
  u <- rvinecop(50, vc, qrng = TRUE)
  set.seed(1)
  expect_true(all(u == rvinecop(50, vc, qrng = TRUE)))
  expect_gte(min(dvinecop(u, vc)), 0)
  pdf_full <- dvinecop(u, vc, keep_all = TRUE)
  expect_named(
    pdf_full,
    c("pdf", "pdf_edges", "hfunc1", "hfunc2", "hfunc1_sub", "hfunc2_sub")
  )
  expect_eql(pdf_full$pdf, dvinecop(u, vc))
  expect_length(pdf_full$pdf_edges, dim(vc)["trunc_lvl"])
  expect_length(pdf_full$pdf_edges[[1]], dim(vc)[1] - 1)
  expect_pdf_full_triangular_vectors(pdf_full, vc, nrow(u))
  expect_equal(
    pdf_full$pdf_edges[[1]][[1]] *
      pdf_full$pdf_edges[[1]][[2]] *
      pdf_full$pdf_edges[[2]][[1]],
    pdf_full$pdf
  )
  expect_equal(dvinecop(u, vc, cores = 2), dvinecop(u, vc))
  expect_equal(dvinecop(u, vc, cores = 2, keep_all = TRUE), pdf_full)
  expect_gte(min(pvinecop(u, vc, 100)), 0)
  expect_lte(max(pvinecop(u, vc, 100)), 1)
})

test_that("scores and hessian work", {
  u <- rvinecop(10, vc)
  expect_silent(s <- scores(u, vc))
  expect_eql(dim(s), c(10L, as.integer(vc$npars)))
  expect_type(s, "double")
  expect_silent(s_full <- scores(u, vc, step_wise = FALSE))
  expect_eql(dim(s_full), dim(s))
  expect_true(all(is.finite(s)))
  expect_true(all(is.finite(s_full)))
  expect_silent(h <- hessian(u, vc))
  expect_eql(dim(h), c(as.integer(vc$npars), as.integer(vc$npars)))
  expect_type(h, "double")
  expect_silent(h_full <- hessian(u, vc, step_wise = FALSE))
  expect_eql(dim(h_full), dim(h))
  expect_true(all(is.finite(h)))
  expect_true(all(is.finite(h_full)))
  expect_equal(scores(u, vc, cores = 2), s)
  expect_equal(scores(u, vc, step_wise = FALSE, cores = 2), s_full)
  expect_equal(hessian(u, vc, cores = 2), h)
  expect_equal(hessian(u, vc, step_wise = FALSE, cores = 2), h_full)
  expect_error(scores(u, vc, step_wise = 1))
  expect_error(scores(u, vc, cores = 0))
  expect_error(hessian(u, vc, step_wise = 1))
  expect_error(hessian(u, vc, cores = 0))
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
  u <- rvinecop(5, trunc_vine)
  pdf_full <- dvinecop(u, trunc_vine, keep_all = TRUE)
  expect_pdf_full_triangular_vectors(pdf_full, trunc_vine, nrow(u))
  expect_eql(pdf_full$pdf, dvinecop(u, trunc_vine))

  # summary table is truncated too
  expect_s3_class(summary(vinecop_dist(pcs[-2], mat)), "summary_df")
  expect_silent(smr <- summary(vinecop_dist(pcs[-2], mat)))
  expect_eql(nrow(smr), 2)
})

test_that("print/summary/dim generics work", {
  expect_output(print(vc))
  expect_s3_class(summary(vc), "summary_df")
  expect_silent(s <- summary(vc))
  expect_is(s, "data.frame")
  expect_eql(nrow(s), 3)
  expect_eql(ncol(s), 10)
  expect_equiv(dim(vc), c(3, 2))
})

test_that("plot functions work", {
  pcs <- lapply(
    1:4,
    function(j) {
      # pair-copulas in tree j
      lapply(runif(5 - j), function(cor) bicop_dist("gaussian", 0, cor))
    }
  )
  mat <- matrix(
    c(
      1,
      2,
      3,
      4,
      5,
      1,
      2,
      3,
      4,
      0,
      1,
      2,
      3,
      0,
      0,
      1,
      2,
      0,
      0,
      0,
      1,
      0,
      0,
      0,
      0
    ),
    5,
    5
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
  expect_equiv(as_rvine_matrix(mat), get_matrix(vc))
  expect_error(get_matrix(12))

  # test get_pair_copulas
  expect_silent(pcc <- get_pair_copula(vc, 1, 1))
  expect_eql(bicop, bicop_dist(pcc$family, pcc$rotation, pcc$parameters))
  expect_error(get_pair_copula(12, 1, 1))
  expect_error(get_pair_copula(vc, 1:2, 1))
  expect_error(get_pair_copula(vc, 1, 1:2))
  expect_error(get_pair_copula(vc, 0, 1))
  expect_error(get_pair_copula(vc, 1, 0))
  expect_error(get_pair_copula(vc, 12, 1))
  expect_error(get_pair_copula(vc, 1, 12))

  # test get_all_pair_copulas
  expect_equiv(pcs, get_all_pair_copulas(vc))
  expect_equiv(pcs[1:2], get_all_pair_copulas(vc, 1:2))
  expect_error(get_all_pair_copulas(12))
  expect_error(get_all_pair_copulas(vc, 0))
  expect_error(get_all_pair_copulas(vc, 12))

  # test other getters
  expect_equiv(get_parameters(vc, 1, 1), coef(pcs[[1]][[1]]))
  expect_equiv(
    get_all_parameters(vc),
    lapply(pcs, function(tree) lapply(tree, coef))
  )
  expect_equiv(get_ktau(vc, 1, 1), par_to_ktau(bicop))
  expect_equiv(
    get_all_ktaus(vc),
    lapply(pcs, function(tree) lapply(tree, function(pc) par_to_ktau(pc)))
  )
  expect_equiv(get_family(vc, 1, 1), "bb1")
  expect_equiv(
    get_all_families(vc),
    lapply(pcs, function(tree) lapply(tree, function(pc) pc$family))
  )

  # test printed output of getters
  expect_output(print(get_all_pair_copulas(vc)))
  expect_output(print(get_all_pair_copulas(vc, 1)))
})


test_that("d = 1 works", {
  vc <- vinecop_dist(list(), rvine_structure(1))
  u <- runif(5)
  expect_eql(dim(summary(vc))[1], 0)

  expect_eql(dvinecop(u, vc), rep(1, 5))
  expect_eql(pvinecop(u, vc), u, tol = 1e-2)
  expect_eql(c(rosenblatt(u, vc)), u)
  expect_eql(c(inverse_rosenblatt(u, vc)), u)
  expect_silent(rvinecop(10, vc))
  expect_error(plot(vc))
  expect_warning(contour(vc))
})
