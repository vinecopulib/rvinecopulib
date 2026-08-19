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
    expected_length <- if (array %in% sub_arrays && all(vc$var_types == "c")) {
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
  expect_equal(dvinecop(u, vc, 2, TRUE), pdf_full)
  expect_gte(min(pvinecop(u, vc, 100)), 0)
  expect_lte(max(pvinecop(u, vc, 100)), 1)
})

test_that("conditional simulation handles R-side ordering and recycling", {
  indep <- bicop_dist()
  pcs_indep <- list(
    rep(list(indep), 3),
    rep(list(indep), 2),
    list(indep)
  )
  vc_cond <- vinecop_dist(pcs_indep, dvine_structure(1:4))
  vc_cond$names <- letters[1:4]
  structure_before <- vc_cond$structure

  # Values follow the explicitly supplied set order (4, 3), even though the
  # current tail order is (3, 4).
  set.seed(12)
  u <- rvinecop(
    20,
    vc_cond,
    u_cond = c(0.8, 0.2),
    conditioning_set = c("d", "c")
  )
  expect_equal(dim(u), c(20L, 4L))
  expect_identical(colnames(u), letters[1:4])
  expect_equal(u[, 4], rep(0.8, 20))
  expect_equal(u[, 3], rep(0.2, 20))
  expect_identical(vc_cond$structure, structure_before)

  set.seed(12)
  u_parallel <- rvinecop(
    20,
    vc_cond,
    u_cond = c(0.8, 0.2),
    conditioning_set = c("d", "c"),
    cores = 2
  )
  expect_equal(u_parallel, u)

  # Observation-specific conditions are accepted without recycling.
  u_cond <- cbind(seq(0.1, 0.5, length.out = 5), rep(0.7, 5))
  u_obs <- rvinecop(5, vc_cond, u_cond = u_cond)
  expect_equal(u_obs[, 3], u_cond[, 1])
  expect_equal(u_obs[, 4], u_cond[, 2])

  # A different admissible set is handled by transient reorientation.
  u_reoriented <- rvinecop(
    5,
    vc_cond,
    u_cond = 0.35,
    conditioning_set = "a"
  )
  expect_equal(u_reoriented[, 1], rep(0.35, 5))

  expect_error(
    rvinecop(5, vc_cond, u_cond = matrix(0.5, 2, 2)),
    "must have one or 'n' rows"
  )
  expect_error(
    rvinecop(
      5,
      vc_cond,
      u_cond = c(0.2, 0.3),
      conditioning_set = "a"
    ),
    "continuous conditioning variables"
  )
  expect_error(
    rvinecop(5, vc_cond, conditioning_set = "a"),
    "requires 'u_cond'"
  )
  expect_error(rvinecop(1.5, vc_cond), "finite positive whole number")
  expect_error(rvinecop(5, vc_cond, qrng = 1), "not a flag")
  expect_error(rvinecop(5, vc_cond, cores = 0), "finite positive whole number")
  expect_error(
    rvinecop(5, vc_cond, u_cond = list(0.5)),
    "must be a vector, matrix, or data frame"
  )
  expect_error(
    rvinecop(5, vc_cond, u_cond = data.frame(value = "a")),
    "must be numeric"
  )
  expect_error(
    rvinecop(5, vc_cond, u_cond = 1.2, conditioning_set = "a"),
    "all data must be contained"
  )
})

test_that("conditional simulation is synchronized with the R seed", {
  indep <- bicop_dist()
  vc_cond <- vinecop_dist(
    list(rep(list(indep), 2), list(indep)),
    dvine_structure(1:3)
  )

  for (use_qrng in c(FALSE, TRUE)) {
    set.seed(314)
    u1 <- rvinecop(
      25,
      vc_cond,
      qrng = use_qrng,
      u_cond = 0.4,
      conditioning_set = 3
    )
    state_after <- .Random.seed

    set.seed(314)
    u2 <- rvinecop(
      25,
      vc_cond,
      qrng = use_qrng,
      u_cond = 0.4,
      conditioning_set = 3
    )
    expect_identical(u2, u1)
    expect_identical(.Random.seed, state_after)

    # get_seeds() advances R's RNG by exactly its 20 seed draws; the backend
    # consumes only the resulting C++ seeds and never touches R's RNG directly.
    set.seed(314)
    runif(20)
    expect_identical(.Random.seed, state_after)

    set.seed(315)
    u3 <- rvinecop(
      25,
      vc_cond,
      qrng = use_qrng,
      u_cond = 0.4,
      conditioning_set = 3
    )
    expect_false(identical(u3[, 1:2], u1[, 1:2]))
  }
})

test_that("conditional simulation accepts discrete left limits", {
  indep <- bicop_dist()
  pcs_indep <- list(
    rep(list(indep), 3),
    rep(list(indep), 2),
    list(indep)
  )
  vc_disc <- vinecop_dist(
    pcs_indep,
    dvine_structure(1:4),
    var_types = c("c", "c", "c", "d")
  )
  vc_disc$names <- letters[1:4]

  # Expanded layout in the supplied set order (4, 3): F4, F3, F4-, F3-.
  set.seed(12)
  u_expanded <- rvinecop(
    20,
    vc_disc,
    u_cond = c(0.8, 0.7, 0.6, 0.7),
    conditioning_set = c("d", "c")
  )
  # Compact layout omits the redundant F3- column.
  set.seed(12)
  u_compact <- rvinecop(
    20,
    vc_disc,
    u_cond = c(0.8, 0.7, 0.6),
    conditioning_set = c("d", "c")
  )
  expect_identical(u_expanded, u_compact)
  expect_equal(dim(u_expanded), c(20L, 4L))
  expect_true(all(u_expanded[, 4] >= 0.6))
  expect_true(all(u_expanded[, 4] <= 0.8))
  expect_equal(u_expanded[, 3], rep(0.7, 20))

  expect_error(
    rvinecop(
      5,
      vc_disc,
      u_cond = c(0.8, 0.7),
      conditioning_set = c(4, 3)
    ),
    "additional left-limit column"
  )
  expect_error(
    rvinecop(
      5,
      vc_disc,
      u_cond = c(0.7, 0.8),
      conditioning_set = 4
    ),
    "left-limit columns.*must not exceed"
  )
})

test_that("dvinecop accepts full-vine parameter matrices", {
  u_vectorized <- rvinecop(20, vc)
  parameters <- matrix(
    rep(unlist(get_all_parameters(vc)), each = nrow(u_vectorized)),
    nrow = nrow(u_vectorized)
  )

  expect_equal(
    dvinecop(u_vectorized, vc, parameters = parameters),
    dvinecop(u_vectorized, vc)
  )
  expect_equal(
    dvinecop(
      u_vectorized,
      vc,
      keep_all = TRUE,
      parameters = parameters
    ),
    dvinecop(u_vectorized, vc, keep_all = TRUE)
  )
})

test_that("dvinecop supports observation-specific parameters", {
  set.seed(17)
  n <- 30
  u_vectorized <- matrix(runif(2 * n, 0.1, 0.9), ncol = 2)
  cop <- bicop_dist("gaussian", parameters = 0)
  vc_vectorized <- vinecop_dist(
    list(list(cop)),
    dvine_structure(1:2)
  )
  parameters <- seq(-0.7, 0.7, length.out = n)

  density <- dvinecop(
    u_vectorized,
    vc_vectorized,
    parameters = parameters
  )
  expect_equal(
    density,
    dbicop(u_vectorized, "gaussian", parameters = parameters)
  )
  expect_equal(
    dvinecop(u_vectorized, vc_vectorized, cores = 2, parameters = parameters),
    density
  )

  full <- dvinecop(
    u_vectorized,
    vc_vectorized,
    keep_all = TRUE,
    parameters = parameters
  )
  expect_equal(full$pdf, density)
  expect_equal(full$pdf_edges[[1]][[1]], density)
  expect_pdf_full_triangular_vectors(full, vc_vectorized, n)
  expect_equal(
    dvinecop(
      u_vectorized,
      vc_vectorized,
      cores = 2,
      keep_all = TRUE,
      parameters = parameters
    ),
    full
  )
})

test_that("dvinecop parameter safeguards are enforced", {
  u_vectorized <- matrix(c(0.2, 0.3, 0.7, 0.8), ncol = 2)
  cop <- bicop_dist("gaussian", parameters = 0.4)
  vc_vectorized <- vinecop_dist(list(list(cop)), dvine_structure(1:2))

  expect_error(
    dvinecop(u_vectorized, vc_vectorized, parameters = "bad"),
    "not a numeric"
  )
  expect_error(
    dvinecop(u_vectorized, vc_vectorized, parameters = matrix(0.2, 1, 1)),
    "one row per row of u"
  )
  expect_error(
    dvinecop(
      u_vectorized,
      vc_vectorized,
      parameters = matrix(rep(0.2, 4), nrow = 2)
    ),
    "get_npars.*columns"
  )
  expect_error(
    dvinecop(
      u_vectorized,
      vc_vectorized,
      parameters = matrix(c(0.2, NA_real_), ncol = 1)
    ),
    "must not contain NaN or Inf"
  )
  expect_error(
    dvinecop(
      u_vectorized,
      vc_vectorized,
      parameters = matrix(c(0.2, 1.1), ncol = 1)
    ),
    "out of bounds"
  )

  vc_discrete <- vc_vectorized
  vc_discrete$var_types[1] <- "d"
  u_discrete <- cbind(u_vectorized, u_vectorized[, 1])
  expect_error(
    dvinecop(u_discrete, vc_discrete, parameters = matrix(c(0.2, 0.3))),
    "continuous"
  )

  vc_nonparametric <- vinecop_dist(
    list(list(bicop_dist("tll"))),
    dvine_structure(1:2)
  )
  expect_error(
    dvinecop(
      u_vectorized,
      vc_nonparametric,
      parameters = matrix(numeric(60), nrow = 2)
    ),
    "parametric pair copulas"
  )
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
  expect_equal(scores(u, vc, FALSE, 2), s_full)
  expect_equal(hessian(u, vc, cores = 2), h)
  expect_equal(hessian(u, vc, step_wise = FALSE, cores = 2), h_full)
  expect_equal(hessian(u, vc, FALSE, 2), h_full)
  expect_error(scores(u, vc, step_wise = 1))
  expect_error(scores(u, vc, cores = 0))
  expect_error(hessian(u, vc, step_wise = 1))
  expect_error(hessian(u, vc, cores = 0))

  parameters <- matrix(
    rep(unlist(get_all_parameters(vc)), each = nrow(u)),
    nrow = nrow(u)
  )
  expect_equal(scores(u, vc, parameters = parameters), s)
  expect_equal(hessian(u, vc, parameters = parameters), h)
  expect_equal(scores(u, vc, parameters = parameters, cores = 2), s)
  expect_equal(hessian(u, vc, parameters = parameters, cores = 2), h)

  expect_error(scores(u, vc, parameters = "bad"), "not a numeric")
  expect_error(
    scores(u, vc, parameters = parameters[-1, , drop = FALSE]),
    "one row per row of u"
  )
  expect_error(
    hessian(u, vc, parameters = parameters[, -1, drop = FALSE]),
    "get_npars.*columns"
  )
  parameters_na <- parameters
  parameters_na[1, 1] <- NA_real_
  expect_error(
    scores(u, vc, parameters = parameters_na),
    "must not contain NaN or Inf"
  )
  vc_discrete <- vc
  vc_discrete$var_types[1] <- "d"
  expect_error(
    hessian(u, vc_discrete, parameters = parameters),
    "continuous"
  )
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
