# fixes problems with change in all.equal() behavior in R 4.1.x
expect_eql <- function(...) expect_equal(..., check.environment = FALSE)
expect_equiv <- function(...) expect_equivalent(..., check.environment = FALSE)

context("Fitting 'vinecop' models")

set.seed(5)
u <- sapply(1:5, function(i) runif(30))
fit <- vinecop(u, family = "nonpar")
fit_with_data <- vinecop(u, family = "nonpar", keep_data = TRUE)

test_that("returns proper 'vinecop' object", {
  expect_s3_class(fit, "vinecop")
  expect_s3_class(fit, "vinecop_dist")
  expect_identical(
    names(fit),
    c(
      "pair_copulas",
      "structure",
      "var_types",
      "npars",
      "loglik",
      "threshold",
      "controls",
      "nobs"
    )
  )
  expect_identical(
    names(fit_with_data),
    c(
      "pair_copulas",
      "structure",
      "var_types",
      "npars",
      "loglik",
      "threshold",
      "data",
      "controls",
      "nobs"
    )
  )

  colnames(u) <- paste(seq_len(ncol(u)))
  expect_identical(
    names(vinecop(u, family = "indep")),
    c(
      "pair_copulas",
      "structure",
      "var_types",
      "npars",
      "loglik",
      "threshold",
      "names",
      "controls",
      "nobs"
    )
  )
})

test_that("works with structure", {
  u <- sapply(1:2, function(i) runif(30))
  expect_silent(fit <- vinecop(u, structure = matrix(c(1:2, 1:0), 2, 2)))
})

if (Sys.info()["sysname"] != "SunOS") {
  test_that("runs in parallel", {
    expect_silent(fit <- vinecop(u, cores = 2))
  })
}

test_that("S3 generics work", {
  expect_eql(predict(fit, u), fitted(fit_with_data))
  expect_eql(
    predict(fit, u, what = "cdf"),
    fitted(fit_with_data, what = "cdf"),
    tolerance = 0.01
  )
  expect_error(predict(fit, u, what = "hfunc1"))
  fit$data <- NULL
  expect_error(fitted(fit))
  expect_length(attr(logLik(fit), "df"), 1)
})

test_that("print/summary generics work", {
  expect_output(print(fit))
  expect_s3_class(s <- summary(fit), "summary_df")
  expect_is(s, "data.frame")

  fit$names <- letters[1:3]
  out <- capture.output(summary(fit))
  expect_eql(out[length(out)], "1 <-> a,   2 <-> b,   3 <-> c ")

  fit$names <- c(
    "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
    "bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb",
    "cccccccccccccccccccccccccccccccccccccccccccccccc"
  )
  out <- capture.output(summary(fit))
  expect_eql(
    out[length(out)],
    "3 <-> cccccccccccccccccccccccccccccccccccccccccccccccc "
  )
})

test_that("truncation works", {
  fit_truncated <- truncate_model(fit, trunc_lvl = 1)
  expect_silent(dvinecop(u, fit_truncated))
  expect_silent(rvinecop(50, fit_truncated))

  fit_truncated <- vinecop(u, family = "par", trunc_lvl = 1)
  expect_silent(dvinecop(u, fit_truncated))
  expect_silent(rvinecop(50, fit_truncated))
})

test_that("partial selection works", {
  fit_partial <- vinecop(
    u[, sample(1:5)],
    structure = truncate_model(fit$structure, 1),
    trunc_lvl = 3
  )
  expect_eql(unname(dim(fit_partial)[2]), 3)

  m_old <- as_rvine_matrix(fit$structure)
  m_new <- as_rvine_matrix(fit_partial$structure)
  tree1_old_edges <- c(
    paste(diag(m_old[5:2, ]), m_old[1, -5]),
    paste(m_old[1, -5], diag(m_old[5:2, ]))
  )
  expect_true(all(paste(diag(m_new[5:2, ]), m_new[1, -5]) %in% tree1_old_edges))
})

test_that("MST algorithms behave as expected", {
  # Generate data
  u <- replicate(7, runif(500))

  # (a) prim and kruskal give the same structure
  set.seed(42)
  fit_prim <- vinecop(u, family_set = "indep", tree_algorithm = "mst_prim")

  set.seed(42)
  fit_kruskal <- vinecop(
    u,
    family_set = "indep",
    tree_algorithm = "mst_kruskal"
  )

  m_prim <- as_rvine_matrix(fit_prim$structure)
  m_kruskal <- as_rvine_matrix(fit_kruskal$structure)

  expect_equal(m_prim, m_kruskal)

  # (b) random gives different structures for seeds 1:10
  structures <- vector("list", 10)

  for (i in 1:10) {
    set.seed(i)
    fit_random <- vinecop(
      u,
      family_set = "indep",
      tree_algorithm = "random_weighted"
    )
    structures[[i]] <- as_rvine_matrix(fit_random$structure)
  }

  # Check uniqueness
  unique_structures <- unique(structures)
  expect_equal(length(unique_structures), 10)
})

test_that("conditioning-aware selection is exposed through the R controls", {
  u_cond <- matrix(runif(400), 100, 4)
  colnames(u_cond) <- letters[1:4]
  fit_cond <- vinecop(
    u_cond,
    family_set = "indep",
    conditioning_set = c("b", "d")
  )

  expect_identical(fit_cond$controls$conditioning_set, c(2L, 4L))
  expect_equal(sort(tail(fit_cond$structure$order, 2)), c(2L, 4L))

  expect_error(
    vinecop(u_cond, conditioning_set = c("b", "b")),
    "must not contain duplicates"
  )
  expect_error(
    vinecop(u_cond, conditioning_set = "unknown"),
    "unknown variable names"
  )
  u_duplicated_names <- u_cond
  colnames(u_duplicated_names)[2] <- colnames(u_duplicated_names)[1]
  expect_error(
    vinecop(u_duplicated_names, conditioning_set = "a"),
    "requires unique variable names"
  )
  expect_error(
    vinecop(unname(u_cond), conditioning_set = "b"),
    "requires named variables"
  )
  for (invalid_set in list(TRUE, NA_real_, Inf, 1.5)) {
    expect_error(
      vinecop(u_cond, conditioning_set = invalid_set),
      "must contain variable indices or names"
    )
  }
  for (invalid_set in list(0, 5)) {
    expect_error(
      vinecop(u_cond, conditioning_set = invalid_set),
      "indices must be between 1 and d"
    )
  }
  expect_error(
    vinecop(u_cond, conditioning_set = 1:4),
    "at most d - 1 variables"
  )
  expect_error(
    vinecop(u_cond, conditioning_set = 2, trunc_lvl = 1),
    "requires a non-truncated vine"
  )
  expect_error(
    vinecop(
      u_cond,
      conditioning_set = 2,
      tree_algorithm = "random_weighted"
    ),
    "requires 'tree_algorithm'"
  )
  expect_error(
    vinecop(u_cond, conditioning_set = 2, vinecop_object = fit_cond),
    "cannot be used when refitting"
  )
})

test_that("custom tree criteria are passed through", {
  set.seed(15)
  u_custom <- matrix(runif(400), ncol = 4)
  n_calls <- 0L
  seen_weights <- NULL
  custom_tau <- function(data, weights) {
    n_calls <<- n_calls + 1L
    seen_weights <<- weights
    -abs(cor(data[, 1], data[, 2], method = "kendall"))
  }

  fit_custom <- vinecop(
    u_custom,
    family_set = "indep",
    tree_crit = custom_tau
  )
  fit_tau <- vinecop(
    u_custom,
    family_set = "indep",
    tree_crit = "tau"
  )

  expect_gt(n_calls, 0L)
  expect_length(seen_weights, 0L)
  expect_equal(
    as_rvine_matrix(fit_custom$structure),
    as_rvine_matrix(fit_tau$structure)
  )
  expect_identical(fit_custom$controls$tree_crit, custom_tau)

  obs_weights <- seq_len(nrow(u_custom))
  vinecop(
    u_custom,
    family_set = "indep",
    weights = obs_weights,
    tree_crit = function(data, weights) {
      seen_weights <<- weights
      0.5
    }
  )
  expect_equal(seen_weights, obs_weights / mean(obs_weights))
})

test_that("custom tree criteria receive filtered data and weights", {
  set.seed(16)
  u_custom <- matrix(runif(60), ncol = 2)
  u_custom[1, 1] <- NA_real_
  obs_weights <- seq_len(nrow(u_custom))
  obs_weights[2] <- 0
  seen_data <- NULL
  seen_weights <- NULL

  vinecop(
    u_custom,
    family_set = "indep",
    weights = obs_weights,
    tree_crit = function(data, weights) {
      seen_data <<- data
      seen_weights <<- weights
      0.5
    }
  )

  standardized_weights <- obs_weights / sum(obs_weights) * length(obs_weights)
  expect_equal(dim(seen_data), c(28L, 2L))
  expect_false(anyNA(seen_data))
  expect_setequal(seen_weights, standardized_weights[-c(1, 2)])
})

test_that("custom tree criterion safeguards are enforced", {
  u_custom <- matrix(runif(90), ncol = 3)
  valid_criterion <- function(data, weights) 0.5

  expect_error(
    vinecop(u_custom, tree_crit = "custom"),
    "requires a function"
  )
  fit_custom <- vinecop(u_custom, family_set = "indep")
  expect_error(
    vinecop(
      u_custom,
      vinecop_object = fit_custom,
      tree_crit = valid_criterion
    ),
    "cannot be used when refitting"
  )
  expect_error(
    vinecop(u_custom, tree_crit = function(data, weights) c(0.1, 0.2)),
    "return one numeric value"
  )
  expect_error(
    vinecop(u_custom, tree_crit = function(data, weights) "bad"),
    "return one numeric value"
  )
  expect_error(
    vinecop(u_custom, tree_crit = function(data, weights) NA_real_),
    "return a finite numeric value"
  )
  expect_error(
    vinecop(u_custom, tree_crit = function(data, weights) Inf),
    "return a finite numeric value"
  )
  expect_error(
    vinecop(u_custom, tree_crit = function(data, weights) stop("boom")),
    "boom"
  )
})

test_that("custom tree criteria support multicore pair-copula fitting", {
  set.seed(17)
  u_custom <- matrix(runif(400), ncol = 4)
  n_calls <- 0L
  criterion <- function(data, weights) {
    n_calls <<- n_calls + 1L
    abs(cor(data[, 1], data[, 2], method = "kendall"))
  }

  fit_parallel <- vinecop(
    u_custom,
    family_set = "indep",
    tree_crit = criterion,
    cores = 2
  )
  parallel_calls <- n_calls
  expect_gt(parallel_calls, 0L)

  n_calls <- 0L
  fit_serial <- vinecop(
    u_custom,
    family_set = "indep",
    tree_crit = criterion,
    cores = 1
  )
  expect_identical(n_calls, parallel_calls)
  expect_equal(
    as_rvine_matrix(fit_parallel$structure),
    as_rvine_matrix(fit_serial$structure)
  )
})

test_that("d = 1 works", {
  vc <- vinecop(runif(20), structure = rvine_structure(1))
  vc2 <- vinecop(runif(20))
  expect_identical(vc, vc2)

  expect_eql(AIC(vc), 0)
  expect_eql(mBICV(vc), 0)
  for (psi0 in c(-1, 0, 1, 2, Inf, NA_real_)) {
    expect_error(mBICV(vc, psi0), "strictly between 0 and 1")
  }
  expect_eql(dim(summary(vc))[1], 0)
})

test_that("fitting only parameters works", {
  vc <- vinecop(u, family = "onepar")
  vc2 <- vinecop(u, vinecop_object = vc, show_trace = TRUE)
  expect_equal(vc[1:6], vc2[1:6])

  vc <- vinecop(u, family = "tll")
  vc2 <- vinecop(u, vinecop_object = vc, mult = 10)
  expect_true(
    TRUE !=
      all.equal(
        vc$pair_copulas[[1]][[1]]$parameters,
        vc2$pair_copulas[[1]][[1]]$parameters
      )
  )

  expect_warning(vinecop(u, structure = dvine_structure(3), vinecop = vc))
  expect_warning(vinecop(u, family_set = "gauss", vinecop = vc))
})
