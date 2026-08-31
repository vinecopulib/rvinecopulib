# fixes problems with change in all.equal() behavior in R 4.1.x
expect_eql <- function(...) expect_equal(..., check.environment = FALSE)
expect_equiv <- function(...) expect_equivalent(..., check.environment = FALSE)

context("(Inverse) Rosenblatt transform")

pc <- bicop_dist("bb1", 90, c(3, 2))
pcs <- list(list(pc, pc), list(pc))
mat <- matrix(c(1, 2, 3, 1, 2, 0, 1, 0, 0), 3, 3)
vc <- vinecop_dist(pcs, mat)
vd <- vine_dist(list(list(distr = "norm")), pcs, mat)

test_that("rosenblatt works with bivariate copulas", {
  u <- rbicop(20, pc)
  expect_eql(inverse_rosenblatt(rosenblatt(u, pc), pc), u)
  pc <- bicop(u, family = "clay")
  expect_eql(inverse_rosenblatt(rosenblatt(u, pc), pc), u)
})

test_that("rosenblatt works with vine copulas", {
  u <- rvinecop(20, vc)
  expect_eql(inverse_rosenblatt(rosenblatt(u, vc), vc), u)
  vc <- vinecop(u, structure = mat, family = "clay")
  expect_eql(inverse_rosenblatt(rosenblatt(u, vc), vc), u)
})

test_that("rosenblatt supports explicit conditioning sets", {
  indep <- bicop_dist()
  vc_cond <- vinecop_dist(
    list(rep(list(indep), 3), rep(list(indep), 2), list(indep)),
    dvine_structure(1:4)
  )
  vc_cond$names <- letters[1:4]
  structure_before <- vc_cond$structure
  u <- rvinecop(20, vc_cond)

  z_names <- rosenblatt(
    u,
    vc_cond,
    conditioning_set = c("d", "c")
  )
  z_indices <- rosenblatt(u, vc_cond, conditioning_set = c(4, 3))
  expect_identical(z_names, z_indices)
  expect_identical(colnames(z_names), letters[1:4])
  expect_equal(
    inverse_rosenblatt(
      z_names,
      vc_cond,
      conditioning_set = c("d", "c")
    ),
    u
  )
  expect_identical(
    rosenblatt(u, vc_cond, cores = 2, conditioning_set = c(4, 3)),
    z_indices
  )
  expect_identical(
    inverse_rosenblatt(
      z_indices,
      vc_cond,
      cores = 2,
      conditioning_set = c(4, 3)
    ),
    u
  )
  expect_identical(vc_cond$structure, structure_before)

  vc_truncated <- truncate_model(vc_cond, 1)
  truncated_structure_before <- vc_truncated$structure
  u_truncated <- rvinecop(20, vc_truncated)
  z_truncated <- rosenblatt(
    u_truncated,
    vc_truncated,
    conditioning_set = c(4, 3)
  )
  expect_equal(
    inverse_rosenblatt(
      z_truncated,
      vc_truncated,
      conditioning_set = c(4, 3)
    ),
    u_truncated
  )
  expect_identical(vc_truncated$structure, truncated_structure_before)

  u_bicop <- rbicop(20, pc)
  expect_equal(
    inverse_rosenblatt(
      rosenblatt(u_bicop, pc, conditioning_set = 1),
      pc,
      conditioning_set = 1
    ),
    u_bicop
  )
})

test_that("conditional rosenblatt safeguards are enforced", {
  indep <- bicop_dist()
  vc_cond <- vinecop_dist(
    list(rep(list(indep), 2), list(indep)),
    dvine_structure(1:3)
  )
  u <- matrix(runif(15), ncol = 3)

  expect_error(
    rosenblatt(u, vc_cond, conditioning_set = "a"),
    "requires named variables"
  )
  vc_cond$names <- letters[1:3]
  expect_error(
    rosenblatt(u, vc_cond, conditioning_set = "unknown"),
    "unknown variable names"
  )
  expect_error(
    inverse_rosenblatt(u, vc_cond, conditioning_set = c(1, 1)),
    "must not contain duplicates"
  )
  expect_error(
    rosenblatt(u, vc_cond, conditioning_set = 0),
    "between 1 and d"
  )
  expect_error(
    inverse_rosenblatt(u, vc_cond, conditioning_set = 1:3),
    "at most d - 1"
  )
  expect_error(
    rosenblatt(u, vc_cond, conditioning_set = 1.5),
    "indices or names"
  )
  expect_error(
    rosenblatt(u, vc_cond, cores = 0),
    "finite positive whole number"
  )
  expect_error(
    inverse_rosenblatt(u, vc_cond, cores = 0),
    "finite positive whole number"
  )
  expect_error(
    rosenblatt(u, vc_cond, randomize_discrete = 1),
    "not a flag"
  )
})

test_that("conditional discrete rosenblatt uses R seeds and both layouts", {
  indep <- bicop_dist()
  vc_disc <- vinecop_dist(
    list(rep(list(indep), 2), list(indep)),
    dvine_structure(1:3),
    var_types = c("c", "d", "c")
  )
  u <- matrix(seq(0.15, 0.85, length.out = 30), ncol = 3)
  u_lower <- pmax(u[, 2] - 0.1, 1e-10)
  compact <- cbind(u, u_lower)
  expanded <- cbind(u, u)
  expanded[, 5] <- u_lower

  set.seed(314)
  z_compact <- rosenblatt(compact, vc_disc, conditioning_set = 2)
  state_after <- .Random.seed
  set.seed(314)
  z_expanded <- rosenblatt(expanded, vc_disc, conditioning_set = 2)
  expect_identical(z_expanded, z_compact)
  expect_identical(.Random.seed, state_after)

  set.seed(314)
  runif(20)
  expect_identical(.Random.seed, state_after)

  set.seed(315)
  z_other_seed <- rosenblatt(compact, vc_disc, conditioning_set = 2)
  expect_false(identical(z_other_seed[, 2], z_compact[, 2]))
})

test_that("discrete rosenblatt works with vine copulas", {
  u <- rvinecop(2000, vc)
  uu <- cbind(u, u[, 2])

  thresh <- 0.05
  uu[u[, 2] < thresh, 2] <- 1e-10
  uu[u[, 2] < thresh, 4] <- thresh

  vc_c <- vc <- vinecop(
    uu,
    var_types = c("c", "d", "c"),
    structure = mat,
    family = "clay"
  )
  vc_c$var_types <- rep("c", 3)

  v <- inverse_rosenblatt(rosenblatt(uu, vc), vc_c)
  expect_eql(v, u, tol = 2 * thresh)

  # other format
  uu <- cbind(u, u)
  uu[u[, 2] < thresh, 2] <- 1e-10
  uu[u[, 2] < thresh, 5] <- thresh
  v <- inverse_rosenblatt(rosenblatt(uu, vc), vc_c)
  expect_eql(v, u, tol = 2 * thresh)
})

test_that("rosenblatt works with vine distribution", {
  u <- rvine(20, vd)
  expect_eql(inverse_rosenblatt(rosenblatt(u, vd), vd), u)
  vd <- vine(u, copula_controls = list(structure = mat, family = "clay"))
  expect_equiv(inverse_rosenblatt(rosenblatt(u, vd), vd), u)

  vd$names <- vd$copula$names <- letters[1:3]
  colnames(u) <- letters[1:3]
  z <- rosenblatt(u, vd, conditioning_set = "a")
  expect_identical(colnames(z), letters[1:3])
  expect_equiv(
    inverse_rosenblatt(z, vd, conditioning_set = "a"),
    u
  )
})

test_that("discrete rosenblatt works with vine distributions", {
  x <- data.frame(
    x1 = as.ordered(sample(1:4, 50, replace = TRUE)),
    x2 = rnorm(50),
    x3 = rbinom(50, 3, 0.5)
  )
  expect_warning(
    vd <- vine(
      x,
      var_types = c("d", "c", "zi"),
      copula_controls = list(family_set = "gauss")
    ),
    "AIC and BIC are unavailable"
  )
  u <- rvine(50, vd)
  expect_true(mean(u[, 3] == 0) > 0.02)
  expect_eql(colnames(rosenblatt(x, vd)), colnames(x))
})
