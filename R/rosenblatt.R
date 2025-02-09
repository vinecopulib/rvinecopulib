#' (Inverse) Rosenblatt transform
#'
#' The Rosenblatt transform takes data generated from a model and turns it into
#' independent uniform variates, The inverse Rosenblatt transform computes
#' conditional quantiles and can be used simulate from a stochastic model,
#' see *Details*.
#'
#' @name rosenblatt
#' @aliases rosenblatt inverse_rosenblatt
#'
#' @param x matrix of evaluation points; must be in \eqn{(0, 1)^d} for copula
#'   models.
#' @param u matrix of evaluation points; must be in \eqn{(0, 1)^d}.
#' @param model a model object; classes currently supported are
#'    `bicop_dist()`, `vinecop_dist()`, and `vine_dist()`.
#' @param cores if `>1`, computation is parallelized over `cores` batches (rows
#'    of `u`).
#' @param randomize_discrete Whether to randomize the transform for discrete
#'   variables; see Details.
#'
#' @details
#' The Rosenblatt transform (Rosenblatt, 1952) \eqn{U = T(V)} of a random vector
#' \eqn{V = (V_1,\ldots,V_d) ~ F} is defined as
#' \deqn{
#'   U_1= F(V_1), U_{2} = F(V_{2}|V_1), \ldots, U_d =F(V_d|V_1,\ldots,V_{d-1}),
#' }
#' where \eqn{F(v_k|v_1,\ldots,v_{k-1})} is the conditional distribution of
#' \eqn{V_k} given \eqn{V_1 \ldots, V_{k-1}, k = 2,\ldots,d}. The vector
#' \eqn{U  = (U_1, \dots, U_d)} then contains independent standard uniform
#' variables. The inverse operation
#' \deqn{
#'   V_1 = F^{-1}(U_1), V_{2} = F^{-1}(U_2|U_1), \ldots,
#'   V_d =F^{-1}(U_d|U_1,\ldots,U_{d-1}),
#' }
#' can be used to simulate from a distribution. For any copula \eqn{F}, if
#' \eqn{U} is a vector of independent random variables, \eqn{V = T^{-1}(U)} has
#' distribution \eqn{F}.
#'
#' The formulas above assume a vine copula model with order \eqn{d, \dots, 1}.
#' More generally, `rosenblatt()` returns the variables
#' \deqn{
#'   U_{M[d + 1- j, j]}= F(V_{M[d - j + 1, j]} | V_{M[d - j, j]}, \dots, V_{M[1, j]}),
#' }
#' where \eqn{M} is the structure matrix. Similarly, `inverse_rosenblatt()`
#' returns
#' \deqn{
#'   V_{M[d + 1- j, j]}= F^{-1}(U_{M[d - j + 1, j]} | U_{M[d - j, j]}, \dots, U_{M[1, j]}).
#' }
#'
#' If some variables have atoms, Brockwell (10.1016/j.spl.2007.02.008) proposed
#' a simple randomization scheme to ensure that output is still independent
#' uniform if the model is correct. The transformation reads
#' \deqn{ U_{M[d - j,
#' j]}= W_{d - j} F(V_{M[d - j, j]} | V_{M[d - j - 1, j - 1]}, \dots, V_{M[0,
#' 0]}) + (1 - W_{d - j}) F^-(V_{M[d - j, j]} | V_{M[d - j - 1, j - 1]}, \dots,
#' V_{M[0, 0]}),
#' }
#' where \eqn{F^-}
#' is the left limit of the conditional cdf
#' and \eqn{W_1, \dots, W_d} are are independent standard uniform random
#' variables. This is used by default. If you are interested in the conditional
#' probabilities
#' \deqn{
#'  F(V_{M[d - j, j]} | V_{M[d - j - 1, j - 1]}, \dots, V_{M[0, 0]}),
#' }
#' set `randomize_discrete = FALSE`.
#'
#' @examples
#' # simulate data with some dependence
#' x <- replicate(3, rnorm(200))
#' x[, 2:3] <- x[, 2:3] + x[, 1]
#' pairs(x)
#'
#' # estimate a vine distribution model
#' fit <- vine(x, copula_controls = list(family_set = "par"))
#'
#' # transform into independent uniforms
#' u <- rosenblatt(x, fit)
#' pairs(u)
#'
#' # inversion
#' pairs(inverse_rosenblatt(u, fit))
#'
#' # works similarly for vinecop models
#' vc <- fit$copula
#' rosenblatt(pseudo_obs(x), vc)
#' @export
rosenblatt <- function(x, model, cores = 1, randomize_discrete = TRUE) {
  assert_that(
    inherits(model, c("bicop_dist", "vinecop_dist", "vine_dist")),
    is.number(cores)
  )

  if (inherits(model, "bicop_dist")) {
    model <- vinecop_dist(
      list(list(model)),
      cvine_structure(1:2),
      var_types = model$var_types
    )
  }

  if (inherits(model, "vine_dist")) {
    x <- expand_factors(x)
    if (!is.null(model$names)) {
      x <- x[, model$names, drop = FALSE]
    }
    x <- compute_pseudo_obs(x, model)
    model <- model$copula
  }

  # model is now a vinecop_dist
  assert_that(all((x >= 0) & (x <= 1)))
  x <- pmin(pmax(x, 1e-10), 1 - 1e-10)
  x <- if_vec_to_matrix(x, dim(model)[1] == 1)
  x <- vinecop_rosenblatt_cpp(x, model, cores, randomize_discrete, get_seeds())
  colnames(x) <- model$names

  x
}

#' @rdname rosenblatt
#' @export
inverse_rosenblatt <- function(u, model, cores = 1) {
  assert_that(
    all((u > 0) & (u < 1)),
    inherits(model, c("bicop_dist", "vinecop_dist", "vine_dist")),
    is.number(cores)
  )

  to_col <- if (inherits(model, "bicop_dist")) FALSE else (dim(model)[1] == 1)
  u <- if_vec_to_matrix(u, to_col)

  if (inherits(model, "bicop_dist")) {
    model <- vinecop_dist(
      list(list(model)),
      cvine_structure(1:2),
      var_types = model$var_types
    )
  }

  if (inherits(model, "vinecop_dist")) {
    u <- vinecop_inverse_rosenblatt_cpp(u, model, cores)
  } else {
    u <- vinecop_inverse_rosenblatt_cpp(u, model$copula, cores)
    u <- dpq_marg(u, model, "q")
  }
  colnames(u) <- model$names

  u
}
