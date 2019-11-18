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
#'
#' @details
#' The Rosenblatt transform (Rosenblatt, 1952) \eqn{U = T(V)} of a random vector
#' \eqn{V = (V_1,\ldots,V_d) ~ F} is defined as
#' \deqn{
#'   U_1= F(V_1), U_{2} = F(V_{2}|V_1), \ldots, U_d =F(V_d|V_1,\ldots,V_{d-1}),
#' }
#' where \eqn{F(v_k|v_1,\ldots,v_{k-1})} is the conditional distribution of
#' \eqn{V_k} given \eqn{V_1 \ldots, V_{k-1}, k = 2,\ldots,d}. The vector \eqn{U}
#' are then independent standard uniform variables. The inverse operation
#' \deqn{
#'   V_1 = F^{-1}(U_1), V_{2} = F^{-1}(U_2|U_1), \ldots,
#'   V_d =F^{-1}(U_d|U_1,\ldots,U_{d-1}),
#' }
#' can can be used to simulate from a distribution. For any copula \eqn{F}, if
#' \eqn{U} is a vector of independent random variables, \eqn{V = T^{-1}(U)} has
#' distribution \eqn{F}.
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
rosenblatt <- function(x, model, cores = 1) {
  assert_that(
    is.matrix(x) | is.data.frame(x),
    inherits(model, c("bicop_dist", "vinecop_dist", "vine_dist")),
    is.number(cores)
  )
  col_names <- colnames(x)

  if (inherits(model, "bicop_dist")) {
    assert_that(ncol(x) == 2, all((x > 0) & (x < 1)))
    x <- cbind(x[, 1], bicop_hfunc1_cpp(x, model))
  } else if (inherits(model, "vinecop_dist")) {
    assert_that(all((x > 0) & (x < 1)))
    x <- vinecop_rosenblatt_cpp(x, model, cores)
  } else {
    # prepare marginals if only one is specified
    if (!inherits(vine, "vine") & depth(model$margins) == 1) {
      model$margins <- replicate(dim(model)[1], model$margins, simplify = FALSE)
    }
    x <- dpq_marg(x, model, "p")
    x <- vinecop_rosenblatt_cpp(x, model$copula, cores)
  }
  colnames(x) <- col_names

  x
}

#' @rdname rosenblatt
#' @export
inverse_rosenblatt <- function(u, model, cores = 1) {
  assert_that(
    is.matrix(u) | is.data.frame(u),
    all((u > 0) & (u < 1)),
    inherits(model, c("bicop_dist", "vinecop_dist", "vine_dist")),
    is.number(cores)
  )
  col_names <- colnames(u)

  if (inherits(model, "bicop_dist")) {
    assert_that(ncol(u) == 2)
    u <- cbind(u[, 1], bicop_hinv1_cpp(u, model))
  } else if (inherits(model, "vinecop_dist")) {
    u <- vinecop_inverse_rosenblatt_cpp(u, model, cores)
  } else {
    u <- vinecop_inverse_rosenblatt_cpp(u, model$copula, cores)
    # prepare marginals if only one is specified
    if (!inherits(vine, "vine") & depth(model$margins) == 1) {
      model$margins <- replicate(dim(model)[1], model$margins, simplify = FALSE)
    }
    u <- dpq_marg(u, model, "q")
  }
  colnames(u) <- col_names

  u
}
