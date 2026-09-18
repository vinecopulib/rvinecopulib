#' Rosenblatt and inverse Rosenblatt transforms
#'
#' The Rosenblatt transform maps observations from a model to independent
#' uniform variates. The inverse transform evaluates conditional quantiles and
#' maps independent uniforms to draws from the model; see *Details*.
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
#' @param conditioning_set optional variable indices or names that define the
#'   conditioning variables. The transform uses an admissible sampling order
#'   whose tail contains exactly this set. The model is not modified. If `NULL`,
#'   the current model order is used. An error is thrown if the requested set
#'   cannot form an admissible tail.
#'
#' @details
#' Let \eqn{(s_1, \ldots, s_d)} be the sampling order of a random vector
#' \eqn{V = (V_1, \ldots, V_d)} with distribution \eqn{F}. For continuous
#' variables, the Rosenblatt transform \eqn{Z = T(V)} is defined by
#' \deqn{
#'   Z_{s_1} = F_{s_1}(V_{s_1}), \qquad
#'   Z_{s_j} = F_{s_j \mid s_1, \ldots, s_{j-1}}
#'   (V_{s_j} \mid V_{s_1}, \ldots, V_{s_{j-1}}),
#'   \quad j = 2, \ldots, d.
#' }
#' If the model is correct, the components of \eqn{Z} are independent standard
#' uniforms. The result is stored in the original variable columns; the
#' sampling order determines only the sequence of conditional distributions.
#'
#' The inverse transform applies the corresponding conditional quantiles:
#' \deqn{
#'   V_{s_1} = F_{s_1}^{-1}(Z_{s_1}), \qquad
#'   V_{s_j} = F_{s_j \mid s_1, \ldots, s_{j-1}}^{-1}
#'   (Z_{s_j} \mid V_{s_1}, \ldots, V_{s_{j-1}}),
#'   \quad j = 2, \ldots, d.
#' }
#' Thus \eqn{T^{-1}(Z)} has distribution \eqn{F} when \eqn{Z} contains
#' independent standard uniforms.
#'
#' If a variable has atoms, its conditional cdf jumps at the observation. Let
#' \eqn{G_j} denote the conditional cdf of \eqn{V_{s_j}} given the preceding
#' variables in the sampling order. Following Brockwell
#' (10.1016/j.spl.2007.02.008), `rosenblatt()` returns
#' \deqn{
#'   Z_{s_j} = W_j G_j(V_{s_j}) + (1 - W_j) G_j(V_{s_j}^{-}),
#' }
#' where \eqn{G_j(V_{s_j}^{-})} is the left limit and the \eqn{W_j} are
#' independent standard uniforms. This randomization is used by default and
#' yields uniform components under the fitted model. Set
#' `randomize_discrete = FALSE` to return the upper endpoint
#' \eqn{G_j(V_{s_j})} instead.
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
rosenblatt <- function(
  x,
  model,
  cores = 1,
  randomize_discrete = TRUE,
  conditioning_set = NULL
) {
  cores <- as_count(cores, "cores")
  assert_that(
    inherits(model, c("bicop_dist", "vinecop_dist", "vine_dist")),
    is.flag(randomize_discrete)
  )

  conditioning_set <- process_conditioning_set(
    conditioning_set,
    model$names,
    if (inherits(model, "bicop_dist")) 2L else dim(model)[1]
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
  x <- vinecop_rosenblatt_cpp(
    x,
    model,
    conditioning_set,
    cores,
    randomize_discrete,
    get_seeds()
  )
  colnames(x) <- model$names

  x
}

#' @rdname rosenblatt
#' @export
inverse_rosenblatt <- function(
  u,
  model,
  cores = 1,
  conditioning_set = NULL
) {
  cores <- as_count(cores, "cores")
  assert_that(
    all((u > 0) & (u < 1)),
    inherits(model, c("bicop_dist", "vinecop_dist", "vine_dist"))
  )

  conditioning_set <- process_conditioning_set(
    conditioning_set,
    model$names,
    if (inherits(model, "bicop_dist")) 2L else dim(model)[1]
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
    u <- vinecop_inverse_rosenblatt_cpp(u, model, conditioning_set, cores)
  } else {
    u <- vinecop_inverse_rosenblatt_cpp(
      u,
      model$copula,
      conditioning_set,
      cores
    )
    u <- dpq_marg(u, model, "q")
  }
  colnames(u) <- model$names

  u
}
