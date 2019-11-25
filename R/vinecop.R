#' Fitting vine copula models
#'
#' Automated fitting and model selection for vine copula models with continuous
#' or discrete data.
#'
#' @inheritParams bicop
#' @param family_set a character vector of families; see [bicop()] for
#'   additional options.
#' @param structure an `rvine_structure` object, namely a compressed
#'   representation of the vine structure, or an object that can be coerced into
#'   one (see [rvine_structure()] and [as_rvine_structure()]). The dimension
#'   must be `length(pair_copulas[[1]]) + 1`; `structure = NA` performs
#'   automatic selection based on Dissman's algorithm. See *Details* for partial
#'   selection of the structure.
#' @param psi0 prior probability of a non-independence copula (only used for
#'   `selcrit = "mbic"` and `selcrit = "mbicv"`).
#' @param trunc_lvl the truncation level of the vine copula; `Inf` means no
#'   truncation, `NA` indicates that the truncation level should be selected
#'   automatically by [mBICV()].
#' @param tree_crit the criterion for tree selection, one of `"tau"`, `"rho"`,
#'   `"hoeffd"`, `"mcor"`, or `"joe"` for Kendall's \eqn{\tau}, Spearman's
#'   \eqn{\rho}, Hoeffding's \eqn{D}, maximum correlation, or logarithm of
#'   the partial correlation, respectively.
#' @param threshold for thresholded vine copulas; `NA` indicates that the
#'   threshold should be selected automatically by [mBICV()].
#' @param show_trace logical; whether a trace of the fitting progress should be
#'   printed.
#' @param cores number of cores to use; if more than 1, estimation of pair
#'   copulas within a tree is done in parallel.
#' @param var_types variable types, a length d vector; e.g., `c("c", "c")` for
#'   two continuous variables, or `c("c", "d")` for first variable continuous
#'   and second discrete.
#'
#' @details
#'
#' ## Discrete variables
#'
#' When at least one variable is discrete, two types of
#' "observations" are required in `data`: the first \eqn{n  x  d} block
#' contains realizations of \eqn{F_{X_j}(X_j)}. The second \eqn{n  x  d}
#' block contains realizations of \eqn{F_{X_j}(X_j^-)}. The minus indicates a
#' left-sided limit of the cdf. For, e.g., an integer-valued variable, it holds
#' \eqn{F_{X_j}(X_j^-) = F_{X_j}(X_j - 1)}. For continuous variables the left
#' limit and the cdf itself coincide. Respective columns can be omitted in the
#' second block.
#'
#' ## Partial structure selection
#'
#' It is possible to fix the vine structure only in the first trees and select
#' the remaining ones automatically. To specify only the first `k` trees, supply
#' a `k`-truncated `rvine_structure()` or `rvine_matrix()`. All trees up to
#' `trunc_lvl` will then be selected automatically.
#'
#' @return
#'
#' Objects inheriting from `vinecop` and `vinecop_dist` for [vinecop()]. In
#' addition to the entries provided by [vinecop_dist()], there are:
#'
#' * `threshold`, the (set or estimated) threshold used for thresholding the
#' vine.
#'
#' * `data` (optionally, if `keep_data = TRUE` was used), the dataset that was
#' passed to [vinecop()].
#'
#' * `controls`, a `list` with fit controls that was passed to [vinecop()].
#'
#' * `nobs`, the number of observations that were used to fit the model.
#'
#' @seealso [vinecop()], [dvinecop()], [pvinecop()], [rvinecop()],
#'   [plot.vinecop()], [contour.vinecop()]
#' @export
#' @examples
#' ## simulate dummy data
#' x <- rnorm(30) * matrix(1, 30, 5) + 0.5 * matrix(rnorm(30 * 5), 30, 5)
#' u <- pseudo_obs(x)
#'
#' ## fit and select the model structure, family and parameters
#' fit <- vinecop(u)
#' summary(fit)
#' plot(fit)
#' contour(fit)
#'
#' ## select by log-likelihood criterion from one-paramter families
#' fit <- vinecop(u, family_set = "onepar", selcrit = "bic")
#' summary(fit)
#'
#' ## Gaussian D-vine
#' fit <- vinecop(u, structure = dvine_structure(1:5), family = "gauss")
#' plot(fit)
#' contour(fit)
#'
#' ## Partial structure selection with only first tree specified
#' structure <- rvine_structure(order = 1:5, list(rep(5, 4)))
#' structure
#' fit <- vinecop(u, structure = structure, family = "gauss")
#' plot(fit)
#'
#' ## 1-truncated model with random structure
#' fit <- vinecop(u, structure = rvine_structure_sim(5), trunc_lvl = 1)
#' contour(fit)
#'
#' ## Model for discrete data
#' x <- qpois(u, 1)  # transform to Poisson margins
#' # we require two types of observations (see Details)
#' u_disc <- cbind(ppois(x, 1), ppois(x - 1, 1))
#' fit <- vinecop(u_disc, var_types = rep("d", 5))
#'
#' ## Model for mixed data
#' x <- qpois(u[, 1], 1)  # transform first variable to Poisson margin
#' # we require two types of observations (see Details)
#' u_disc <- cbind(ppois(x, 1), u[, 2:5], ppois(x - 1, 1))
#' fit <- vinecop(u_disc, var_types = c("d", rep("c", 4)))
vinecop <- function(data, var_types = rep("c", ncol(data)), family_set = "all",
                    structure = NA, par_method = "mle",
                    nonpar_method = "constant", mult = 1, selcrit = "bic",
                    weights = numeric(), psi0 = 0.9, presel = TRUE,
                    trunc_lvl = Inf, tree_crit = "tau", threshold = 0,
                    keep_data = FALSE, show_trace = FALSE, cores = 1) {
  assert_that(
    is.character(family_set),
    inherits(structure, "matrix") ||
      inherits(structure, "rvine_structure") ||
      (is.scalar(structure) && is.na(structure)),
    is.string(par_method),
    is.string(nonpar_method),
    is.number(mult), mult > 0,
    is.string(selcrit),
    is.numeric(weights),
    is.number(psi0), psi0 > 0, psi0 < 1,
    is.flag(presel),
    is.scalar(trunc_lvl),
    is.string(tree_crit),
    is.scalar(threshold),
    is.flag(keep_data),
    is.number(cores), cores > 0,
    correct_var_types(var_types)
  )

  # check if families known (w/ partial matching) and expand convenience defs
  family_set <- process_family_set(family_set, par_method)

  ## pre-process input
  if (is.scalar(structure) && is.na(structure)) {
    structure <- rvine_structure(seq_along(var_types))
  } else {
    structure <- as_rvine_structure(structure)
  }

  ## fit and select copula model
  vinecop <- vinecop_select_cpp(
    data = data,
    structure = structure,
    family_set = family_set,
    par_method = par_method,
    nonpar_method = nonpar_method,
    mult = mult,
    selection_criterion = selcrit,
    weights = weights,
    psi0 = psi0,
    preselect_families = presel,
    truncation_level = ifelse( # Inf cannot be passed to C++
      is.finite(trunc_lvl),
      trunc_lvl,
      .Machine$integer.max
    ),
    tree_criterion = tree_crit,
    threshold = threshold,
    select_truncation_level = is.na(trunc_lvl),
    select_threshold = is.na(threshold),
    show_trace = show_trace,
    num_threads = cores,
    var_types = var_types
  )

  ## make all pair-copulas bicop objects
  vinecop$pair_copulas <- lapply(
    vinecop$pair_copulas,
    function(tree) lapply(tree, as.bicop)
  )

  ## add information about the fit
  vinecop$names <- colnames(data)
  if (keep_data) {
    vinecop$data <- data
  }
  vinecop$controls <- list(
    family_set = family_set,
    par_method = par_method,
    nonpar_method = nonpar_method,
    mult = mult,
    selcrit = selcrit,
    weights = weights,
    presel = presel,
    trunc_lvl = trunc_lvl,
    tree_crit = tree_crit,
    threshold = threshold
  )
  vinecop$nobs <- nrow(data)
  vinecop
}

#' Vine copula models
#'
#' Create custom vine copula models by specifying the pair-copulas, structure,
#' and variable types.
#'
#' @return
#'
#' Object of class `vinecop_dist`, i.e., a list containing:
#'
#' * `pair_copulas`, a list of lists. Each element of `pair_copulas` corresponds
#' to a tree, which is itself a list of [bicop_dist()] objects.
#'
#' * `structure`, a compressed representation of the vine structure, or an
#' object that can be coerced into one (see [rvine_structure()] and
#' [as_rvine_structure()]).
#'
#' * `npars`, a `numeric` with the number of (effective) parameters.
#'
#' * `var_types` the variable types.
#'
#' @inheritParams vinecop
#' @param pair_copulas A nested list of '[bicop_dist()]' objects, where
#'   \code{pair_copulas[[t]][[e]]} corresponds to the pair-copula at edge `e` in
#'   tree `t`.
#' @seealso [rvine_structure()], [rvine_matrix()], [vinecop()],
#'   [plot.vinecop_dist()], [contour.vinecop_dist()], [dvinecop()],
#'   [pvinecop()], [rvinecop()]
#' @export
#' @examples
#' # specify pair-copulas
#' bicop <- bicop_dist("bb1", 90, c(3, 2))
#' pcs <- list(
#'   list(bicop, bicop), # pair-copulas in first tree
#'   list(bicop) # pair-copulas in second tree
#' )
#'
#' # specify R-vine matrix
#' mat <- matrix(c(1, 2, 3, 1, 2, 0, 1, 0, 0), 3, 3)
#'
#' # set up vine copula model
#' vc <- vinecop_dist(pcs, mat)
#'
#' # visualization
#' plot(vc)
#' contour(vc)
#'
#' # simulate from the model
#' pairs(rvinecop(200, vc))
vinecop_dist <- function(pair_copulas, structure,
                         var_types = rep("c", length(pair_copulas[[1]]) + 1)) {
  # create object
  vinecop <- structure(
    list(
      pair_copulas = pair_copulas,
      structure = as_rvine_structure(structure)
    ),
    class = "vinecop_dist"
  )

  # sanity checks
  assert_that(is.list(pair_copulas), correct_var_types(var_types))
  if (length(pair_copulas) > length(pair_copulas[[1]])) {
    stop("'pair_copulas' has more trees than variables.")
  }

  pc_lst <- unlist(pair_copulas, recursive = FALSE)
  if (!all(sapply(pc_lst, function(x) inherits(x, "bicop_dist")))) {
    stop("some objects in pair_copulas aren't of class 'bicop_dist'")
  }

  vinecop$structure <- truncate_model(
    vinecop$structure,
    length(vinecop$pair_copulas)
  )
  vinecop$var_types <- var_types
  vinecop_check_cpp(vinecop)
  vinecop$npars <- sum(sapply(pc_lst, function(x) x[["npars"]]))
  vinecop$loglik <- NA

  vinecop
}
