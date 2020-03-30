#' Fit and select bivariate copula models
#'
#' Fit a bivariate copula model for continuous or discrete data. The family
#' can be selected automatically from a vector of options.
#'
#' @param data a matrix or data.frame with at least two columns, containing the
#'   (pseudo-)observations for the two variables (copula data should have
#'   approximately uniform margins). More columns are required for discrete
#'   models, see *Details*.
#' @param family_set a character vector of families; see *Details* for
#'   additional options.
#' @param par_method the estimation method for parametric models, either `"mle"`
#'   for maximum likelihood or `"itau"` for inversion of Kendall's tau (only
#'   available for one-parameter families and `"t"`.
#' @param nonpar_method the estimation method for nonparametric models, either
#'   `"constant"` for the standard transformation estimator, or
#'   `"linear"`/`"quadratic"` for the local-likelihood approximations of order
#'   one/two.
#' @param mult multiplier for the smoothing parameters of nonparametric
#'   families. Values larger than 1 make the estimate more smooth, values less
#'   than 1 less smooth.
#' @param selcrit criterion for family selection, either `"loglik"`, `"aic"`,
#'   `"bic"`, `"mbic"`. For `vinecop()` there is the additional option
#'   `"mbicv"`.
#' @param weights optional vector of weights for each observation.
#' @param psi0 see [mBICV()].
#' @param presel whether the family set should be thinned out according to
#'   symmetry characteristics of the data.
#' @param keep_data whether the data should be stored (necessary for using
#'   [fitted()]).
#' @param cores number of cores to use; if more than 1, estimation for multiple
#'   families is done in parallel.
#' @param var_types variable types, a length 2 vector; e.g., `c("c", "c")` for
#'   both continuous (default), or `c("c", "d")` for first variable continuous
#'   and second discrete.
#'
#' @details
#'
#' ## Discrete variables
#'
#' When at least one variable is discrete, more than two columns are required
#' for `data`: the first \eqn{n \times 2} block contains realizations of
#' \eqn{F_{X_1}(x_1), F_{X_2}(x_2)}. The second \eqn{n \times 2} block contains
#' realizations of \eqn{F_{X_1}(x_1^-), F_{X_1}(x_1^-)}. The minus indicates a
#' left-sided limit of the cdf. For, e.g., an integer-valued variable, it holds
#' \eqn{F_{X_1}(x_1^-) = F_{X_1}(x_1 - 1)}. For continuous variables the left
#' limit and the cdf itself coincide. Respective columns can be omitted in the
#' second block.
#'
#' ## Family collections
#'
#' The `family_set` argument accepts all families in `bicop_dist()` plus the
#' following convenience definitions:
#'
#' * `"all"` contains all the families,
#'
#' * `"parametric"` contains the parametric families (all except `"tll"`),
#'
#' * `"nonparametric"` contains the nonparametric families (`"indep"` and
#' `"tll"`)
#'
#' * `"onepar"` contains the parametric families with a single parameter,
#'
#' (`"gaussian"`, `"clayton"`, `"gumbel"`, `"frank"`, and `"joe"`),
#'
#' * `"twopar"` contains the parametric families with two parameters
#' (`"student"`, `"bb1"`, `"bb6"`, `"bb7"`, and `"bb8"`),
#'
#' * `"elliptical"` contains the elliptical families,
#'
#' * `"archimedean"` contains the archimedean families,
#'
#' * `"BB"` contains the BB families,
#'
#' * `"itau"` families for which estimation by Kendall's tau inversion is
#' available (`"indep"`,`"gaussian"`, `"student"`,`"clayton"`, `"gumbel"`,
#' `"frank"`, `"joe"`).
#'
#' @return
#' An object inheriting from classes `bicop` and  `bicop_dist` . In addition to
#' the entries contained in `bicop_dist()`, objects from the `bicop` class
#' contain:
#'
#' * `data` (optionally, if `keep_data = TRUE` was used), the dataset that was
#' passed to [bicop()].
#'
#' * `controls`, a `list` with the set of fit controls that was passed to
#' [bicop()].
#'
#' * `loglik` the log-likelihood.
#'
#' * `nobs`, an `integer` with the number of observations that was used to fit
#' the model.
#'
#' @seealso [bicop_dist()], [plot.bicop()], [contour.bicop()], [dbicop()],
#'   [pbicop()], [hbicop()], [rbicop()]
#'
#' @examples
#' ## fitting a continuous model from simulated data
#' u <- rbicop(100, "clayton", 90, 3)
#' fit <- bicop(u, family_set = "par")
#' summary(fit)
#'
#' ## compare fit with true model
#' contour(fit)
#' contour(bicop_dist("clayton", 90, 3), col = 2, add = TRUE)
#'
#' ## fit a model from discrete data
#' x_disc <- qpois(u, 1)  # transform to Poisson margins
#' plot(x_disc)
#' udisc <- cbind(ppois(x_disc, 1), ppois(x_disc - 1, 1))
#' fit_disc <- bicop(udisc, var_types = c("d", "d"))
#' summary(fit_disc)
#' @export
bicop <- function(data, var_types = c("c", "c"), family_set = "all",
                  par_method = "mle", nonpar_method = "quadratic", mult = 1,
                  selcrit = "bic", weights = numeric(), psi0 = 0.9,
                  presel = TRUE, keep_data = FALSE, cores = 1) {
  assert_that(
    is.character(family_set),
    is.string(par_method),
    is.string(nonpar_method),
    is.number(mult), mult > 0,
    is.string(selcrit),
    is.numeric(weights),
    is.number(psi0), psi0 > 0, psi0 < 1,
    is.flag(presel),
    is.flag(keep_data),
    is.number(cores), cores > 0,
    correct_var_types(var_types)
  )

  # check if families known (w/ partial matching) and expand convenience defs
  family_set <- process_family_set(family_set, par_method)

  ## fit and select copula model
  data <- if_vec_to_matrix(data)
  bicop <- bicop_select_cpp(
    data = data,
    family_set = family_set,
    par_method = par_method,
    nonpar_method = nonpar_method,
    mult = mult,
    selcrit = selcrit,
    weights = weights,
    psi0 = psi0,
    presel = presel,
    num_threads = cores,
    var_types = var_types
  )

  ## add information about the fit
  bicop$names <- colnames(data)
  if (keep_data) {
    bicop$data <- data
  }
  bicop$controls <- list(
    family_set = family_set,
    par_method = par_method,
    nonpar_method = nonpar_method,
    mult = mult,
    selcrit = selcrit,
    weights = weights,
    psi0 = psi0,
    presel = presel
  )
  bicop$nobs <- nrow(data)

  as.bicop(bicop)
}

as.bicop <- function(object) {
  if (!all(c("family", "rotation", "parameters", "npars") %in% names(object))) {
    stop("object cannot be coerced to class 'bicop'")
  }
  structure(object, class = c("bicop", "bicop_dist"))
}

#' Bivariate copula models
#'
#' Create custom bivariate copula models by specifying the family, rotation,
#' parameters, and variable types.
#'
#' ## Implemented families
#'
#' | type          | name                  | name in R     |
#' |---------------|-----------------------|---------------|
#' | -             | Independence          | "indep"       |
#' | Elliptical    | Gaussian              | "gaussian"    |
#' | "             | Student t             | "student"     |
#' | Archimedean   | Clayton               | "clayton"     |
#' | "             | Gumbel                | "gumbel"      |
#' | "             | Frank                 | "frank"       |
#' | "             | Joe                   | "joe"         |
#' | "             | Clayton-Gumbel (BB1)  | "bb1"         |
#' | "             | Joe-Gumbel (BB6)      | "bb6"         |
#' | "             | Joe-Clayton (BB7)     | "bb7"         |
#' | "             | Joe-Frank (BB8)       | "bb8"         |
#' | Nonparametric | Transformation kernel | "tll"         |
#'
#' @return
#'
#' An object of class `bicop_dist`, i.e., a list containing:
#'
#' * `family`, a `character` indicating the copula family.
#'
#' * `rotation`, an `integer` indicating the rotation (i.e., either 0, 90, 180,
#' or 270).
#'
#' * `parameters`, a `numeric` vector or matrix of parameters.
#'
#' * `npars`, a `numeric` with the (effective) number of parameters.
#'
#' * `var_types`, the variable types.
#'
#' @param family the copula family, a string containing the family name (see
#' *Details* for all possible families).
#' @param rotation the rotation of the copula, one of `0`, `90`, `180`, `270`.
#' @param parameters a vector or matrix of copula parameters.
#' @param var_types variable types, a length 2 vector; e.g., `c("c", "c")` for
#'   both continuous (default), or `c("c", "d")` for first variable continuous
#'   and second discrete.
#'
#' @seealso [bicop_dist()], [plot.bicop()], [contour.bicop()], [dbicop()],
#'   [pbicop()], [hbicop()], [rbicop()]
#' @export
#' @examples
#' ## Clayton 90° copula with parameter 3
#' cop <- bicop_dist("clayton", 90, 3)
#' cop
#' str(cop)
#'
#' ## visualization
#' plot(cop)
#' contour(cop)
#' plot(rbicop(200, cop))
#'
#' ## BB8 copula model for discrete data
#' cop_disc <- bicop_dist("bb8", 0, c(2, 0.5), var_types = c("d", "d"))
#' cop_disc
#'
bicop_dist <- function(family = "indep", rotation = 0, parameters = numeric(0),
                       var_types = c("c", "c")) {
  assert_that(is.string(family), is.number(rotation), is.numeric(parameters))
  assert_that(correct_var_types(var_types))
  if (family %in% setdiff(family_set_nonparametric, "indep")) {
    stop("bicop_dist should not be used directly with nonparametric families.")
  }

  family <- family_set_all[pmatch(family, family_set_all)]
  dist <- list(
    family = family,
    rotation = rotation,
    parameters = as.matrix(parameters),
    var_types = var_types,
    npars = length(parameters)
  )
  bicop_check_cpp(dist)
  structure(dist, class = "bicop_dist")
}
