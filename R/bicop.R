#' Bivariate copula models
#' 
#' @aliases bicop_dist
#'
#' @param data a matrix or data.frame (copula data should have approximately
#'  uniform margins).
#' @param family_set a character vector of families; see *Details* for 
#' additional options.
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
#'   `"bic"`, `"mbic"`. For `vinecop()` there is the additional option `"mbicv"`.
#' @param weights optional vector of weights for each observation.
#' @param psi0 see [mBICV()].
#' @param presel whether the family set should be thinned out according to
#'   symmetry characteristics of the data.
#' @param keep_data whether the data should be stored (necessary for using 
#'   [fitted()]).
#' @param cores number of cores to use; if more than 1, estimation for multiple
#'   families is done in parallel.
#' 
#' @details
#' 
#' The implemented families are:\cr
#' 
#' `"indep"` = Independence copula.\cr
#' `"gaussian"` = Gaussian copula.\cr
#' `"t"` = Student t copula.\cr
#' `"clayton"` = Clayton copula.\cr
#' `"gumbel"` = Gumbel copula.\cr
#' `"frank"` = Frank copula.\cr
#' `"joe"` = Joe copula.\cr
#' `"bb1"` = BB1 copula.\cr
#' `"bb6"` = BB6 copula.\cr
#' `"bb7"` = BB7 copula.\cr
#' `"bb8"` = BB8 copula.\cr
#' `"tll"` = transformation kernel local likelihood, only for `bicop()`.\cr
#' 
#' In addition, the following convenience definitions can be used (and combined) 
#' with `bicop`:\cr
#' 
#' `"all"` =  all families.\cr
#' `"parametric"` =  parametric families.\cr
#' `"nonparametric"` =  nonparametric families.\cr
#' `"archimedean"` =  archimedean families.\cr
#' `"elliptical"` =  elliptical families.\cr
#' `"bbs"` =  BB families.\cr
#' `"oneparametric"` =  one parameter families.\cr
#' `"twoparametric"` =  two parameter families.\cr
#' `"itau"` =  one parameter families and Student t copula.\cr
#' Partial matching is activated. For example, `"gauss"` is equivalent to 
#' `"gaussian"`, or you can write  `"nonpar"` instead of `"nonparametric"`.
#'
#'
#' @return Objects inheriting from `bicop_dist` for `bicop_dist()`, and
#' `bicop` and `bicop_dist` for `bicop()`.
#' 
#' Object from the `bicop_dist` class are lists containing:
#' 
#' * `family`, a `character` indicating the copula family.
#' * `rotation`, an `integer` indicating the rotation (i.e., either 0, 90, 180, 
#' or 270).
#' * `parameters`, a `numeric` vector or matrix of parameters.
#' * `npars`, a `numeric` with the (effective) number of parameters.
#' 
#' Additionally, objects from the `bicop` class contain:
#' 
#' * `data` (optionally, if `keep_data = TRUE` was used), the dataset that was 
#' passed to [bicop()].
#' * `controls`, a `list` with the set of fit controls that was passed to [bicop()].
#' * `nobs`, an `integer` with the number of observations that was used 
#' to fit the model.
#' 
#' @examples
#' ## bicop_dist objects
#' bicop_dist("gaussian", 0, 0.5)
#' str(bicop_dist("gauss", 0, 0.5))
#' bicop <- bicop_dist("clayton", 90, 3)
#' 
#' ## bicop objects
#' u <- rbicop(500, "gauss", 0, 0.5)
#' fit1 <- bicop(u, "par")
#' fit1
#' 
#' @export
bicop <- function(data, family_set = "all", par_method = "mle",
                  nonpar_method = "quadratic", mult = 1, selcrit = "bic", 
                  weights = numeric(), psi0 = 0.9, presel = TRUE, 
                  keep_data = FALSE, cores = 1) {
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
        is.number(cores), cores > 0
    )
    
    stopifnot(ncol(data) == 2)
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
        num_threads = cores
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
    if (!all(c("family", "rotation", "parameters", "npars") %in% names(object)))
        stop("object cannot be coerced to class 'bicop'")
    structure(object, class = c("bicop", "bicop_dist"))
}

#' @param family the copula family, a string containing the family name (see
#' *Details* for all possible families).
#' @param rotation the rotation of the copula, one of `0`, `90`, `180`, `270`.
#' @param parameters a vector or matrix of copula parameters.
#' @rdname bicop
#' @export
bicop_dist <- function(family = "indep", rotation = 0, parameters = numeric(0)) {
    assert_that(is.string(family), is.number(rotation), is.numeric(parameters))
    if (family %in% setdiff(family_set_nonparametric, "indep"))
        stop("bicop_dist should not be used directly with nonparametric families.")
    
    family <- family_set_all[pmatch(family, family_set_all)]
    dist <- list(family     = family,
                 rotation   = rotation,
                 parameters = as.matrix(parameters),
                 npars      = length(parameters))
    bicop_check_cpp(dist)
    structure(dist, class = "bicop_dist")
}