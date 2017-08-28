#' Bivariate copula models
#' 
#' @aliases bicop bicop
#'
#' @param data a matrix or data.frame (copula data should have approximately
#'  uniform margins).
#' @param family_set a character vector of families; as in `bicop_dist()`,
#'   see *Details* for additional options.
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
#' @param selcrit criterion for family selection, either `"loglik"`, `"aic"`, or 
#'   `"bic"`.
#' @param presel whether the family set should be thinned out according to
#'   symmetry characteristics of the data.
#' @param keep_data whether the data should be stored (necessary for computing
#'   fit statistics and using `fitted()`).
#' 
#' @details
#' In addition, to the families in `bicop_dist()`, the following convenience
#' defintion can be used (and combined):
#' \describe{
#' \item{"all"}{all families}.
#' \item{"parametric"}{parametric families.}
#' \item{"nonparametric"}{nonparametric families.}
#' \item{"archimedean"}{archimedean families.}
#' \item{"elliptical"}{elliptical families.}
#' \item{"bbs"}{BB families.}
#' \item{"oneparametric "}{one parameter families.}
#' \item{"twoparametric "}{two parameter families.}
#' }
#' Partial matching is activated. For example, you can write `"nonpar"` instead 
#' of the full name `"nonparametric"`.
#'
#'
#' @return An object inherting from `bicop` and `bicop_dist`.
#'
#' @examples
#' u <- rbicop(500, "gauss", 0, 0.5)
#' fit1 <- bicop(u, "par")
#' fit1
#' 
#' @export
bicop <- function(data, family_set = "all", par_method = "mle",
                  nonpar_method = "quadratic", mult = 1, selcrit = "bic", 
                  presel = TRUE, keep_data = TRUE) {
    stopifnot(ncol(data) == 2)
    # family_set can only use standard family names in cpp
    family_set <- family_set_all_defs[pmatch(family_set, family_set_all_defs)]
    family_set <- expand_family_set(family_set)
    
    ## fit and select copula model
    data <- if_vec_to_matrix(data)
    bicop <- bicop_select_cpp(
        data = data, 
        family_set = family_set,
        par_method = par_method,
        nonpar_method = nonpar_method,
        mult = mult,
        selcrit = selcrit,
        presel = presel
    )
    
    ## add information about the fit
    if (keep_data) {
        bicop$data <- data
    }
    bicop$controls <- list(
        family_set = family_set,
        par_method = par_method,
        nonpar_method = nonpar_method,
        mult = mult,
        selcrit = selcrit,
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


#' Predictions and fitted values for a bivariate copula model
#'
#' @param object a `bicop` object.
#' @param newdata points where the fit shall be evaluated.
#' @param what what to predict, one of `"pdf"`, `"cdf"`, `"hfunc1"`, `"hfunc2"`, 
#'    `"hinv1"`, `"hinv2"`.
#' @param ... unused.
#'
#' @export
#'
#' @examples
#' u <- rbicop(500, "gauss", 0, 0.5)
#' fit <- bicop(u, "par")
#' all.equal(predict(fit, u, "hfunc1"), fitted(fit, "hfunc1"))
predict.bicop <- function(object, newdata, what = "pdf", ...) {
    stopifnot(what %in% c("pdf", "cdf", "hfunc1", "hfunc2", "hinv1", "hinv2"))
    newdata <- if_vec_to_matrix(newdata)
    switch(
        what,
        "pdf"    = bicop_pdf_cpp(newdata, object),
        "cdf"    = bicop_cdf_cpp(newdata, object),
        "hfunc1" = bicop_hfunc1_cpp(newdata, object),
        "hfunc2" = bicop_hfunc2_cpp(newdata, object),
        "hinv1"  = bicop_hinv1_cpp(newdata, object),
        "hinv2"  = bicop_hinv2_cpp(newdata, object)
    )
}

#' @rdname predict.bicop
#' @export
fitted.bicop <- function(object, what = "pdf", ...) {
    if (is.null(object$data))
        stop("data have not been stored, use keep_data = TRUE when fitting.")
    stopifnot(what %in% c("pdf", "cdf", "hfunc1", "hfunc2", "hinv1", "hinv2"))
    switch(
        what,
        "pdf"    = bicop_pdf_cpp(object$data, object),
        "cdf"    = bicop_cdf_cpp(object$data, object),
        "hfunc1" = bicop_hfunc1_cpp(object$data, object),
        "hfunc2" = bicop_hfunc2_cpp(object$data, object),
        "hinv1"  = bicop_hinv1_cpp(object$data, object),
        "hinv2"  = bicop_hinv2_cpp(object$data, object)
    )
}

#' @importFrom stats logLik
logLik.bicop <- function(object, ...) {
    if (is.null(object$data))
        stop("data have not been stored, use keep_data = TRUE when fitting.")
    structure(bicop_loglik_cpp(object$data, object), "df" = object$npars)
}

print.bicop <- function(x, ...) {
    info <- bicop_fit_info(x)
    if (x$family %in% setdiff(family_set_nonparametric, "indep")) {
        x$parameters <- "[30x30 grid]"
    }
    cat("Bivariate copula fit ('bicop'): ",
        "family = ", x$family,
        ", rotation = ", x$rotation,
        ", parameters = ", x$parameters,
        "\n",
        sep = "")
    cat("nobs =", info$nobs, "  ")
    cat("logLik =", round(info$logLik, 2), "  ")
    cat("npars =", round(info$npars, 2), "  ")
    cat("AIC =", round(info$AIC, 2), "  ")
    cat("BIC =", round(info$BIC, 2), "  ")
    
    attr(x, "info") <- info
    invisible(x)
}

bicop_fit_info <- function(bc) {
    ll <- logLik(bc)
    list(
        nobs   = bc$nobs,
        logLik = ll[1],
        npars  = attr(ll, "df"),
        AIC    = -2 * ll[1] + 2 * attr(ll, "df"),
        BIC    = -2 * ll[1] + log(bc$nobs) * attr(ll, "df")
    )
}