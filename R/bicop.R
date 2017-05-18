#' Bivariate copula models
#' 
#' @aliases bicop_fit bicop
#'
#' @param data a matrix or data.frame (copula data should have approximately
#'  uniform margins).
#' @param family_set a character vector of families; as in `bicop_dist()`,
#'   see *Details* for additional options.
#' @param method the estimation method for parametric models, either `"mle"` for
#'   maximum likelihood or `"itau"` for inversion of Kendall's tau (only 
#'   available for one-parametr families and `"t"`.
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
#' \item{"par"}{parametric families.}
#' \item{"nonpar"}{nonparametric families.}
#' \item{"arch"}{archimedean families.}
#' \item{"ell"}{elliptical families.}
#' \item{"bb"}{BB families.}
#' \item{"onepar"}{one parameter families.}
#' \item{"twopar"}{two parameter families.}
#' }
#'
#' @return An object inherting from `bicop_fit` and `bicop_dist`.
#'
#' @examples
#' u <- rbicop(500, "gauss", 0, 0.5)
#' fit1 <- bicop_fit(u, "par")
#' fit1
#' 
#' @export
bicop_fit <- function(data, family_set = "all", method = "mle", mult = 1, 
                      selcrit = "bic", presel = TRUE, keep_data = TRUE) {
    
    # family_set can only use standard family names in cpp
    family_set <- family_set_all_defs[pmatch(family_set, family_set_all_defs)]
    family_set <- expand_family_set(family_set)
    
    ## fit and select copula model
    data <- if_vec_to_matrix(data)
    bicop <- bicop_select_cpp(
        data = data, 
        family_set = family_set,
        method = method,
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
        method = method,
        mult = mult,
        selcrit = selcrit,
        presel = presel
    )
    
    as.bicop_fit(bicop)
}

as.bicop_fit <- function(object) {
    if (!all(c("family", "rotation", "parameters", "npars") %in% names(object)))
        stop("object cannot be coerced to class 'bicop_fit'")
    structure(object, class = c("bicop_fit", "bicop_dist"))
}


#' Predictions and fitted values for a bivariate copula model
#'
#' @param object a `bicop_fit` object.
#' @param newdata points where the fit shall be evaluated.
#' @param what what to predict, one of `"pdf"`, `"hfunc1"`, `"hfunc2"`, 
#'    `"hinv1"`, `"hinv2"`.
#' @param ... unused.
#'
#' @export
#'
#' @examples
#' u <- rbicop(500, "gauss", 0, 0.5)
#' fit <- bicop_fit(u, "par")
#' all.equal(predict(fit, u, "hfunc1"), fitted(fit, "hfunc1"))
predict.bicop_fit <- function(object, newdata, what = "pdf", ...) {
    stopifnot(what %in% c("pdf", "hfunc1", "hfunc2", "hinv1", "hinv2"))
    newdata <- if_vec_to_matrix(newdata)
    switch(
        what,
        "pdf"    = bicop_pdf_cpp(newdata, object),
        "hfunc1" = bicop_hfunc1_cpp(newdata, object),
        "hfunc2" = bicop_hfunc2_cpp(newdata, object),
        "hinv1"  = bicop_hinv1_cpp(newdata, object),
        "hinv2"  = bicop_hinv2_cpp(newdata, object)
    )
}

#' @rdname predict.bicop_fit
#' @export
fitted.bicop_fit <- function(object, what = "pdf", ...) {
    if (is.null(object$data))
        stop("data have not been stored, use keep_data = TRUE when fitting.")
    stopifnot(what %in% c("pdf", "hfunc1", "hfunc2", "hinv1", "hinv2"))
    switch(
        what,
        "pdf"    = bicop_pdf_cpp(object$data, object),
        "hfunc1" = bicop_hfunc1_cpp(object$data, object),
        "hfunc2" = bicop_hfunc2_cpp(object$data, object),
        "hinv1"  = bicop_hinv1_cpp(object$data, object),
        "hinv2"  = bicop_hinv2_cpp(object$data, object)
    )
}

logLik.bicop_fit <- function(object, ...) {
    if (is.null(object$data))
        stop("data have not been stored, use keep_data = TRUE when fitting.")
    structure(bicop_loglik_cpp(object$data, object), "df" = object$npars)
}