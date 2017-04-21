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
#'   available for one-parametr families and '"t"`.
#' @param mult multiplier for the smoothing parameters of nonparametric 
#'   families. Values larger than 1 make the estimate more smooth, values less
#'   than 1 less smooth.
#' @param selcrit criterion for family selection, either "loglik", "aic", or 
#'   "bic".
#' @param presel whether the familyset should be thinned out according to
#'   symmetry characteristics of the data.
#' @param save_data whether the data should be stored (necessary for computing
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
#' @export
#'
#' @examples
#' u <- rbicop(500, "gauss", 0, 0.5)
#' fit1 <- bicop_fit(object$data, "par")
#' fit2 <- bicop(~., as.data.frame(u), "par")
#' all.equal(fit1, fit2)
bicop_fit <- function(data, family_set = "all", method = "mle", mult = 1, 
                      selcrit = "bic",  presel = TRUE, keep_data = TRUE) {
    ## fit and select copula model
    data <- if_vec_to_matrix(data)
    bicop <- bicop_select_cpp(
        data = data, 
        family_set = expand_familyset(family_set),
        method = method,
        mult = mult,
        selcrit = selcrit,
        presel = presel
    )
    
    ## add information about the fit
    if (keep_data)
        bicop$data <- data
    bicop$controls <- list(
        family_set = family_set,
        method = method,
        mult = mult,
        selcrit = selcrit,
        presel = presel
    )
    
    structure(bicop, class = c("bicop_fit", "bicop_dist"))
}

#' @rdname bicop_fit
#' @param formula a formula with empty LHS and the two variables to model on the
#'    RHS, e.g., `~ u1 + u2`, or `~ .` if the data only has two columns.
#'  
#' @export
bicop <- function(formula, data, family_set = "all", method = "mle", mult = 1, 
                  selcrit = "bic",  presel = TRUE, keep_data = TRUE) {
    ## extract data matrix
    if (!inherits(data, "data.frame"))
        data <- as.data.frame(data)
    data <- as.matrix(model.frame(formula, data))
    if (all(colnames(data) == c("V1", "V2")))
        attr(data, "dimnames") <- NULL
    ## call default fit method
    bicop_fit(
        data = data, 
        family_set = family_set,
        method = method, 
        mult = mult, 
        selcrit = selcrit, 
        presel = presel, 
        keep_data = keep_data
    )
}

#' Predictions and fitted values for a bivariate copula model
#'
#' @param object a `bicop_fit` object.
#' @param newdata points where the fit shall be evaluated.
#' @param what what to predict, one of `"pdf"`, `"hfunc1"``, `"hfunc2"``, 
#'    `"hinv1"``, `"hinv2"`.
#' @param ... unused.
#'
#' @export
#'
#' @examples
#' u <- rbicop(500, "gauss", 0, 0.5)
#' fit <- bicop_fit(object$data, "par")
#' all.equal(predict(fit, u, "hfunc1"), fitted(fit, "hfunc1"))
predict.bicop_fit <- function(object, newdata, what = "pdf", ...) {
    stopifnot(what %in% c("pdf", "hfunc1", "hfunc2", "hinv1", "hinv2"))
    u <- if_vec_to_matrix(newdata)
    switch(
        what,
        "pdf"    = bicop_pdf_cpp(object$data, object),
        "hfunc1" = bicop_hfunc1_cpp(object$data, object),
        "hfunc2" = bicop_hfunc2_cpp(object$data, object),
        "hinv1"  = bicop_hinv1_cpp(object$data, object),
        "hinv2"  = bicop_hinv2_cpp(object$data, object)
    )
}

#' @rdname predict.bicop_fit
#' @export
fitted.bicop_fit <- function(object, what = "pdf", ...) {
    if (is.null(object$data))
        stop("data have not been stored, use keep_data = TRUE when fitting.")
    stopifnot(what %in% c("pdf", "hfunc1", "hfunc2", "hinv1", "hinv2"))
    
    u <- if_vec_to_matrix(newdata)
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