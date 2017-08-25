#' Fit vine copula models
#' 
#' Automated fitting for vine copula models
#' 
#' @inheritParams bicop
#' @param matrix an R-vine matrix specifying the structure matrix (see 
#'   [`check_rvine_matrix()`]), or `NA` for
#'   automatic structure selection (default).
#' @param trunc_lvl the truncation level of the vine copula; `Inf` means no
#'   truncation, `NA` indicates that the truncation level should be selected
#'   automatically by GIC.
#' @param tree_crit the criterion for tree selection, one of `"tau"`, `"rho"`,
#'    `"hoeffd"` for Kendall's \eqn{tau}, Spearman's \eqn{rho}, and Hoeffding's
#'    \eqn{D}, respectively.
#' @param threshold for thresholded vine copulas; `NA` indicates that the 
#'   threshold should be selected automatically by GIC.
#' @param show_trace logical; whether a trace of the fitting progress should be 
#'    printed.
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
#' @return An object inherting from `bicop` and `bicop_dist`.
#'
#' @examples
#' u <- sapply(1:3, function(i) runif(50))
#' fit <- vinecop(u, "par")
#' fit
#' summary(fit)
#' str(fit, 3)
#' 
#' @export
vinecop <- function(data, family_set = "all", matrix = NA, 
                    par_method = "mle", nonpar_method = "constant",
                    mult = 1, selcrit = "bic", presel = TRUE, 
                    trunc_lvl = Inf, tree_crit = "tau", threshold = 0, 
                    keep_data = TRUE, show_trace = FALSE) {
    
    ## family_set can only use standard family names in cpp
    family_set <- family_set_all_defs[pmatch(family_set, family_set_all_defs)]
    family_set <- expand_family_set(family_set)
    
    ## pre-process input
    data <- if_vec_to_matrix(data)
    if (any(is.na(matrix)))
        matrix <- as.matrix(0)
    
    ## fit and select copula model
    vinecop <- vinecop_select_cpp(
        data = data, 
        matrix = matrix,
        family_set = family_set,
        par_method = par_method,
        nonpar_method = nonpar_method,
        mult = mult,
        selection_criterion = selcrit,
        preselect_families = presel,
        truncation_level = ifelse(  # Inf cannot be passed to C++
            is.finite(trunc_lvl),
            trunc_lvl, 
            .Machine$integer.max
        ),
        tree_criterion = tree_crit,
        threshold = threshold,
        select_truncation_level = is.na(trunc_lvl),
        select_threshold = is.na(threshold),
        show_trace = show_trace
    )
    
    ## make all pair-copulas bicop objects
    vinecop$pair_copulas <- lapply(
        vinecop$pair_copulas, 
        function(tree) lapply(tree, as.bicop)
    )
    
    ## add information about the fit
    if (keep_data) {
        vinecop$data <- data
    }
    vinecop$controls <- list(
        family_set = family_set,
        par_method = par_method,
        nonpar_method = nonpar_method,
        mult = mult,
        selcrit = selcrit,
        presel = presel,
        trunc_lvl = trunc_lvl,
        tree_crit = tree_crit,
        threshold = threshold
    )
    vinecop$nobs <- nrow(data)
    
    structure(vinecop, class = c("vinecop", "vinecop_dist"))
}

#' Predictions and fitted values for a vine copula model
#'
#' @param object a `vinecop` object.
#' @param newdata points where the fit shall be evaluated.
#' @param what what to predict, either `"pdf"` or `"cdf"`.
#' @param n_mc number of samples used for quasi Monte Carlo integration when
#'    `what = "cdf"`.
#' @param ... unused.
#'
#' @export
#'
#' @examples
#' u <- sapply(1:5, function(i) runif(50))
#' fit <- vinecop(u, "par")
#' all.equal(predict(fit, u), fitted(fit))
#' logLik(fit)
predict.vinecop <- function(object, newdata, what = "pdf", n_mc = 10^4, ...) {
    stopifnot(what %in% c("pdf", "cdf"))
    newdata <- if_vec_to_matrix(newdata)
    switch(
        what,
        "pdf" = vinecop_pdf_cpp(newdata, object),
        "cdf" = vinecop_cdf_cpp(object$data, object, n_mc)
    )
}

#' @rdname predict.vinecop
#' @export
fitted.vinecop <- function(object, what = "pdf", n_mc = 10^4, ...) {
    if (is.null(object$data))
        stop("data have not been stored, use keep_data = TRUE when fitting.")
    stopifnot(what %in% c("pdf", "cdf"))
    switch(
        what,
        "pdf" = vinecop_pdf_cpp(object$data, object),
        "cdf" = vinecop_cdf_cpp(object$data, object, n_mc)
    )
}

#' @rdname predict.vinecop
#' @export
logLik.vinecop <- function(object, ...) {
    if (is.null(object$data))
        stop("data have not been stored, use keep_data = TRUE when fitting.")
    pc_lst <- unlist(object$pair_copulas, recursive = FALSE)
    npars <- sum(sapply(pc_lst, function(x) x[["npars"]]))
    structure(vinecop_loglik_cpp(object$data, object), "df" = npars)
}

#' @export
print.vinecop <- function(x, ...) {
    info <- vinecop_fit_info(x)
    print.vinecop_dist(x)
    cat(" fit\n")
    cat("nobs =", info$nobs, "  ")
    cat("logLik =", round(info$logLik, 2), "  ")
    cat("npars =", round(info$npars, 2), "  ")
    cat("AIC =", round(info$AIC, 2), "  ")
    cat("BIC =", round(info$BIC, 2), "  ")
    attr(x, "info") <- info
    invisible(x)
}

summary.vinecop <- function(object, ...) {
    info <- attr(print.vinecop(object), "info")
    cat("\n----\n")
    s <- summary.vinecop_dist(object)
    attr(s, "info") <- info
    invisible(s)
}


vinecop_fit_info <- function(vc) {
    stopifnot(inherits(vc, "vinecop"))
    ll <- logLik(vc)
    list(
        nobs   = vc$nobs,
        logLik = ll[1],
        npars  = attr(ll, "df"),
        AIC    = -2 * ll[1] + 2 * attr(ll, "df"),
        BIC    = -2 * ll[1] + log(vc$nobs) * attr(ll, "df")
    )
}
