#' Vine copula distributions
#' 
#' Density, distribution function and random generation 
#' for the vine copula distribution.
#' 
#' @name vinecop_distributions
#' @aliases dvinecop pvinecop rvinecop dvinecop_dist pvinecop_dist rvinecop_dist
#' @param u evaluation points, either a length d vector or a d-column matrix,
#'   where d is the number of variables in the vine.
#' @param vinecop an object of class `"vinecop_dist"`.
#' @details 
#' See [vinecop] for the estimation and construction of vine copula models. 
#' Here, the density, distribution function and random generation 
#' for the vine copulas are standard.
#' 
#' The Rosenblatt transform (Rosenblatt, 1952) \eqn{U = T(V)} of a random vector
#' \eqn{V = (V_1,\ldots,V_d) ~ C} is defined as
#' \deqn{ 
#'   U_1 = V_1, U_2 = C(V_2|V_1), \ldots, U_d =C(V_d|V_1,\ldots,V_{d-1}), 
#' } 
#' where \eqn{C(v_k|v_1,\ldots,v_{k-1})} is the conditional distribution of 
#' \eqn{V_k} given \eqn{V_1 \ldots, V_{k-1}, k = 2,\ldots,d}. The vector \eqn{V}
#' are then independent standard uniform variables. The inverse operation 
#' \deqn{ 
#'   V_1 = U_1, V_2 = C^{-1}(U_2|U_1), \ldots, V_d =C^{-1}(U_d|U_1,\ldots,U_{d-1}), 
#' } 
#' can can be used to simulate from a copula. For any copula \eqn{C}, if 
#' \eqn{U} is a vector of independent random variables, \eqn{V = T^{-1}(U)} has 
#' distribution \eqn{C}.
#' @return 
#' `dvinecop()` gives the density, `pvinecop()` gives the distribution function, 
#' and `rvinecop()` generates random deviates.
#' 
#' The length of the result is determined by `n` for `rvinecop()`, and 
#' the number of rows in `u` for the other functions.
#' 
#' The `vinecop` object is recycled to the length of the 
#' result.
#' @examples
#' # specify pair-copulas
#' bicop <- bicop_dist("bb1", 90, c(3, 2))
#' pcs <- list(
#'     list(bicop, bicop),  # pair-copulas in first tree 
#'     list(bicop)          # pair-copulas in second tree 
#'  )
#'  
#' # specify R-vine matrix
#' mat <- matrix(c(1, 2, 3, 1, 2, 0, 1, 0, 0), 3, 3) 
#' 
#' # set up vine copula model
#' vc <- vinecop_dist(pcs, mat)
#' 
#' # simulate from the model
#' u <- rvinecop(200, vc)
#' pairs(u)
#' 
#' # evaluate the density and cdf
#' dvinecop(u[1, ], vc)
#' pvinecop(u[1, ], vc)
#' @rdname vinecop_methods
#' @export
dvinecop <- function(u, vinecop) {
    stopifnot(inherits(vinecop, "vinecop_dist"))
    vinecop_pdf_cpp(if_vec_to_matrix(u), vinecop)
}

#' @rdname vinecop_methods
#' @param n_mc number of samples used for quasi Monte Carlo integration.
#' @export
pvinecop <- function(u, vinecop, n_mc = 10^4) {
    stopifnot(inherits(vinecop, "vinecop_dist"))
    vinecop_cdf_cpp(if_vec_to_matrix(u), vinecop, n_mc)
}

#' @rdname vinecop_methods
#' @param n number of observations.
#' @param U optionally, an \eqn{n \times d} matrix of values in \eqn{(0,1)}.
#'    The result is then the inverse Rosenblatt transform of `U`; if `U` is a
#'    matrix of independent \eqn{U(0, 1)} variables, this simulates data 
#'    from `vinecop`.
#' @export
rvinecop <- function(n, vinecop, U = NULL) {
    stopifnot(inherits(vinecop, "vinecop_dist"))
    d <- ncol(vinecop$matrix)
    U <- prep_uniform_data(n, d, U)
    stopifnot(inherits(vinecop, "vinecop_dist"))
    U <- vinecop_inverse_rosenblatt_cpp(U, vinecop)
    if (!is.null(vinecop$names))
        colnames(U) <- vinecop$names
    
    U
}

#' @export
print.vinecop_dist <- function(x, ...) {
    d <- nrow(x$matrix)
    cat(d, "-dimensional vine copula model", sep = "")
    n_trees <- length(x$pair_copulas)
    if (n_trees < d - 1)
        cat(" (", n_trees, "-truncated)", sep = "")
    cat("\n")
}

#' @importFrom utils capture.output
#' @export
summary.vinecop_dist <- function(object, ...) {
    mat <- object$matrix
    d <- nrow(mat)
    n_trees <- length(object$pair_copulas)
    n_pcs <- length(unlist(object$pair_copulas, recursive = FALSE))
    mdf <- as.data.frame(matrix(NA, n_pcs, 7))
    names(mdf) <- c("tree", "edge", 
                    "conditioned", "conditioning", 
                    "family", "rotation", "parameters")
    k <- 1
    for (t in seq_len(n_trees)) {
        for (e in seq_len(d - t)) {
            mdf$tree[k] <- t
            mdf$edge[k] <- e
            mdf$conditioned[k]  <- list(c(mat[d - e + 1, e], mat[t, e]))
            mdf$conditioning[k] <- list(mat[rev(seq_len(t - 1)), e])
            pc <- object$pair_copulas[[t]][[e]]
            mdf$family[k]     <- pc$family
            mdf$rotation[k]   <- pc$rotation
            mdf$parameters[k] <- list(pc$parameters)
            k <- k + 1
        }
    }
    pr <- capture.output(print(object))
    cat(pr,  "\n", strrep("-", nchar(pr)), "\n", sep = "")
    print(mdf, digits = 2, row.names = FALSE)

    invisible(mdf)
}

#' Predictions and fitted values for a vine copula model
#'
#' Predictions of the density and distribution function
#' for a vine copula model.
#'
#' @name vinecop_predict_and_fitted
#' @aliases fitted.vinecop predict.vinecop
#' @param object a `vinecop` object.
#' @param newdata points where the fit shall be evaluated.
#' @param what what to predict, either `"pdf"` or `"cdf"`.
#' @param n_mc number of samples used for quasi Monte Carlo integration when
#'    `what = "cdf"`.
#' @param ... unused.
#'
#' @return 
#' `fitted()` and `predict()` have return values similar to [dvinecop()] 
#' and [pvinecop()].
#' @export
#' @rdname predict_vinecop
#' @examples
#' u <- sapply(1:5, function(i) runif(50))
#' fit <- vinecop(u, "par")
#' all.equal(predict(fit, u), fitted(fit))
predict.vinecop <- function(object, newdata, what = "pdf", n_mc = 10^4, ...) {
    stopifnot(what %in% c("pdf", "cdf"))
    newdata <- if_vec_to_matrix(newdata)
    switch(
        what,
        "pdf" = vinecop_pdf_cpp(newdata, object),
        "cdf" = vinecop_cdf_cpp(object$data, object, n_mc)
    )
}

#' @rdname predict_vinecop
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


#' @export
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
