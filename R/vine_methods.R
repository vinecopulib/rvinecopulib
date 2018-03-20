#' Vine copula distributions
#' 
#' Density, distribution function and random generation 
#' for the vine copula distribution.
#' 
#' @name vine_distributions
#' @aliases dvine pvine rvine dvine_dist pvine_dist rvine_dist
#' @param u evaluation points, either a length d vector or a d-column matrix,
#'   where d is the number of variables in the vine.
#' @param vine an object of class `"vine_dist"`.
#' @details 
#' See [vine] for the estimation and construction of vine copula models. 
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
#' `dvine()` gives the density, `pvine()` gives the distribution function, 
#' and `rvine()` generates random deviates.
#' 
#' The length of the result is determined by `n` for `rvine()`, and 
#' the number of rows in `u` for the other functions.
#' 
#' The `vine` object is recycled to the length of the 
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
#' vc <- vine_dist(pcs, mat)
#' 
#' # simulate from the model
#' u <- rvine(200, vc)
#' pairs(u)
#' 
#' # evaluate the density and cdf
#' dvine(u[1, ], vc)
#' pvine(u[1, ], vc)
#' @rdname vine_methods
#' @export
dvine <- function(u, vine) {
    stopifnot(inherits(vine, "vine_dist"))
    vine_pdf_cpp(if_vec_to_matrix(u), vine)
}

#' @rdname vine_methods
#' @param n_mc number of samples used for quasi Monte Carlo integration.
#' @export
pvine <- function(u, vine, n_mc = 10^4) {
    stopifnot(inherits(vine, "vine_dist"))
    vine_cdf_cpp(if_vec_to_matrix(u), vine, n_mc)
}

#' @rdname vine_methods
#' @param n number of observations.
#' @param U optionally, an \eqn{n \times d} matrix of values in \eqn{(0,1)}.
#'    The result is then the inverse Rosenblatt transform of `U`; if `U` is a
#'    matrix of independent \eqn{U(0, 1)} variables, this simulates data 
#'    from `vine`.
#' @export
rvine <- function(n, vine, U = NULL) {
    stopifnot(inherits(vine, "vine_dist"))
    d <- ncol(vine$matrix)
    U <- prep_uniform_data(n, d, U)
    stopifnot(inherits(vine, "vine_dist"))
    U <- vine_inverse_rosenblatt_cpp(U, vine)
    if (!is.null(vine$names))
        colnames(U) <- vine$names
    
    U
}

#' @export
print.vine_dist <- function(x, ...) {
    d <- nrow(x$matrix)
    cat(d, "-dimensional vine copula model ('vine_dist')", sep = "")
    n_trees <- length(x$pair_copulas)
    if (n_trees < d - 1)
        cat(", ", n_trees, "-truncated", sep = "")
    cat("\n")
}

#' @importFrom utils capture.output
#' @export
summary.vine_dist <- function(object, ...) {
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
            if (pc$family %in% setdiff(family_set_nonparametric, "indep")) {
                pc$parameters <- paste0(round(pc$npars, 2), sep = " d.f.")
            }
            mdf$family[k]     <- pc$family
            mdf$rotation[k]   <- pc$rotation
            mdf$parameters[k] <- list(pc$parameters)
            k <- k + 1
        }
    }
    print.vine_dist(object)
    cat(strrep("-", 63), "\n", sep = "")
    mdf
}

#' Predictions and fitted values for a vine copula model
#'
#' Predictions of the density and distribution function
#' for a vine copula model.
#'
#' @name vine_predict_and_fitted
#' @aliases fitted.vine predict.vine
#' @param object a `vine` object.
#' @param newdata points where the fit shall be evaluated.
#' @param what what to predict, either `"pdf"` or `"cdf"`.
#' @param n_mc number of samples used for quasi Monte Carlo integration when
#'    `what = "cdf"`.
#' @param ... unused.
#'
#' @return 
#' `fitted()` and `predict()` have return values similar to [dvine()] 
#' and [pvine()].
#' @export
#' @rdname predict_vine
#' @examples
#' u <- sapply(1:5, function(i) runif(50))
#' fit <- vine(u, "par")
#' all.equal(predict(fit, u), fitted(fit))
predict.vine <- function(object, newdata, what = "pdf", n_mc = 10^4, ...) {
    stopifnot(what %in% c("pdf", "cdf"))
    newdata <- if_vec_to_matrix(newdata)
    switch(
        what,
        "pdf" = vine_pdf_cpp(newdata, object),
        "cdf" = vine_cdf_cpp(object$data, object, n_mc)
    )
}

#' @rdname predict_vine
#' @export
fitted.vine <- function(object, what = "pdf", n_mc = 10^4, ...) {
    if (is.null(object$data))
        stop("data have not been stored, use keep_data = TRUE when fitting.")
    stopifnot(what %in% c("pdf", "cdf"))
    switch(
        what,
        "pdf" = vine_pdf_cpp(object$data, object),
        "cdf" = vine_cdf_cpp(object$data, object, n_mc)
    )
}

#' @export
logLik.vine <- function(object, ...) {
    if (is.null(object$data))
        stop("data have not been stored, use keep_data = TRUE when fitting.")
    pc_lst <- unlist(object$pair_copulas, recursive = FALSE)
    npars <- ifelse(length(pc_lst) == 0, 0, 
                    sum(sapply(pc_lst, function(x) x[["npars"]])))
    structure(vine_loglik_cpp(object$data, object), "df" = npars)
}

#' @export
print.vine <- function(x, ...) {
    d <- nrow(x$matrix)
    cat(d, "-dimensional vine copula fit ('vine')", sep = "")
    n_trees <- length(x$pair_copulas)
    if (n_trees < d - 1)
        cat(", ", n_trees, "-truncated", sep = "")
    cat("\n")
    cat("nobs =", x$nobs, "  ")
    if (!is.null(x$data)) {
        info <- vine_fit_info(x)
        cat("logLik =", round(info$logLik, 2), "  ")
        cat("npars =", round(info$npars, 2), "  ")
        cat("AIC =", round(info$AIC, 2), "  ")
        cat("BIC =", round(info$BIC, 2), "  ")
        attr(x, "info") <- info
    } else {
        cat("(for mor information, fit model with keep_data = TRUE)")
    }
    cat("\n")
    invisible(x)
}


#' @export
summary.vine <- function(object, ...) {
    info <- attr(print.vine(object), "info")
    capture.output(s <- summary.vine_dist(object))
    cat(strrep("-", 63), "\n", sep = "")
    attr(s, "info") <- info
    s
}


vine_fit_info <- function(vc) {
    stopifnot(inherits(vc, "vine"))
    ll <- logLik(vc)
    list(
        nobs   = vc$nobs,
        logLik = ll[1],
        npars  = attr(ll, "df"),
        AIC    = -2 * ll[1] + 2 * attr(ll, "df"),
        BIC    = -2 * ll[1] + log(vc$nobs) * attr(ll, "df")
    )
}
