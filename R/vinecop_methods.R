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
#' 
#' # get single pair copula
#' get_pair_copula(vc, 1, 1)
#' 
#' # get all pair copulas
#' get_all_pair_copulas(vc)
#' 
#' # get vine matrix
#' get_matrix(vc)
#' 
#' # extract a truncated sub-vine based on truncation level supplied by user
#' truncate_model(vc, 1)
#' 
#' @rdname vinecop_methods
#' @export
dvinecop <- function(u, vinecop) {
    assert_that(inherits(vinecop, "vinecop_dist"))
    vinecop_pdf_cpp(if_vec_to_matrix(u), vinecop)
}

#' @rdname vinecop_methods
#' @param n_mc number of samples used for quasi Monte Carlo integration.
#' @export
pvinecop <- function(u, vinecop, n_mc = 10^4) {
    assert_that(inherits(vinecop, "vinecop_dist"), is.number(n_mc))
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
    assert_that(inherits(vinecop, "vinecop_dist"))
    d <- ncol(vinecop$matrix)
    U <- prep_uniform_data(n, d, U)
    U <- vinecop_inverse_rosenblatt_cpp(U, vinecop)
    if (!is.null(vinecop$names))
        colnames(U) <- vinecop$names
    
    U
}

#' @export
print.vinecop_dist <- function(x, ...) {
    d <- nrow(x$matrix)
    cat(d, "-dimensional vine copula model ('vinecop_dist')", sep = "")
    n_trees <- length(x$pair_copulas)
    if (n_trees < d - 1)
        cat(", ", n_trees, "-truncated", sep = "")
    cat("\n")
    invisible(x)
}

#' @importFrom utils capture.output
#' @export
summary.vinecop_dist <- function(object, ...) {
    mat <- object$matrix
    d <- nrow(mat)
    n_trees <- length(object$pair_copulas)
    n_pcs <- length(unlist(object$pair_copulas, recursive = FALSE))
    mdf <- as.data.frame(matrix(NA, n_pcs, 9))
    names(mdf) <- c("tree", "edge", 
                    "conditioned", "conditioning", 
                    "family", "rotation", "parameters", "df", "tau")
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
            mdf$df[k] <- pc$npars
            mdf$tau[k] <- par_to_ktau(pc)
            k <- k + 1
        }
    }
    class(mdf) <- c("vinecop_dist_summary", class(mdf))
    mdf
}

#' @export
print.vinecop_dist_summary <- function(x, ...) {
    x_print <- x[1:min(nrow(x), 10), ]
    x_print[x_print$family == "tll", "parameters"] <- list("[30x30 grid]")
    cat("# A data.frame:", nrow(x), "x", ncol(x), "\n")
    print.data.frame(x_print, digits = 2)
    if (nrow(x) > 10)
        cat("# ... with", nrow(x) - 10, "more rows\n")
    invisible(x)
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
    assert_that(in_set(what, c("pdf", "cdf")), is.number(n_mc))
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
    assert_that(in_set(what, c("pdf", "cdf")), is.number(n_mc))
    switch(
        what,
        "pdf" = vinecop_pdf_cpp(object$data, object),
        "cdf" = vinecop_cdf_cpp(object$data, object, n_mc)
    )
}

#' @export
logLik.vinecop <- function(object, ...) {
    structure(object$loglik, "df" = object$npars)
}

#' calculates the vine copula Bayesian information criterion (vBIC), which is 
#' defined as
#' \deqn{\mathrm{BIC} = -2\, \mathrm{loglik} +  \nu \ln(n), - 2 * 
#' \sum_{t=1}^(d - 1) \{q_t log(\psi_0^t) - (d - t - q_t) log(1 - \psi_0^t)\}
#' }
#' where \eqn{\mathrm{loglik}} is the log-likelihood and \eqn{\nu} is the
#' (effective) number of parameters of the model, \eqn{t} is the tree level 
#' \eqn{\psi_0} is the prior probability of having a non-independence copula and 
#' \eqn{q_t} is the number of non-independence copulas in tree \eqn{t}.
#' The vBIC is a consistent model 
#' selection criterion for parametric sparse vine copula models.
#'
#' @param object a fitted `vinecop` object.
#' @param psi0 baseline prior probability of a non-independence copula.
#' @export mBICV
mBICV <- function(object, psi0 = 0.9) {
    assert_that(inherits(object, "vinecop_dist"), is.number(psi0))
    if (is.null(object$data))
        stop("data have not been stored, use keep_data = TRUE when fitting.")
    vinecop_mbicv_cpp(object$data, object, psi0)
}

#' @export
print.vinecop <- function(x, ...) {
    d <- nrow(x$matrix)
    cat(d, "-dimensional vine copula fit ('vinecop')", sep = "")
    n_trees <- length(x$pair_copulas)
    if (n_trees < d - 1)
        cat(", ", n_trees, "-truncated", sep = "")
    cat("\n")
    cat("nobs =", x$nobs, "  ")
    if (!is.null(x$data)) {
        info <- vinecop_fit_info(x)
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
summary.vinecop <- function(object, ...) {
    mdf <- summary.vinecop_dist(object)
    
    d <- dim(object)
    n_trees <- length(object$pair_copulas)
    k <- 1
    for (t in seq_len(n_trees)) {
        for (e in seq_len(d - t)) {
            mdf$loglik[k] <- object$pair_copulas[[t]][[e]]$loglik
            k <- k + 1
        }
    }
    
    mdf
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

#' extract a truncated sub-vine based on truncation level supplied by user.
#' @param object a `vinecop` or a `vine` object.
#' @param trunc_lvl truncation level for the vine copula.
#'
#' @export
truncate_model <- function(object, trunc_lvl = NA) { 
    assert_that(inherits(object, "vinecop_dist") || 
                    inherits(object, "vine_dist"))
    is_vinecop <- inherits(object, "vinecop_dist")
    if (is_vinecop) {
        pcs <- object$pair_copulas
    } else {
        pcs <- object$copula$pair_copulas
    }
    d <- length(pcs[[1]]) + 1
    if (!all(is.na(trunc_lvl)))
        assert_that(is.number(trunc_lvl), trunc_lvl <= d - 1, trunc_lvl > 0)
    
    if (!is.na(trunc_lvl)) {
        n_trees <- length(pcs)
        trunc_lvl <- min(trunc_lvl, n_trees)
        
        # truncate pair_copulas
        if (is_vinecop) {
            object$pair_copulas <- pcs[seq_len(trunc_lvl)]
        } else {
            object$copula$pair_copulas <- pcs[seq_len(trunc_lvl)]
        }
        
        # adjust npars for truncation
        pcs <- unlist(pcs[min(trunc_lvl+1,n_trees):n_trees], recursive = FALSE)
        npars <- ifelse(length(pcs) == 0, 0, 
                        sum(sapply(pcs, function(x) x[["npars"]])))
        object$npars <- object$npars - npars
        
        # adjust loglik for truncation
        if (!is.na(object$loglik)) {
            loglik <- ifelse(length(pcs) == 0, 0, 
                             sum(sapply(pcs, function(x) x[["loglik"]])))
            object$loglik <- object$loglik - loglik
        }
        
        # adjust copula object for truncation
        if (!is_vinecop) {
            object$copula$npars <- object$copula$npars - npars
            if (!is.na(object$loglik)) {
                object$copula$loglik <- object$copula$loglik - loglik
            }
        }
    } else {
        if (is_vinecop) {
            object$pair_copulas <- pcs
        } else {
            object$copula$pair_copulas <- pcs
        }
    }
    return(object)
}

#' @export
dim.vinecop_dist <- function(x) {
    ncol(x$matrix)
}