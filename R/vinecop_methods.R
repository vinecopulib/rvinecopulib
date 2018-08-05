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
#' @param cores number of cores to use; if larger than one, computations are
#'   done in parallel on `cores` batches .
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
#' @rdname vinecop_methods
#' @export
dvinecop <- function(u, vinecop, cores = 1) {
    assert_that(inherits(vinecop, "vinecop_dist"))
    vinecop_pdf_cpp(if_vec_to_matrix(u), vinecop, cores)
}

#' @rdname vinecop_methods
#' @param n_mc number of samples used for quasi Monte Carlo integration.
#' @importFrom assertthat is.count
#' @export
pvinecop <- function(u, vinecop, n_mc = 10^4, cores = 1) {
    assert_that(inherits(vinecop, "vinecop_dist"), 
                is.number(n_mc), is.count(cores))
    vinecop_cdf_cpp(if_vec_to_matrix(u), vinecop, n_mc, cores, get_seeds())
}

#' @rdname vinecop_methods
#' @param n number of observations.
#' @param U optionally, an \eqn{n \times d} matrix of values in \eqn{(0,1)}.
#'    The result is then the inverse Rosenblatt transform of `U`; if `U` is a
#'    matrix of independent \eqn{U(0, 1)} variables, this simulates data 
#'    from `vinecop`.
#' @param qrng if `TRUE`, generates quasi-random numbers using the multivariate 
#' Generalized Halton sequence up to dimension 300 and the Generalized Sobol 
#' sequence in higher dimensions (default `qrng = FALSE`).
#' @export
rvinecop <- function(n, vinecop, U = NULL, qrng = FALSE, cores = 1) {
    assert_that(inherits(vinecop, "vinecop_dist"))
    check_u_and_qrng(U, qrng, n, ncol(vinecop$matrix))
    
    U <- vinecop_sim_cpp(vinecop, n, qrng, cores, get_seeds())
    if (!is.null(vinecop$names))
        colnames(U) <- vinecop$names
    
    U
}

#' @export
print.vinecop_dist <- function(x, ...) {
    cat(dim(x), "-dimensional vine copula model ('vinecop_dist')", sep = "")
    print_truncation_info(x)
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
    class(mdf) <- c("summary_df", class(mdf))
    mdf
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
#' @param cores number of cores to use; if larger than one, computations are
#'   done in parallel on `cores` batches .
#' @param ... unused.
#' 
#' @details `fitted()` can only be called if the model was fit with the
#'    `keep_data = TRUE` option.
#' 
#' @return 
#' `fitted()` and `predict()` have return values similar to [dvinecop()] 
#' and [pvinecop()].
#' @export
#' @rdname predict_vinecop
#' @examples
#' u <- sapply(1:5, function(i) runif(50))
#' fit <- vinecop(u, "par", keep_data = TRUE)
#' all.equal(predict(fit, u), fitted(fit))
predict.vinecop <- function(object, newdata, what = "pdf", n_mc = 10^4, 
                            cores = 1, ...) {
    assert_that(
        in_set(what, c("pdf", "cdf")), 
        is.number(n_mc),
        is.number(cores), cores > 0
    )
    newdata <- if_vec_to_matrix(newdata)
    switch(
        what,
        "pdf" = vinecop_pdf_cpp(newdata, object, cores),
        "cdf" = vinecop_cdf_cpp(newdata, object, n_mc, cores, get_seeds())
    )
}

#' @rdname predict_vinecop
#' @export
fitted.vinecop <- function(object, what = "pdf", n_mc = 10^4, cores = 1, ...) {
    if (is.null(object$data))
        stop("data have not been stored, use keep_data = TRUE when fitting.")
    assert_that(
        in_set(what, c("pdf", "cdf")), 
        is.number(n_mc), 
        is.number(cores), cores > 0
    )
    switch(
        what,
        "pdf" = vinecop_pdf_cpp(object$data, object, cores),
        "cdf" = vinecop_cdf_cpp(object$data, object, n_mc, cores, get_seeds())
    )
}

#' @export
logLik.vinecop <- function(object, ...) {
    structure(object$loglik, "df" = object$npars)
}

#' Modified vine copula Bayesian information criterion (mBICv)
#' 
#' Calculates the modified vine copula Bayesian information criterion.
#' 
#' The modified vine copula Bayesian information criterion (mBICv) is defined as
#' 
#' \deqn{\mathrm{BIC} = -2\, \mathrm{loglik} +  \nu \log(n) - 2 * 
#' \sum_{t=1}^{d - 1} \{q_t \log(\psi_0^t) - (d - t - q_t) \log(1 - \psi_0^t)\}
#' }
#' 
#' where \eqn{\mathrm{loglik}} is the log-likelihood and \eqn{\nu} is the
#' (effective) number of parameters of the model, \eqn{t} is the tree level 
#' \eqn{\psi_0} is the prior probability of having a non-independence copula and 
#' \eqn{q_t} is the number of non-independence copulas in tree \eqn{t}.
#' The mBICv is a consistent model selection criterion for parametric sparse 
#' vine copula models.
#'
#' @param object a fitted `vinecop` object.
#' @param psi0 baseline prior probability of a non-independence copula.
#' @param newdata optional; a new data set.
#' @export mBICV
#' @examples
#' u <- sapply(1:5, function(i) runif(50))
#' fit <- vinecop(u, "par", keep_data = TRUE)
#' mBICV(fit, 0.9) # with a 0.9 prior probability of a non-independence copula
#' mBICV(fit, 0.1) # with a 0.1 prior probability of a non-independence copula
mBICV <- function(object, psi0 = 0.9, newdata = NULL) {
    assert_that(inherits(object, "vinecop_dist"), is.number(psi0))
    ll <- ifelse(is.null(newdata), 
                 object$loglik, 
                 sum(log(dvinecop(newdata, object))))
    - 2 * ll + compute_mBICV_penalty(object, psi0) 
}

compute_mBICV_penalty <- function(object, psi0) {
    d <- dim(object)
    smr <- summary(object)
    q_m <- tapply(smr$family, smr$tree, function(x) sum(x == "indep"))
    q_m <- c(q_m, rep(0, d - 1 - length(q_m)))
    m_seq <- seq_len(d - 1)
    pen <- object$npars * log(object$nobs)
    pen - 2 * sum(q_m * log(psi0^m_seq) + (d - 1 - q_m) * log(1 - psi0^m_seq))
}

#' @export
print.vinecop <- function(x, ...) {
    cat(dim(x), "-dimensional vine copula fit ('vinecop')", sep = "")
    print_truncation_info(x)
    print_fit_info(x)
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

#' Truncate a vine copula model
#' 
#' Extracts a truncated sub-vine based on a truncation level supplied by user.
#' 
#' While a vine model for a `d` dimensional random vector contains at most `d-1` 
#' nested trees, this function extracts a sub-model based on a given truncation 
#' level. 
#' 
#' For instance, `truncate_model(object, 1)` results in a 1-truncated 
#' vine (i.e., a vine with a single tree). Similarly `truncate_model(object, 2)` 
#' results in a 2-truncated vine (i.e., a vine with two trees). Note that 
#' `truncate_model(truncate_model(object, 1), 2)` returns a 1-truncated vine.
#' 
#' @param object a `vinecop` or a `vine` object.
#' @param trunc_lvl truncation level for the vine copula (`NA` corresponds to 
#' no truncation).
#'
#' @export
#' @examples
#' # specify pair-copulas
#' bicop <- bicop_dist("bb1", 90, c(3, 2))
#' pcs <- list(
#'     list(bicop, bicop),  # pair-copulas in first tree 
#'     list(bicop)          # pair-copulas in second tree 
#' )
#' 
#' # specify R-vine matrix
#' mat <- matrix(c(1, 2, 3, 1, 2, 0, 1, 0, 0), 3, 3) 
#' 
#' # set up vine copula model
#' vc <- vinecop_dist(pcs, mat)
#' 
#' # truncate the model
#' truncate_model(vc, 1)
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