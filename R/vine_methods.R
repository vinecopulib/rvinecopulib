#' Vine based distributions
#' 
#' Density, distribution function and random generation 
#' for the vine based distribution.
#' 
#' @name vine_distributions
#' @aliases dvine pvine rvine dvine_dist pvine_dist rvine_dist
#' @param x evaluation points, either a length d vector or a d-column matrix,
#'   where d is the number of variables in the vine.
#' @param vine an object of class `"vine_dist"`.
#' @details 
#' See [vine] for the estimation and construction of vine models. 
#' Here, the density, distribution function and random generation 
#' for the vine distributions are standard.
#' 
#' The functions are based on [dvinecop()], [pvinecop()] and [rvinecop()] for 
#' [vinecop] objects, and either [kde1d::dkde1d()], [kde1d::pkde1d()] and 
#' [kde1d::qkde1d()] for estimated vines (i.e., output of [vine()]), or the 
#' standard *d/p/q-xxx* from [stats::Distributions] for customly created vines 
#' (i.e., output of [vine_dist()])
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
#' vc <- vine_dist(list(name = "norm"), pcs, mat)
#' 
#' # simulate from the model
#' x <- rvine(200, vc)
#' pairs(x)
#' 
#' # evaluate the density and cdf
#' dvine(x[1, ], vc)
#' pvine(x[1, ], vc)
#' @rdname vine_methods
#' @importFrom cctools expand_as_numeric
#' @export
dvine <- function(x, vine) {
    stopifnot(inherits(vine, "vine_dist"))
    if (NCOL(x) == 1)
        x <- t(x)
    nms <- colnames(x)
    # must be numeric, factors are expanded
    x <- expand_as_numeric(x)
    # variables must be in same order
    if (!is.null(nms))
        x <- x[, colnames(vine$data), drop = FALSE]
    
    # prepare marginals if only one is specified
    d <- ncol(x)
    if (!inherits(vine, "vine") & depth(vine$marg) == 1) 
        vine$marg <- replicate(d, vine$marg, simplify = FALSE)
    
    ## evaluate marginal densities
    margvals <- u <- x
    for (k in 1:d) {
        x_k <- x[, k]
        if (inherits(vine, "vine")) {
            if (k %in% attr(vine$data, "i_disc")) {
                # use normalization if discrete
                attr(x_k, "i_disc") <- 1
                vine$marg[[k]]$levels <- attr(vine$data, "levels")[[k]]
            }
            margvals[, k] <- dkde1d(x_k, vine$marg[[k]])
        } else {
            dfun <- get(paste0("d", vine$marg[[k]]$name))
            par <- vine$marg[[k]][names(vine$marg[[k]]) != "name"]
            par$x <- x_k
            margvals[, k] <- do.call(dfun, par)
        }
    }
    
    if (!is.null(vine$cop)) {
        # PIT to copula data
        u <- get_u(x, vine)
        # evaluate vine density
        vinevals <- vinecop_pdf_cpp(u, vine$cop)
    } else {
        vinevals <- rep(1, nrow(x))
    }
    
    ## final density estimate is product of marginals and copula density
    apply(cbind(margvals, vinevals), 1, prod)
}

#' @rdname vine_methods
#' @param n_mc number of samples used for quasi Monte Carlo integration.
#' @export
pvine <- function(x, vine, n_mc = 10^4) {
    
    stopifnot(inherits(vine, "vine_dist"))
    
    if (NCOL(x) == 1)
        x <- t(x)
    nms <- colnames(x)
    # must be numeric, factors are expanded
    x <- expand_as_numeric(x)
    # variables must be in same order
    if (!is.null(nms))
        x <- x[, colnames(vine$data), drop = FALSE]
    
    # prepare marginals if only one is specified
    if (!inherits(vine, "vine") & depth(vine$marg) == 1) 
        vine$marg <- replicate(ncol(x), vine$marg, simplify = FALSE)
    
    # PIT to copula data
    u <- get_u(x, vine)
    
    # Evaluate copula if needed
    if (!is.null(vine$cop)) {
        vals <- vinecop_cdf_cpp(u, vine$cop, n_mc)
    } else {
        vals <- u
    }
    
    return(vals)
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
    
    # prepare uniform data
    d <- ncol(vine$cop$matrix)
    U <- prep_uniform_data(n, d, U)

    # simulate from copula
    U <- vinecop_inverse_rosenblatt_cpp(U, vine$cop)
    
    # prepare marginals if only one is specified
    if (!inherits(vine, "vine") & depth(vine$marg) == 1) 
        vine$marg <- replicate(d, vine$marg, simplify = FALSE)
    
    # use quantile transformation for marginals
    if (inherits(vine, "vine")) {
        U <- sapply(seq_len(d), function(i) qkde1d(U[, i], vine$marg[[i]]))
    } else {
        U <- sapply(seq_len(d), function(i) {
            qfun <- get(paste0("q", vine$marg[[i]]$name))
            par <- vine$marg[[i]][names(vine$marg[[i]]) != "name"]
            par$p <- U[, i]
            do.call(qfun, par)
        })
    }

    U
}

#' @export
print.vine_dist <- function(x, ...) {
    print(x$cop)
}

#' @export
summary.vine_dist <- function(object, ...) {
    summary(object$cop)
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
#' x <- sapply(1:5, function(i) rnorm(50))
#' fit <- vine(x, cop_controls = list(family_set = "par"))
#' all.equal(predict(fit, x), fitted(fit))
predict.vine <- function(object, newdata, what = "pdf", n_mc = 10^4, ...) {
    stopifnot(what %in% c("pdf", "cdf"))
    switch(
        what,
        "pdf" = dvine(newdata, object),
        "cdf" = pvine(newdata, object, n_mc)
    )
}

#' @rdname predict_vine
#' @export
fitted.vine <- function(object, what = "pdf", n_mc = 10^4, ...) {
    if (all(is.na(object$data)))
        stop("data have not been stored, use keep_data = TRUE when fitting.")
    stopifnot(what %in% c("pdf", "cdf"))
    switch(
        what,
        "pdf" = dvine(object$data, object),
        "cdf" = pvine(object$data, object, n_mc)
    )
}

#' @export
logLik.vine <- function(object, ...) {
    if (all(is.na(object$data)))
        stop("data have not been stored, use keep_data = TRUE when fitting.")
    
    ll_marg <- lapply(object$marg, logLik)
    npars_marg <- sapply(ll_marg, function(x) attr(x, "df"))
    u <- get_u(object$data, object)
    ll_cop <- vinecop_loglik_cpp(u, object$cop)
    pc_lst <- unlist(object$cop$pair_copulas, recursive = FALSE)
    npars_cop <- ifelse(length(pc_lst) == 0, 0, 
                    sum(sapply(pc_lst, function(x) x[["npars"]])))
    structure(sum(unlist(ll_marg)) + ll_cop, "df" = sum(npars_marg) + npars_cop)
}

#' @export
print.vine <- function(x, ...) {
    print(collate_u(x)$cop)
}


#' @export
summary.vine <- function(object, ...) {
    summary(collate_u(object)$cop)
}

# PIT to copula level
get_u <- function(x, vine) {
    d <- ncol(x)
    u <- matrix(NA, nrow(x), ncol(x))
    
    for (k in 1:d) {
        x_k <- x[, k]
        if (inherits(vine, "vine")) {
            if (k %in% attr(vine$data, "i_disc")) {
                # use continuous variant for PIT
                attr(x_k, "i_disc") <- integer(0)
                vine$marg[[k]]$levels <- NULL
            }
            u[, k] <- pkde1d(x_k, vine$marg[[k]])
        } else {
            pfun <- get(paste0("p", vine$marg[[k]]$name))
            par <- vine$marg[[k]][names(vine$marg[[k]]) != "name"]
            par$q <- x_k
            u[, k] <- do.call(pfun, par)
        }
    }
    return(u)
}

collate_u <- function(x) {
    if (!all(is.na(x$data))) {
        u <- get_u(x$data, x)
        x$cop$data <- u
    }
    x
}