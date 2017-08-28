#' @param pair_copulas A nested list of 'bicop_dist' objects, where 
#'    \code{pair_copulas[[t]][[e]]} corresponds to the pair-copula at edge `e` in
#'    tree `t`.
#' @param matrix A quadratic structure matrix of dimension 
#'   `length(pair_copulas) + 1` (see *Examples* and [`check_rvine_matrix()`]). 
#' @rdname vinecop
#' @export
vinecop_dist <- function(pair_copulas, matrix) {
    # create object
    vinecop <- structure(
        list(pair_copulas = pair_copulas, matrix = matrix),
        class = "vinecop_dist"
    )
    
    # sanity checks
    stopifnot(is.list(pair_copulas))
    pc_lst <- unlist(pair_copulas, recursive = FALSE)
    if (!all(sapply(pc_lst, function(x) inherits(x, "bicop_dist")))) {
        stop("some objects in pair_copulas aren't of class 'bicop_dist'")
    }
    vinecop_check_cpp(vinecop)
    check_rvine_matrix(matrix)
    vinecop$npars <- sum(sapply(pc_lst, function(x) x[["npars"]]))
    
    vinecop
}

#' Vine copula distributions
#' 
#' Density, distribution function and random generation 
#' for the vine copula distribution.
#' 
#' @aliases vinecop_methods dvinecop pvinecop rvinecop dvinecop_dist pvinecop_dist rvinecop_dist
#' @param u evaluation points, either a length d vector or a d-column matrix,
#'   where d is the number of variables in the vine.
#' @param vinecop an object of class `"vinecop_dist"`.
#' @return 
#' `dvinecop` gives the density, `pvinecop` gives the distribution function, 
#' and `rvinecop` generates random deviates.
#' 
#' The length of the result is determined by `n` for `rvinecop`, and 
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
#' @export
rvinecop <- function(n, vinecop) {
    stopifnot(inherits(vinecop, "vinecop_dist"))
    vinecop_sim_cpp(n, vinecop)
}

#' @export
print.vinecop_dist <- function(x, ...) {
    d <- nrow(x$matrix)
    cat(d, "-dimensional vine copula model", sep = "")
}

#' @export
summary.vinecop_dist <- function(object, ...) {
    mat <- object$matrix
    d <- nrow(mat)
    mdf <- as.data.frame(matrix(NA, choose(d, 2), 7))
    names(mdf) <- c("tree", "edge", 
                    "conditioned", "conditioning", 
                    "family", "rotation", "parameters")
    k <- 1
    for (t in seq.int(d - 1)) {
        for (e in seq.int(d - t)) {
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
    print(mdf, digits = 2, row.names = FALSE)
    
    invisible(mdf)
}
