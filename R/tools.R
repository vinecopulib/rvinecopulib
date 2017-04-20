#' Internal: Turns vector input into a matrix with two columns
#'
#' @param u input data
#'
#' @return either a matrix with two columns, or an error if u is neither a
#' matrix, data.frame, or a length two vector
#'
#' @noRd
if_vec_to_matrix <- function(u) {
    if (NCOL(u) == 1) {
        if (length(u) != 2)
            stop("u must be a mx2 matrix, data.frame or a length two vector.")
        u <- matrix(u, 1, 2)
    }
    if (!NCOL(u) == 2)
        stop("u must be a mx2 matrix, data.frame or a length two vector.")
    
    as.matrix(u)
}

args2bicop <- function(args) {
    if (all(inherits(args$family, "bicop_dist"))) {
        return(args$family)
    } else {
        return(bicop_dist(args$family, args$rotation, args$parameters))
    }
}

expand_familyset <- function(familyset) {
    unique(sapply(familyset, expand_family))
}

expand_family <- function(family) {
    switch(
        family,
        "onepar" = family_set_onepar,
        "bb"     = family_set_bb,
        "twopar" = family_set_twopar,
        "par"    = family_set_parametric,
        "nonpar" = family_set_nonpar,
        "all" = family_set_all,
        family  # default is no expansion
    )
}