#' Internal: Turn vector input into a matrix with two columns
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
    if (!is.matrix(u))
        u <- as.matrix(u)

    u
}

#' Internal: Convert arguments to `bicop_dist` object.
#' @param family the family as passed in function call.
#' @param rotation the rotation as passed in function call.
#' @param parameters the parameters as passed in function call.
#' @return A `bicop_dist` object.
#' @noRd
args2bicop <- function(family, rotation, parameters) {
    if (all(inherits(family, "bicop_dist"))) {
        return(family)
    } else {
        return(bicop_dist(family, rotation, parameters))
    }
}

#' Internal: Expand shortcuts in the familyset.
#' @noRd
expand_familyset <- function(familyset) {
    unique(sapply(familyset, expand_family))
}

expand_family <- function(family) {
    fams <- switch(
        family,
        "arch"   = family_set_archimedean,
        "ellip"  = family_set_elliptical,
        "bb"     = family_set_bb,
        "onepar" = family_set_onepar,
        "twopar" = family_set_twopar,
        "par"    = family_set_parametric,
        "nonpar" = family_set_nonparametric,
        "all"    = family_set_all,
        family  # default is no expansion
    )
}
