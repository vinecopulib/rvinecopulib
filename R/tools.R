#' Internal: Turn vector input into a matrix with two columns
#'
#' @param u input data
#'
#' @return either a matrix with two columns, or an error if u is neither a
#' matrix, data.frame, or a length two vector
#'
#' @noRd
if_vec_to_matrix <- function(u) {
    if (NCOL(u) == 1)
        u <- matrix(u, 1, length(u))
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
        if (missing(rotation))
            rotation <- 0
        if (missing(parameters))
            parameters <- numeric(0)
        return(bicop_dist(family, rotation, parameters))
    }
}

#' Internal: Expand shortcuts in the familyset.
#' @noRd
expand_family_set <- function(family_set) {
    unique(unlist(lapply(family_set, expand_family)))
}

expand_family <- function(family) {
    switch(
        family,
        "archimedean"   = family_set_archimedean,
        "ellipiltical"  = family_set_elliptical,
        "bbs"           = family_set_bb,
        "oneparametric" = family_set_onepar,
        "twoparametric" = family_set_twopar,
        "parametric"    = family_set_parametric,
        "nonparametric" = family_set_nonparametric,
        "all"           = family_set_all,
        family  # default is no expansion
    )
}

check_family_set <- function(family_set) {
    i_wrong <- which(!(family_set %in% family_set_all_defs))
    if (length(i_wrong) > 0) {
        stop(
            "unknown families in family_set: ",
            paste0('"', family_set[i_wrong], '"', collapse = ", ")
        )
    }
}

prep_uniform_data <- function(n, d, U) {
    if (is.null(U)) {
        U <- matrix(runif(n * d), n, d)
    } else {
        stopifnot(is.matrix(U))
        stopifnot(nrow(U) == n)
        stopifnot(ncol(U) == d)
    }
}
