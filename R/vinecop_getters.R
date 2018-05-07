#' Extracts components of the vinecop object
#' 
#' Extracts either the structure matrix, or pair-copulas,
#' their parameters, Kendall's taus, or families
#' 
#' @name vinecop_getters
#' @aliases get_pair_copula get_all_pair_copulas get_parameters get_all_parameters
#' get_kendall_tau get_all_kendall_taus get_matrix
#' @param object a `vinecop` or `vine` object.
#' @details 
#' The [get_matrix] method extracts the structure matrix (see 
#' [check_rvine_matrix] for more details).
#' 
#' The other `get_xyz` methods return the entries corresponding to the 
#' pair-copula indexed by its `tree` and `edge`: \cr
#' 
#' [get_pair_copula] = the pair-copula itself (see [bicop]).\cr
#' [get_parameters] = the parameters of the pair-copula (i.e., a `numeric` 
#'  scalar, vector, or matrix).\cr
#' [get_family] = a character for the family (see [bicop] for implemented families).\cr
#' [get_kendall_tau] = a `numeric` scalar with the pair-copula Kendall's tau.\cr
#' 
#' The `get_all_xyz` methods return lists of lists, with each element 
#' corresponding to a tree in `trees`, and then elements of the sublists 
#' correspond to edges.
#' 
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
#' # get the structure matrix
#' all.equal(get_matrix(vc), mat)
#' 
#' # get pair-copulas
#' get_pair_copula(vc, 1, 1)
#' get_all_pair_copulas(vc)
#' all.equal(get_all_pair_copulas(vc), pcs, check.attributes = FALSE)
#' 
#' @return 
#' The structure matrix, or pair-copulas, their parameters, Kendall's taus, 
#' or families.
#' @rdname vinecop_getters
#' @export
get_matrix <- function(object) {
    assert_that(inherits(object, "vinecop_dist") || 
                    inherits(object, "vine_dist"))
    if (inherits(object, "vinecop_dist")) {
        return(object$matrix)
    } else {
        return(object$copula$matrix)
    }  
    
}

#' @param tree tree index.
#' @param edge edge index.
#'
#' @rdname vinecop_getters
#' @export
get_pair_copula <- function(object, tree, edge) {

    ## sanity checks
    assert_that(inherits(object, "vinecop_dist") || 
                    inherits(object, "vine_dist"))
    
    d <- dim(object)
    assert_that(is.numeric(tree), 
                is.scalar(tree),
                tree >= 1, 
                tree <= d - 1)
    assert_that(is.numeric(edge), 
                is.scalar(edge),
                edge >= 1, 
                edge <= d - tree)
    
    ## return pair-copula
    if (inherits(object, "vinecop_dist")) {
        return(object$pair_copulas[[tree]][[edge]])
    } else {
        return(object$copula$pair_copulas[[tree]][[edge]])
    }  
}

#' @rdname vinecop_getters
#' @export
get_parameters <- function(object, tree, edge) {
    get_pair_copula(object, tree, edge)$parameters
}

#' @rdname vinecop_getters
#' @export
get_tau <- function(object, tree, edge) {
    pc <- get_pair_copula(object, tree, edge)
    par_to_tau(pc$family, pc$rotation, pc$parameters)
}

#' @rdname vinecop_getters
#' @export
get_family <- function(object, tree, edge) {
    get_pair_copula(object, tree, edge)$family
}

#' @param trees the trees to extract from `object` (`trees = NA` extracts all 
#' trees).
#' @rdname vinecop_getters
#' @export
get_all_pair_copulas <- function(object, trees = NA) {
    
    assert_that(inherits(object, "vinecop_dist") || 
                    inherits(object, "vine_dist"))
    d <- dim(object)
    
    if (inherits(object, "vinecop_dist")) {
        pcs <- object$pair_copulas
    } else {
        pcs <- object$copula$pair_copulas
    }
    
    if (!any(is.na(trees)))
        assert_that(is.numeric(trees), all(trees >= 1), all(trees <= d - 1))
    
    t <- length(pcs)
    if (any(is.na(trees))) {
        trees <- seq_len(t)
    } else {
        if (any(trees > t)) {
            warning("vine copula is ", t, "-truncated; ",
                    "only returning available trees.")
            trees <- trees[trees <= t]
        }
    }
    res <- structure(pcs[trees], class = c("rvine_list", "pair_copulas"))
    attr(res, "d") <- d
    attr(res, "trees") <- trees
    return(res)
}

get_all_xyz <- function(object, trees, func, list_name) {
    pcs <- get_all_pair_copulas(object, trees)
    res <- structure(lapply(pcs, function(tree) lapply(tree, func)), 
                     class = c("rvine_list", list_name))
    attr(res, "d") <- attr(pcs, "d")
    attr(res, "trees") <- attr(pcs, "trees")
    return(res)
}

#' @rdname vinecop_getters
#' @export
get_all_parameters <- function(object, trees = NA) {
    get_all_xyz(object, trees, coef, "parameters")
}

#' @rdname vinecop_getters
#' @export
get_all_kendall_taus <- function(object, trees = NA) {
    get_all_xyz(object, trees, par_to_tau, "kendall_taus")
}

#' @rdname vinecop_getters
#' @export
get_all_families <- function(object, trees = NA) {
    get_all_xyz(object, trees, function(pc) pc$family, "families")
}

#' @export
print.rvine_list <- function(x, ...) {
    what <- switch(class(x)[2],
                   pair_copulas = "pair-copulas",
                   parameters = "copula parameters",
                   kendall_taus = "Kendall's taus",
                   families = "copula families")

    d <- attr(x, "d")
    trees <- attr(x, "trees")
    msg <- paste("Nested list of lists for the", what, "of a", 
                 d, "dimensional vine with")
    ntrees <- length(trees)
    if (ntrees == d - 1) {
        cat(paste(msg, "all trees. \n"))
    } else if (ntrees == 1) {
        cat(paste(msg, " tree ", trees, ". \n", sep = ""))
    } else {
        trees <- paste(paste(trees[1:(ntrees - 1)], collapse = ", "),
                       "and", trees[ntrees])
        cat(paste(msg, " trees ", trees, ". \n", sep = ""))
    }

}
