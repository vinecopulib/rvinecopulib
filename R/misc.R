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
#' @param object an `rvine_structure`, `vinecop` or `vine` object.
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
#' structure <- to_rvine_structure(mat)
#' 
#' # truncate the model
#' truncate_model(structure, 1)
#' 
#' # set up vine copula model
#' vc <- vinecop_dist(pcs, mat)
#' 
#' # truncate the model
#' truncate_model(vc, 1)
truncate_model <- function(object, trunc_lvl = NA) { 
    
    assert_that(inherits(object, "rvine_structure") || 
                    inherits(object, "vinecop_dist") || 
                    inherits(object, "vine_dist"))
    
    assert_that((is.number(trunc_lvl) && 
                     trunc_lvl <= dim(object)[2] && 
                     trunc_lvl >= 1) ||
                    (is.scalar(trunc_lvl) && is.na(trunc_lvl)),
                msg = "trunc_lvl should be NA or a number between 1 and 
                the object's truncation level.")
    
    if (!is.na(trunc_lvl) || trunc_lvl < dim(object)[2]) {
        
        if (inherits(object, "rvine_structure")) {
            object$struct_array <- lapply(object$struct_array, 
                                          truncate_column,
                                          trunc_lvl = trunc_lvl)
            object$trunc_lvl <- trunc_lvl
            return(object)
        }
        
        pcs <- get_all_pair_copulas(object)
        n_trees <- dim(object)[2]
        
        # truncate pair_copulas and structure
        if (inherits(object, "vinecop_dist")) {
            object$structure <- truncate_model(object$structure, trunc_lvl)
            object$pair_copulas <- pcs[seq_len(trunc_lvl)]
        } else {
            object$copula$structure <- truncate_model(object$copula$structure, trunc_lvl)
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
        if (!inherits(object, "vinecop_dist")) {
            object$copula$npars <- object$copula$npars - npars
            if (!is.na(object$loglik)) {
                object$copula$loglik <- object$copula$loglik - loglik
            }
        }
    }
    
    return(object)
}

#' internal function
#' @noRd
truncate_column <- function(column, trunc_lvl) {
    column[1:min(length(column), trunc_lvl)]
}