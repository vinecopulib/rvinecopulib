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
#' @param object a model object.
#' @param ... further arguments passed to specific methods.
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
#' # set up vine structure
#' structure <- as_rvine_structure(mat)
#' 
#' # truncate the model
#' truncate_model(structure, 1)
#' 
#' # set up vine copula model
#' vc <- vinecop_dist(pcs, mat)
#' 
#' # truncate the model
#' truncate_model(vc, 1)
truncate_model <- function(object, ...) {
    UseMethod("truncate_model")
}

#' @export
#' @param trunc_lvl tree level after which the vine copula should be truncated.
#' @rdname truncate_model
truncate_model.rvine_structure <- function(object, trunc_lvl, ...) {
    check_trunc_lvl(object, trunc_lvl)
    if (trunc_lvl < dim(object)["trunc_lvl"]) {
        object$struct_array <- lapply(object$struct_array, 
                                      truncate_column,
                                      trunc_lvl = trunc_lvl)
        object$trunc_lvl <- trunc_lvl
    }
    object
}

#' @export
#' @rdname truncate_model
truncate_model.rvine_matrix <- function(object, trunc_lvl, ...) {
    check_trunc_lvl(object, trunc_lvl)
    d <- dim(object)["dim"]
    trees_to_truncate <- setdiff(seq_len(d - 1), seq_len(trunc_lvl))
    for (tree in trees_to_truncate) {
        object[tree, seq_len(d - tree)] <- 0
    }
    attr(object, "trunc_lvl") <- trunc_lvl
    object
}

#' @export
#' @rdname truncate_model
truncate_model.vinecop_dist <- function(object, trunc_lvl, ...) {
    check_trunc_lvl(object, trunc_lvl)
    if (trunc_lvl < dim(object)["trunc_lvl"]) {
        object <- adjust_fit_stats(object, trunc_lvl)
        object$structure <- truncate_model(object$structure, trunc_lvl)
        object$pair_copulas <- object$pair_copulas[seq_len(trunc_lvl)]
    }
    object
}

#' @export
#' @rdname truncate_model
truncate_model.vine_dist <- function(object, trunc_lvl, ...) {
    check_trunc_lvl(object, trunc_lvl)
    if (trunc_lvl < dim(object)["trunc_lvl"]) {
        object <- adjust_fit_stats(object, trunc_lvl)
        object$copula <- truncate_model(object$copula, trunc_lvl)
    }
    object
}

#' internal function
#' @noRd
truncate_column <- function(column, trunc_lvl) {
    column[1:min(length(column), trunc_lvl)]
}

#' internal function
#' @noRd
check_trunc_lvl <- function(object, trunc_lvl) {
    msg <- paste0(
        "trunc_lvl should be a number between 1 and the number of trees (",
        dim(object)["dim"] - 1, ")."
    )
    assert_that(
        is.count(trunc_lvl), 
        trunc_lvl >= 1, 
        trunc_lvl < dim(object)["dim"],
        msg = msg
    )
    if (trunc_lvl > dim(object)["trunc_lvl"])
        warning("truncation has no effect; vine is already ", 
                dim(object)["trunc_lvl"], "-truncated.", call. = FALSE)
    invisible(TRUE)
}

#' internal function
#' @noRd
get_truncated_pcs <- function(object, trunc_lvl) {
    check_trunc_lvl(object, trunc_lvl)
    if (!is.null(object$copula)) {
        pcs <- object$copula$pair_copulas[-seq_len(trunc_lvl)]
    } else {
        pcs <- object$pair_copulas[-seq_len(trunc_lvl)]
    }
    unlist(pcs, recursive = FALSE)
}

#' internal function
#' @noRd
adjust_fit_stats <- function(object, trunc_lvl) {
    trunc_pcs <- get_truncated_pcs(object, trunc_lvl)
    if (length(trunc_pcs) == 0)  # model is unchanged
        return(object)
    
    # adjust npars for truncation
    trunc_npars <- sum(sapply(trunc_pcs, function(x) x[["npars"]]))
    object$npars <- object$npars - trunc_npars
    
    # adjust loglik for truncation
    if (!is.na(object$loglik)) {
        trunc_loglik <- sum(sapply(trunc_pcs, function(x) x[["loglik"]]))
        object$loglik <- object$loglik - trunc_loglik
    }
    
    object
}