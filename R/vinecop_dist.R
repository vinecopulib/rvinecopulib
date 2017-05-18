#' Vine copula models
#'
#' @param pair_copulas A nested list of 'bicop_dist' objects, where 
#'    `pair_copulas[[t]][[e]]` corresponds to the pair-copula at edge `e` in
#'    tree `t`.
#' @param matrix A quadratic structure matrix of dimension 
#'   `length(pair_copulas) + 1`. 
#'
#' @examples
#' 
#' @export
vinecop_dist <- function(pair_copulas = NULL, matrix = NULL) {
    # create object
    vinecop <- structure(
        list(pair_copulas = pair_copulas, matrix = matrix),
        class = "vinecop_dist"
    )
    
    # sanity checks
    if (!is.null(pair_copulas)) {
        pc_lst <- unlist(pair_copulas, recursive = FALSE)
        if (!all(sapply(pc_lst, function(x) inherits(x, "bicop_dist")))) {
            stop("some objects in pair_copulas aren't of class 'bicop_dist'")
        }
    }
    vinecop_check_cpp(vinecop)
    
    vinecop
}

#' @rdname vinecop_dist
#' @param u evaluation points, either a length d vector or a d-column matrix,
#'   where d is the number of variables in the vine.
#' @param vinecop an object of class `"vinecop"`.
#' @export
dvinecop <- function(u, vinecop) {
    vinecop_pdf_cpp(if_vec_to_matrix(u), vinecop)
}

#' @rdname vinecop_dist
#' @param n number of observations.
#' @export
rvinecop <- function(n, vinecop) {
    vinecop_sim_cpp(n, vinecop)
}
