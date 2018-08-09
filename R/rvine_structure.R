#' R-vine structure
#' 
#' R-vine structures are compressed representations encoding the tree 
#' structure of the vine, i.e. the conditioned/conditioning 
#' variables of each edge.
#' 
#' The R-vine structure is essentially a lower-triangular matrix/triangular array, 
#' with a notation that differs from the one in the VineCopula package. 
#' An example array is
#' ```
#' 1 1 1 1
#' 2 2 2
#' 3 3
#' 4
#' ```
#' which encodes the following pair-copulas:
#'   
#' \tabular{lll}{
#' tree \tab  edge \tab pair-copulas   \cr
#' 0    \tab  0    \tab `(4, 1)`       \cr
#'      \tab  1    \tab `(3, 1)`       \cr
#'      \tab  2    \tab `(2, 1)`       \cr
#' 1    \tab  0    \tab `(4, 2; 1)`    \cr
#'      \tab  1    \tab `(3, 2; 1)`    \cr
#' 2    \tab  0    \tab `(4, 3; 2, 1)` 
#' }
#' 
#' An R-vine structure can be converted to an R-vine matrix using 
#' [as_rvine_matrix()], which encodes the same model with a square matrix 
#' filled with zeros. For instance, the matrix corresponding to the structure 
#' above is:
#' ```
#' 1 1 1 1
#' 2 2 2 0
#' 3 3 0 0
#' 4 0 0 0
#' ```
#' Similarly, an R-vine matrix can be converted to an R-vine structure using 
#' [as_rvine_structure()].
#' 
#' Denoting by `M[i, j]` the array entry in row `i` and column `j` (the 
#' pair-copula index for edge `e` in tree `t` of a `d` dimensional vine is
#' `(M[d - 1 - t, e], M[t, e]; M[t - 1, e], ..., M[0, e])`. Less formally,
#' 1. Start with the counter-diagonal element of column `e` (first conditioned
#'                                                           variable).
#' 2. Jump up to the element in row `t` (second conditioned variable).
#' 3. Gather all entries further up in column `e` (conditioning set).
#' 
#' Internally, the diagonal is stored separately from the off-diagonal elements, 
#' which are stored as a triangular array. For instance, the off-diagonal elements 
#' off the structure above are stored as
#' ```
#' 1 1 1
#' 2 2
#' 3
#' ```
#' for the structure above. The reason is that it allows for parsimonious 
#' representations of truncated models. For instance, the 2-truncated model 
#' is represented by the same diagonal and the following truncated triangular 
#' array:
#' ```
#' 1 1 1
#' 2 2
#' ```
#' 
#' A valid R-vine structure or matrix must satisfy several conditions which 
#' are checked when [rvine_structure()], [rvine_matrix()], or some coercion 
#' methods (see [as_rvine_structure()] and `as_rvine_matrix(`) are called:
#' 1. It can only contain numbers between 1 and d (and additionally zeros for 
#' R-vine matrices).
#' 3. The anti-diagonal must contain the numbers 1, ..., d.
#' 4. The anti-diagonal entry of a column must not be contained in any
#' column further to the right.
#' 5. The entries of a column must be contained in all columns to the left.
#' 6. The proximity condition must hold: For all t = 1, ..., d - 2 and
#' e = 0, ..., d - t - 1 there must exist an index j > d, such that
#' `(M[t, e], {M[0, e], ..., M[t-1, e]})` equals either
#' `(M[d-j-1, j], {M[0, j], ..., M[t-1, j]})` or
#' `(M[t-1, j], {M[d-j-1, j], M[0, j], ..., M[t-2, j]})`.
#' 
#' Condition 6 already implies conditions 2-5, but is more difficult to
#' check by hand. 
#'
#' @param order a vector of positive integers.
#' @param struct_array a list of vectors of positive integers. The vectors 
#' represent rows of the r-rvine structure and the number of elements have to 
#' be compatible with the `order` vector.
#' @param is_natural_order whether `struct_array` is assumed to be provided 
#' in natural order already.
#'
#' @return Either an `rvine_structure` or an `rvine_matrix`.
#' @export
#' @seealso as_rvine_structure
#' @examples
#' 
#' # R-vine structures can be constructed from the order vector and struct_array
#' rvine_structure(order = 1:4, struct_array = list(c(1, 1, 1), 
#'                                                  c(2, 2), 
#'                                                  3))
#' 
#' # R-vine matrices can be constructed from standard matrices
#' mat <- matrix(c(1, 2, 3, 4, 1, 2, 3, 0, 1, 2, 0, 0, 1, 0, 0, 0), 4, 4)
#' rvine_matrix(mat)
#' 
#' # coerce to R-vine structure
#' str(as_rvine_structure(mat))
#' 
#' # truncate and construct the R-vine matrix
#' mat[3, 1] <- 0
#' rvine_matrix(mat)
#' 
#' # or use directly the R-vine structure constructor
#' rvine_structure(order = 1:4, struct_array = list(c(1, 1, 1), 
#'                                                  c(2, 2)))
#' 
#' # throws an error
#' mat[3, 1] <- 5
#' try(rvine_matrix(mat))
#' 
#' @name rvine_structure
#' @aliases rvine_matrix is.rvine_structure is.rvine_matrix
rvine_structure <- function(order, struct_array, is_natural_order = TRUE) {
    
    # sanity checks and extract dimension/trunc_lvl
    assert_that(is.vector(order) && all(sapply(order, is.count)),
                msg = "Order should be a vector of positive integers.")
    assert_that(is.vector(struct_array) || is.list(struct_array))
    if (is.list(struct_array))
        struct_array <- unlist(struct_array)
    assert_that(all(sapply(struct_array, is.count)), 
                msg = "All elements of struct_array should be positive integers.")
    d <- length(order)
    dd <- cumsum((d-1):1)
    assert_that(length(struct_array) %in% dd, 
                msg = "The number of elements in struct_array is incompatible 
                with order.")
    trunc_lvl <- which(dd == length(struct_array))
    
    # create column-wise structure array
    dd <- c(1, 1 + dd[-(d-1)])
    struct_array <- lapply(1:(d - 1), function(i) 
        struct_array[dd[1:min((d - i), trunc_lvl)] + (i - 1)])
    
    # create and check output
    output <- structure(list(order = order, 
                             struct_array = struct_array,
                             d = d,
                             trunc_lvl = trunc_lvl),
                        class = c("rvine_structure", "list"))
    validate_rvine_structure(output, is_natural_order)
    
    # return output
    output
}

#' @param matrix an R-vine matrix, see *Details*.
#' @rdname rvine_structure
#' @export
rvine_matrix <- function(matrix) {
    validate_rvine_matrix(matrix)
    rvine_matrix_nocheck(matrix)
}

rvine_matrix_nocheck <- function(matrix) {
    d <- ncol(matrix)
    class(matrix) <- c("rvine_matrix", class(matrix))
    attr(matrix, "d") <- d
    attr(matrix, "trunc_lvl") <- ifelse(any(matrix[, 1] == 0), 
                                        min(which(matrix[, 1] == 0)),
                                        d - 1)
    matrix
}

validate_rvine_structure <- function(structure, is_natural_order = TRUE) {
    assert_that(is.rvine_structure(structure))
    assert_that(is.flag(is_natural_order))
    rvine_structure_check_cpp(structure, is_natural_order)
}

validate_rvine_matrix <- function(matrix) {
    assert_that(is.rvine_matrix(matrix) || 
                    (is.matrix(matrix) && is.numeric(matrix)))
    rvine_matrix_check_cpp(matrix)
}

#' @export
dim.rvine_matrix <- function(x) {
    output <- c(attr(x, "d"), attr(x, "trunc_lvl"))
    names(output) <- c("dim", "trunc_lvl")
    output
}

#' @export
#' @importFrom utils write.table
print.rvine_matrix <- function(x, ..., zero.print = " ", n = 10, 
                               structure = FALSE) {
    if (structure) {
        cat(dim(x)[1], "-dimensional R-vine structure ('rvine_structure')", sep = "")
    } else {
        cat(dim(x)[1], "-dimensional R-vine matrix ('rvine_matrix')", sep = "")
    }
    print_truncation_info(x)
    write.table(format(x, zero.print = zero.print), quote = FALSE, 
                row.names = FALSE, col.names = FALSE)
    invisible(x)
}

#' @export
print.rvine_structure <- function(x, ...) {
    print(as_rvine_matrix(x), structure = TRUE)
}

#' @export
dim.rvine_structure <- function(x) {
    output <- c(x$d, x$trunc_lvl)
    names(output) <- c("dim", "trunc_lvl")
    output
}

#' @export
is.rvine_structure <- function(structure) {
    inherits(structure, "rvine_structure")
}

#' @export
is.rvine_matrix <- function(matrix) {
    inherits(matrix, "rvine_matrix")
}
