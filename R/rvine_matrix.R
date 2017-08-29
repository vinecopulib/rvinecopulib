#' R-vine matrices
#' 
#' R-vine matrices are compressed representations of the vine structure. It 
#' needs to satisfy several properties that can be checked by
#' `check_rvine_matrix()`, see *Details*.
#' 
#' The R-vine matrix notation in vinecopulib is different from the one in 
#' the VineCopula package. An example matrix is
#' ```
#' 1 1 1 1
#' 2 2 2 0
#' 3 3 0 0
#' 4 0 0 0
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
#' Denoting by `M[i, j]` the matrix entry in row `i` and column `j` (the 
#' pair-copula index for edge `e` in tree `t` of a `d` dimensional vine is
#' `(M[d - 1 - t, e], M[t, e]; M[t - 1, e], ..., M[0, e])`. Less formally,
#' 1. Start with the counter-diagonal element of column `e` (first conditioned
#'                                                           variable).
#' 2. Jump up to the element in row `t` (second conditioned variable).
#' 3. Gather all entries further up in column `e` (conditioning set).
#' 
#' A valid R-vine matrix must satisfy several conditions which are checked
#' when `RVineMatrix()` is called:
#' 1. The lower right triangle must only contain zeros.
#' 2. The upper left triangle can only contain numbers between 1 and d.
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
#' @param matrix a quadratic matrix, see *Details*.
#'
#' @return Throws an error if `matrix` is not a valid R-vine matrix, otherwise
#'     `TRUE` is returned invisibly.
#' @export
#'
#' @examples
#' mat <- matrix(c(1, 2, 3, 4, 1, 2, 3, 0, 1, 2, 0, 0, 1, 0, 0, 0), 4, 4)
#' check_rvine_matrix(mat)
#' 
#' # throws an error
#' mat[4, 4] <- 5
#' try(check_rvine_matrix(mat))
#' 
check_rvine_matrix <- function(matrix) {
    stopifnot(is.matrix(matrix))
    rvine_matrix_check_cpp(matrix)
    invisible(TRUE)
}



