#' R-vine structure
#'
#' R-vine structures are compressed representations encoding the tree structure
#' of the vine, i.e. the conditioned/conditioning variables of each edge. The
#' functions `[cvine_structure()]` or `[dvine_structure()]` give a simpler way
#' to construct C-vines (every tree is a star) and D-vines (every tree is a
#' path), respectively (see *Examples*).
#'
#' The R-vine structure is essentially a lower-triangular matrix/triangular
#' array, with a notation that differs from the one in the VineCopula package.
#' An example array is
#' ```
#' 4 4 4 4
#' 3 3 3
#' 2 2
#' 1
#' ```
#' which encodes the following pair-copulas:
#'
#' \tabular{lll}{
#' tree \tab  edge \tab pair-copulas   \cr
#' 0    \tab  0    \tab `(1, 4)`       \cr
#'      \tab  1    \tab `(2, 4)`       \cr
#'      \tab  2    \tab `(3, 4)`       \cr
#' 1    \tab  0    \tab `(1, 3; 4)`    \cr
#'      \tab  1    \tab `(2, 3; 4)`    \cr
#' 2    \tab  0    \tab `(1, 2; 3, 4)`
#' }
#'
#' An R-vine structure can be converted to an R-vine matrix using
#' [as_rvine_matrix()], which encodes the same model with a square matrix filled
#' with zeros. For instance, the matrix corresponding to the structure above is:
#' ```
#' 4 4 4 4
#' 3 3 3 0
#' 2 2 0 0
#' 1 0 0 0
#' ```
#' Similarly, an R-vine matrix can be converted to an R-vine structure using
#' [as_rvine_structure()].
#'
#' Denoting by `M[i, j]` the array entry in row `i` and column `j` (the
#' pair-copula index for edge `e` in tree `t` of a `d` dimensional vine is
#' `(M[d + 1 - e, e], M[t, e]; M[t - 1, e], ..., M[1, e])`. Less formally,
#'
#' 1. Start with the counter-diagonal element of column `e` (first conditioned
#' variable).
#'
#' 2. Jump up to the element in row `t` (second conditioned variable).
#'
#' 3. Gather all entries further up in column `e` (conditioning set).
#'
#' Internally, the diagonal is stored separately from the off-diagonal elements,
#' which are stored as a triangular array. For instance, the off-diagonal
#' elements off the structure above are stored as
#' ```
#' 4 4 4
#' 3 3
#' 2
#' ```
#' for the structure above. The reason is that it allows for parsimonious
#' representations of truncated models. For instance, the 2-truncated model is
#' represented by the same diagonal and the following truncated triangular
#' array:
#' ```
#' 4 4 4
#' 3 3
#' ```
#'
#' A valid R-vine structure or matrix must satisfy several conditions which are
#' checked when [rvine_structure()], [rvine_matrix()], or some coercion methods
#' (see [as_rvine_structure()] and `as_rvine_matrix(`) are called:
#'
#' 1. It can only contain numbers between 1 and d (and additionally zeros for
#' R-vine matrices).
#'
#' 2. The anti-diagonal must contain the numbers 1, ..., d.
#'
#' 3. The anti-diagonal entry of a column must not be contained in any column
#' further to the right.
#'
#' 4. The entries of a column must be contained in all columns to the left.
#'
#' 5. The proximity condition must hold: For all t = 1, ..., d - 2 and e = 1,
#' ..., d - t there must exist an index j > d, such that
#' `(M[t, e], {M[1, e], ..., M[t - 1, e]})` equals either
#' `(M[d + 1 - j, j], {M[1, j], ..., M[t - 1, j]})` or
#' `(M[t - 1, j], {M[d + 1 - j, j], M[1, j], ..., M[t - 2, j]})`.
#'
#'
#' Condition 5 already implies conditions 2-4, but is more difficult to check by
#' hand.
#'
#' @param order a vector of positive integers.
#' @param struct_array a list of vectors of positive integers. The vectors
#'   represent rows of the r-rvine structure and the number of elements have to
#'   be compatible with the `order` vector. If empty, the model is 0-truncated.
#' @param is_natural_order whether `struct_array` is assumed to be provided in
#'   natural order already (a structure is in natural order if the anti-
#'   diagonal is 1, .., d from bottom left to top right).
#'
#' @return Either an `rvine_structure` or an `rvine_matrix`.
#' @export
#' @seealso [as_rvine_structure()], [as_rvine_matrix()],
#'   [plot.rvine_structure()], [plot.rvine_matrix()],
#'   [rvine_structure_sim()], [rvine_matrix_sim()]
#' @examples
#'
#' # R-vine structures can be constructed from the order vector and struct_array
#' rvine_structure(order = 1:4, struct_array = list(
#'   c(4, 4, 4),
#'   c(3, 3),
#'   2
#' ))
#'
#' # R-vine matrices can be constructed from standard matrices
#' mat <- matrix(c(4, 3, 2, 1, 4, 3, 2, 0, 4, 3, 0, 0, 4, 0, 0, 0), 4, 4)
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
#' rvine_structure(order = 1:4, struct_array = list(
#'   c(4, 4, 4),
#'   c(3, 3)
#' ))
#'
#' # throws an error
#' mat[3, 1] <- 5
#' try(rvine_matrix(mat))
#'
#' # C-vine structure
#' cvine <- cvine_structure(1:5)
#' cvine
#' plot(cvine)
#'
#' # D-vine structure
#' dvine <- dvine_structure(c(1, 4, 2, 3, 5))
#' dvine
#' plot(dvine)
#'
#' @name rvine_structure
#' @aliases rvine_matrix is.rvine_structure is.rvine_matrix dvine_structure
#'   cvine_structure
rvine_structure <- function(order, struct_array = list(), is_natural_order = FALSE) {

  # sanity checks and extract dimension/trunc_lvl
  assert_that(is.vector(order) && all(sapply(order, is.count)),
    msg = "Order should be a vector of positive integers."
  )
  assert_that(is.vector(struct_array) || is.list(struct_array))
  assert_that(all(rapply(struct_array, function(i) sapply(i, is.count))),
    msg = "All elements of struct_array should be positive integers."
  )
  d <- length(order)
  assert_that(length(struct_array) < length(order),
    msg = paste("The length of struct array is incompatible with order,",
                "it should be in {1, ..., length(order) - 1}")
  )

  trunc_lvl <- length(struct_array)
  lens <- sapply(seq_len(trunc_lvl), function(j) length(struct_array[[j]]))
  assert_that(all(lens == seq(d - 1, d - trunc_lvl)),
    msg = paste(
      "The number of elements in struct_array is incompatible with order.",
      "The first entry of struct_array should be a vectors of size" ,
      "length(order) - 1, the second of size length(order) - 2), etc.")
  )

  # create and check output
  output <- rvine_structure_cpp(
    list(
      order = order,
      struct_array = struct_array,
      d = d,
      trunc_lvl = trunc_lvl
    ),
    TRUE, is_natural_order
  )
  output <- structure(output, class = c("rvine_structure", "list"))

  # return output
  output
}

#' Plotting R-vine structures
#'
#' Plot one or all trees of an R-vine structure.
#'
#' @param x an `rvine_structure` or `rvine_matrix` object.
#' @param ... passed to `plot.vinecop_dist()`.
#' @aliases plot.rvine_matrix
#' @export
#' @examples
#' plot(cvine_structure(1:5))
#' plot(rvine_structure_sim(5))
plot.rvine_structure <- function(x, ...) {
  d <- dim(x)[1]
  trunc_lvl <- dim(x)[2]
  pcs <- lapply(seq_len(min(d - 1, trunc_lvl)),
                function(i) lapply(seq_len(d - i), function(j) bicop_dist()))
  plot.vinecop_dist(vinecop_dist(pcs, x), ...)
}

#' @rdname plot.rvine_structure
#' @examples
#' mat <- rbind(c(1, 1, 1), c(2, 2, 0), c(3, 0, 0))
#' plot(rvine_matrix(mat))
#' plot(rvine_matrix_sim(5))
#' @export
plot.rvine_matrix <- function(x, ...) {
  plot(as_rvine_structure(x), ...)
}

#' @rdname rvine_structure
#' @param trunc_lvl the truncation level
#' @export
cvine_structure <- function(order, trunc_lvl = Inf) {
  if (is.count(order))
    order <- seq_len(order)
  assert_that(is.vector(order) && all(sapply(order, is.count)),
              msg = "Order should be a vector of positive integers.")
  assert_that(is.scalar(trunc_lvl) & is.number(trunc_lvl))

  d <- length(order)
  trunc_lvl <- min(trunc_lvl, d - 1)
  vars <- rev(seq_len(d))
  struct_array <- lapply(seq(1, trunc_lvl), function(i) rep(vars[i], d - i))

  output <- rvine_structure_cpp(
    list(
      order = order,
      struct_array = struct_array,
      d = d,
      trunc_lvl = trunc_lvl
    ),
    FALSE, TRUE
  )
  structure(output, class = c("rvine_structure", "list"))
}

#' @rdname rvine_structure
#' @export
dvine_structure <- function(order, trunc_lvl = Inf) {
  if (is.count(order))
    order <- seq_len(order)
  assert_that(is.vector(order) && all(sapply(order, is.count)),
              msg = "Order should be a vector of positive integers.")
  assert_that(is.scalar(trunc_lvl) & is.number(trunc_lvl))

  d <- length(order)
  trunc_lvl <- min(trunc_lvl, d - 1)
  vars <- rev(seq_len(d))
  struct_array <- lapply(seq(1, trunc_lvl), function(i) seq(i + 1, d))

  output <- rvine_structure_cpp(
    list(
      order = order,
      struct_array = struct_array,
      d = d,
      trunc_lvl = trunc_lvl
    ),
    FALSE, TRUE
  )
  structure(output, class = c("rvine_structure", "list"))
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
    min(which(matrix[, 1] == 0)) - 1,
    d - 1
  )
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
  write.table(format(x, zero.print = zero.print),
    quote = FALSE,
    row.names = FALSE, col.names = FALSE
  )
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


#' Simulate R-vine structures
#'
#' Simulates from a uniform distribution over all R-vine structures on d
#' variables. `rvine_structure_sim()` returns an [rvine_structure()] object,
#' `rvine_matrix_sim()` an [rvine_matrix()].
#'
#' @aliases rvine_matrix_sim
#'
#' @param d the number of variables
#' @param natural_order boolean; whether the structures should be in natural
#'   order (counter-diagonal is `1:d`).
#'
#' @seealso [rvine_structure()], [rvine_matrix()],
#'    [plot.rvine_structure()], [plot.rvine_matrix()]
#' @export
#' @examples
#' rvine_structure_sim(10)
#'
#' rvine_structure_sim(10, natural_order = TRUE)  # counter-diagonal is 1:d
#'
#' rvine_matrix_sim(10)
rvine_structure_sim <- function(d, natural_order = FALSE) {
  assert_that(is.count(d), is.flag(natural_order))
  rvine_structure_sim_cpp(d, natural_order, get_seeds())
}

#' @rdname rvine_structure_sim
#' @export
rvine_matrix_sim <- function(d, natural_order = FALSE) {
  as_rvine_matrix(rvine_structure_sim(d, natural_order))
}
