#' Coerce various kind of objects to R-vine structures and matrices
#'
#' `as_rvine_structure` and `as_rvine_matrix` are new S3 generics allowing to
#' coerce objects into R-vine structures and matrices (see [rvine_structure()]
#' and [rvine_matrix()]).
#'
#' The coercion to `rvine_structure` and `rvine_matrix` can be applied to
#' different kind of objects Currently, `rvine_structure`, `rvine_matrix`,
#' `matrix` and `list` are supported.
#'
#' For `as_rvine_structure`:
#'
#'  * `rvine_structure` : the main use case is to re-check an object via
#'  `validate = TRUE`.
#'  * `rvine_matrix` and `matrix` : allow to coerce matrices into R-vine
#'  structures (see [rvine_structure()] for more details). The main difference
#'  between `rvine_matrix` and `matrix` is the nature of the validity
#'  checks.
#'  * `list` : must contain named elements `order` and `struct_array` to be
#'  coerced into an R-vine structure (see [rvine_structure()] for more details).
#'
#' #' For `as_rvine_matrix`:
#'
#'  * `rvine_structure` : allow to coerce an `rvine_structure` into an
#'  R-vine matrix (useful e.g. for printing).
#'  * `rvine_matrix`: similar to `as_rvine_structure` for `rvine_structure`,
#'  the main use case is to re-check an object via `validate = TRUE`.
#'  * `matrix` : allow to coerce matrices into R-vine
#'  matrices (mainly by checking that the matrix defines a valid
#'  R-vine, see [rvine_matrix()] for more details).
#'  * `list` : must contain named elements `order` and `struct_array` to be
#'  coerced into an R-vine matrix (see [rvine_structure()] for more details).
#'
#' @param x An object of class `rvine_structure`, `rvine_matrix`, `matrix` or
#' `list` that can be coerced into an R-vine structure or R-vine matrix
#' (see *Details*).
#' @param ... Other arguments passed on to individual methods.
#' @return Either an object of class `rvine_structure` or of class
#' `rvine_matrix` (see [rvine_structure()] or [rvine_matrix()]).
#' @examples
#' # R-vine structures can be constructed from the order vector and struct_array
#' rvine_structure(order = 1:4, struct_array = list(
#'   c(4, 4, 4),
#'   c(3, 3),
#'   2
#' ))
#' 
#' # ... or a similar list can be coerced into an R-vine structure
#' as_rvine_structure(list(order = 1:4, struct_array = list(
#'   c(4, 4, 4),
#'   c(3, 3),
#'   2
#' )))
#' 
#' # similarly, standard matrices can be coerced into R-vine structures
#' mat <- matrix(c(4, 3, 2, 1, 4, 3, 2, 0, 4, 3, 0, 0, 4, 0, 0, 0), 4, 4)
#' as_rvine_structure(mat)
#' 
#' # or truncate and construct the structure
#' mat[3, 1] <- 0
#' as_rvine_structure(mat)
#' 
#' # throws an error
#' mat[3, 1] <- 5
#' try(as_rvine_structure(mat))
#' @seealso rvine_structure rvine_matrix
#' @rdname as_rvine_structure
#' @aliases as_rvine_matrix
#' @export
as_rvine_structure <- function(x, ...) {
  UseMethod("as_rvine_structure")
}

#' @export
#' @rdname as_rvine_structure
as_rvine_matrix <- function(x, ...) {
  UseMethod("as_rvine_matrix")
}


#' @param validate When `TRUE``, verifies that the input is a valid
#' rvine-structure (see *Details*). You may want to suppress this when you
#' know that you already have a valid structure and you want to save some time,
#' or to explicitly enable it if you have a structure that you want to re-check.
#' @export
#' @rdname as_rvine_structure
as_rvine_structure.rvine_structure <- function(x, ..., validate = FALSE) {
  assert_that(is.flag(validate))
  if (validate) {
    validate_rvine_structure(x)
  } else {
    assert_that(is.rvine_structure(x))
  }

  x
}

#' @export
#' @rdname as_rvine_structure
as_rvine_matrix.rvine_structure <- function(x, ..., validate = FALSE) {
  assert_that(is.flag(validate))
  if (validate) {
    validate_rvine_structure(x)
  } else {
    assert_that(is.rvine_structure(x))
  }

  # extract order and dimension
  order <- x$order
  d <- dim(x)[1]

  # set-up output
  matrix <- matrix(0, d, d)

  # fill output
  diag(matrix[d:1, ]) <- order
  for (i in 1:(d - 1)) {
    newcol <- order[x[["struct_array"]][[i]]]
    matrix[1:length(newcol), i] <- newcol
  }

  class(matrix) <- c("rvine_matrix", class(matrix))
  attr(matrix, "d") <- d
  attr(matrix, "trunc_lvl") <- dim(x)[2]
  matrix
}

#' @param is_natural_order A flag indicating whether the `struct_array` element
#' of `x` is assumed to be provided in natural order already (a structure is in 
#' natural order if the anti-diagonal is 1, .., d from bottom left to top 
#' right).
#' @param byrow whether the element of the list named `struct_array` 
#' is assumed to be provided by column or by row.
#' @export
#' @rdname as_rvine_structure
as_rvine_structure.list <- function(x, ..., is_natural_order = FALSE, byrow = TRUE) {
  assert_that(
    is.list(x),
    all(c("order", "struct_array") %in% names(x)),
    is.flag(is_natural_order),
    is.flag(byrow)
  )

  rvine_structure(
    x[["order"]],
    x[["struct_array"]],
    is_natural_order,
    byrow
  )
}

#' @export
#' @rdname as_rvine_structure
as_rvine_matrix.list <- function(x, ..., is_natural_order = FALSE, byrow = TRUE) {
  assert_that(
    is.list(x),
    all(c("order", "struct_array") %in% names(x)),
    is.flag(is_natural_order),
    is.flag(byrow)
  )

  rvine_struct <- rvine_structure(
    x[["order"]],
    x[["struct_array"]],
    is_natural_order,
    byrow
  )

  as_rvine_matrix.rvine_structure(rvine_struct)
}


#' @export
#' @rdname as_rvine_structure
as_rvine_structure.rvine_matrix <- function(x, ..., validate = FALSE) {
  assert_that(is.flag(validate))
  if (validate) {
    validate_rvine_matrix(x)
  } else {
    assert_that(is.rvine_matrix(x))
  }

  # compute structure array in natural order
  d <- dim(x)[1]
  order <- order(diag(x[d:1, ]))
  struct_array <- lapply(1:(d - 1), function(i) order[x[1:(d - i), i]])

  # create and return x
  structure(list(
    order = diag(x[d:1, ]),
    struct_array = struct_array,
    d = d,
    trunc_lvl = dim(x)[2]
  ),
  class = c("rvine_structure", "list")
  )
}

#' @export
#' @rdname as_rvine_structure
as_rvine_matrix.rvine_matrix <- function(x, ..., validate = FALSE) {
  assert_that(is.flag(validate))
  if (validate) {
    validate_rvine_matrix(x)
  } else {
    assert_that(is.rvine_matrix(x))
  }

  x
}


#' @export
#' @rdname as_rvine_structure
as_rvine_structure.matrix <- function(x, ..., validate = TRUE) {
  x <- as_rvine_matrix(x, validate)
  as_rvine_structure.rvine_matrix(x)
}

#' @export
#' @rdname as_rvine_structure
as_rvine_matrix.matrix <- function(x, ..., validate = TRUE) {
  assert_that(is.flag(validate))
  if (validate) {
    x <- rvine_matrix(x)
  } else {
    assert_that(is.matrix(x) && is.numeric(x))
    x <- rvine_matrix_nocheck(x)
  }

  x
}
