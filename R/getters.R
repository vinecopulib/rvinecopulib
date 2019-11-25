#' Extracts components of `bicop_dist` and `vinecop_dist` objects
#'
#' Extracts either the structure matrix  (for `vinecop_dist` only), or
#' pair-copulas, their parameters, Kendall's taus, or families (for `bicop_dist`
#' and `vinecop_dist`).
#'
#' @name getters
#' @aliases get_pair_copula get_all_pair_copulas get_parameters
#'   get_all_parameters get_ktau get_all_ktaus get_matrix get_structure
#' @param object a `bicop_dist`, `vinecop_dist` or `vine_dist` object.
#' @details #' The [get_structure] method (for `vinecop_dist` or `vine_dist`
#' objects only) extracts the structure (see [rvine_structure] for more
#' details).
#'
#' The [get_matrix] method (for `vinecop_dist` or `vine_dist` objects only)
#' extracts the structure matrix (see [rvine_structure] for more details).
#'
#' The other `get_xyz` methods for `vinecop_dist` or `vine_dist` objects return
#' the entries corresponding to the pair-copula indexed by its `tree` and
#' `edge`. When `object` is of class `bicop_dist`, `tree` and `edge` are not
#' required. \cr
#'
#' * [get_pair_copula()] = the pair-copula itself (see [bicop]).
#'
#' * [get_parameters()] = the parameters of the pair-copula (i.e., a `numeric`
#' scalar, vector, or matrix).
#'
#' * [get_family()] = a character for the family (see [bicop] for
#' implemented families).
#'
#' * [get_ktau()] = a `numeric` scalar with the pair-copula Kendall's tau.
#'
#' The `get_all_xyz` methods (for `vinecop_dist` or `vine_dist` objects only)
#' return lists of lists, with each element corresponding to a tree in `trees`,
#' and then elements of the sublists correspond to edges. The returned lists
#' have two additional attributes:
#'
#' * `"d"` = the dimension of the model.
#'
#' * `"trees"` = the extracted trees.
#'
#' @examples
#' # specify pair-copulas
#' bicop <- bicop_dist("bb1", 90, c(3, 2))
#' pcs <- list(
#'   list(bicop, bicop), # pair-copulas in first tree
#'   list(bicop) # pair-copulas in second tree
#' )
#'
#' # specify R-vine matrix
#' mat <- matrix(c(1, 2, 3, 1, 2, 0, 1, 0, 0), 3, 3)
#'
#' # set up vine copula model
#' vc <- vinecop_dist(pcs, mat)
#'
#' # get the structure
#' get_structure(vc)
#' all.equal(get_matrix(vc), mat, check.attributes = FALSE)
#'
#' # get pair-copulas
#' get_pair_copula(vc, 1, 1)
#' get_all_pair_copulas(vc)
#' all.equal(get_all_pair_copulas(vc), pcs, check.attributes = FALSE)
#' @return
#' The structure matrix, or pair-copulas, their parameters, Kendall's taus,
#' or families.
#' @rdname getters
#' @export
get_structure <- function(object) {
  assert_that(inherits(object, "vinecop_dist") ||
    inherits(object, "vine_dist"))
  if (inherits(object, "vinecop_dist")) {
    return(object$structure)
  } else {
    return(object$copula$structure)
  }
}

#' @export
get_matrix <- function(object) {
  as_rvine_matrix(get_structure(object))
}


#' @param tree tree index (not required if `object` is of class `bicop_dist`).
#' @param edge edge index (not required if `object` is of class `bicop_dist`).
#'
#' @rdname getters
#' @export
get_pair_copula <- function(object, tree = NA, edge = NA) {

  ## sanity checks
  assert_that(inherits(object, "bicop_dist") ||
    inherits(object, "vinecop_dist") ||
    inherits(object, "vine_dist"))

  if (inherits(object, "bicop_dist")) {
    if (!is.scalar(tree) || !is.na(tree)) {
      warning("tree argument not used for bicop_dist objects")
    }
    if (!is.scalar(edge) || !is.na(edge)) {
      warning("edge argument not used for bicop_dist objects")
    }
    return(object)
  } else {
    d <- dim(object)
    assert_that(is.number(tree),
      tree >= 1,
      tree <= d[2],
      msg = "tree should be a number between 1 and
                    the truncation level."
    )
    assert_that(is.number(edge),
      edge >= 1,
      edge <= d[1] - tree,
      msg = "tree should be a number between 1 and
                    dimension minus tree."
    )

    ## return pair-copula
    if (inherits(object, "vinecop_dist")) {
      return(object$pair_copulas[[tree]][[edge]])
    } else {
      return(object$copula$pair_copulas[[tree]][[edge]])
    }
  }
}

#' @rdname getters
#' @export
get_parameters <- function(object, tree = NA, edge = NA) {
  coef(get_pair_copula(object, tree, edge))
}

#' @rdname getters
#' @export
get_ktau <- function(object, tree = NA, edge = NA) {
  par_to_ktau(get_pair_copula(object, tree, edge))
}

#' @rdname getters
#' @export
get_family <- function(object, tree = NA, edge = NA) {
  get_pair_copula(object, tree, edge)$family
}

#' @param trees the trees to extract from `object` (`trees = NA` extracts all
#' trees).
#' @rdname getters
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

  if (!any(is.na(trees))) {
    assert_that(is.numeric(trees),
      all(trees >= 1),
      all(trees <= d[2]),
      msg = "the elements of trees should be numbers between 1 and
                    the truncation level."
    )
  }

  t <- length(pcs)
  if (any(is.na(trees))) {
    trees <- seq_len(t)
  } else {
    if (any(trees > t)) {
      warning(
        "vine copula is ", t, "-truncated; ",
        "only returning available trees."
      )
      trees <- trees[trees <= t]
    }
  }
  res <- structure(pcs[trees], class = c("rvine_list", "pair_copulas"))
  attr(res, "d") <- d[1]
  attr(res, "trunc_lvl") <- d[2]
  attr(res, "trees") <- trees
  return(res)
}

get_all_xyz <- function(object, trees, func, list_name) {
  pcs <- get_all_pair_copulas(object, trees)
  res <- structure(lapply(pcs, function(tree) lapply(tree, func)),
    class = c("rvine_list", list_name)
  )
  attr(res, "d") <- attr(pcs, "d")
  attr(res, "trees") <- attr(pcs, "trees")
  return(res)
}

#' @rdname getters
#' @export
get_all_parameters <- function(object, trees = NA) {
  get_all_xyz(object, trees, coef, "parameters")
}

#' @rdname getters
#' @export
get_all_ktaus <- function(object, trees = NA) {
  get_all_xyz(object, trees, par_to_ktau, "ktaus")
}

#' @rdname getters
#' @export
get_all_families <- function(object, trees = NA) {
  get_all_xyz(object, trees, function(pc) pc$family, "families")
}

#' @export
print.rvine_list <- function(x, ...) {
  what <- switch(class(x)[2],
    pair_copulas = "pair-copulas",
    parameters = "copula parameters",
    ktaus = "Kendall's taus",
    families = "copula families"
  )

  d <- attr(x, "d")
  trees <- attr(x, "trees")
  msg <- paste(
    "Nested list of lists for the", what, "of a",
    d, "dimensional vine with"
  )

  ntrees <- length(trees)
  if (ntrees == d - 1) {
    cat(paste(msg, "all trees: \n"))
  } else if (ntrees == 1) {
    cat(paste(msg, " tree ", trees, ". \n", sep = ""))
  } else {
    trees_print <- paste(
      paste(trees[1:(ntrees - 1)], collapse = ", "),
      "and", trees[ntrees]
    )
    cat(paste(msg, " trees ", trees_print, ". \n", sep = ""))
  }
  cat(rep("-", 50), "\n")
  arg <- deparse(substitute(x))
  for (t in seq_along(trees)) {
    if (t > 10) {
      if (length(trees) - 10 > 1) {
        trees_print <- "more trees.\n"
      } else {
        trees_print <- "more tree.\n"
      }
      cat("# ... with", length(trees) - 10, trees_print)
      break
    }
    cat(paste(arg, "[[", t, "]] -> a list with the ", d - trees[t],
      " ", what, " of tree ", trees[t], ". \n",
      sep = ""
    ))
  }
  invisible(x)
}
