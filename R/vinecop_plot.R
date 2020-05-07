#' Plotting \code{vinecop_dist} and `vinecop` objects.
#'
#' There are two plotting generics for \code{vinecop_dist} objects.
#' \code{plot.vinecop_dist} plots one or all trees of a given R-vine copula
#' model. Edges can be labeled with information about the corresponding
#' pair-copula. \code{contour.vinecop_dist} produces a matrix of contour plots
#' (using \code{\link[rvinecopulib:plot.bicop]{plot.bicop}}).
#'
#' If you want the contour boxes to be perfect squares, the plot height should
#' be \code{1.25/length(tree)*(d - min(tree))} times the plot width.
#'
#' @method plot vinecop_dist
#'
#' @aliases contour.vinecop_dist plot.vinecop contour.vinecop
#'
#' @param x \code{vinecop_dist} object.
#' @param tree \code{"ALL"} or integer vector; specifies which trees are
#' plotted.
#' @param var_names integer; specifies how to make use of variable names:
#' * `"ignore"`` = variable names are ignored,
#' * `"use"`` = variable names are used to annotate vertices,
#' * `"legend"`` = uses numbers in plot and adds a legend for variable names,
#' * `"hide"`` = no numbers or names, just the node.
#' @param edge_labels character; options are:
#' * `"family"` = pair-copula family (see `[bicop_dist()]`),
#' * `"tau"`` = pair-copula Kendall's tau
#' * `"family_tau"`` = pair-copula family and Kendall's tau,
#' * `"pair"`` = the name of the involved variables.
#' @param cex.nums numeric; expansion factor for font of the numbers.
#' @param \dots Unused for \code{plot} and passed to
#' \code{\link[rvinecopulib:contour.bicop]{contour.bicop}} for \code{contour}.
#'
#' @author Thomas Nagler, Thibault Vatter
#'
#' @seealso \code{\link[rvinecopulib:vinecop_dist]{vinecop_dist}},
#' \code{\link[rvinecopulib:plot.bicop]{plot.bicop}}
#'
#' @keywords plot
#'
#' @examples
#' # set up vine copula model
#' d <- 20
#' n <- 2e2
#' u <- matrix(runif(n * d), n, d)
#' vc <- vinecop(u, family = "indep")
#'
#' # plot
#' plot(vc, tree = c(1, 2))
#' plot(vc, edge_labels = "pair")
#'
#' # set up another vine copula model
#' pcs <- lapply(1:3, function(j) # pair-copulas in tree j
#'   lapply(runif(4 - j), function(cor) bicop_dist("gaussian", 0, cor)))
#' mat <- rvine_matrix_sim(4)
#' vc <- vinecop_dist(pcs, mat)
#'
#' # contour plot
#' contour(vc)
#' @export
plot.vinecop_dist <- function(x, tree = 1, var_names = "ignore",
                              edge_labels = NULL, ...) {
  if (!requireNamespace("ggraph", quietly = TRUE)) {
    stop("The 'ggraph' package must be installed to plot.")
  }
  if (!requireNamespace("grid", quietly = TRUE)) {
    stop("The 'grid' package must be installed to plot.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The 'ggplot2' package must be installed to plot.")
  }
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("The 'igraph' package must be installed to plot.")
  }

  d <- dim(x)[1]
  trunc_lvl <- dim(x)[2]

  ## sanity checks
  if (!identical(tree, "ALL") && any(tree > trunc_lvl)) {
    stop("Selected tree does not exist.")
  }
  if (any(tree == "ALL")) {
    if (d > 5) {
      warning(paste(
        "tree = 'ALL' is not recommended for d > 5 and",
        " it is set as c(1,2), please use tree = 1:d"
      ))
      tree <- c(1, 2)
    } else {
      tree <- seq_len(trunc_lvl)
    }
  }
  assert_that(in_set(var_names, c("ignore", "use", "legend", "hide")))
  if (!is.null(edge_labels)) {
    assert_that(in_set(edge_labels, c("pair", "tau", "family", "family_tau")))
  }

  ## set names if empty
  if (is.null(x$names)) {
    x$names <- as.character(1:d)
  }
  if (var_names %in% c("ignore", "legend")) {
    names <- x$names
    x$names <- as.character(1:d)
  }

  #### loop through the trees and create graph objects
  g <- lapply(tree, get_graph, vc = x, edge_labels = edge_labels)

  plots <- vector("list", length(tree))
  name <- NULL # for the CRAN check
  for (i in seq_along(tree)) {
    p <- ggraph::ggraph(g[[i]], "igraph",
      algorithm = "tree", circular = TRUE
    )
    if (!is.null(edge_labels)) {
      p <- p +
        ggraph::geom_edge_link(ggplot2::aes(label = name),
          colour = "#000000",
          angle_calc = "along",
          label_dodge = grid::unit(7, "points")
        )
    } else {
      p <- p + ggraph::geom_edge_link(colour = "#000000")
    }
    p <- p +
      ggraph::geom_node_point(col = "#56B4E9", size = 3) +
      ggplot2::theme_void() +
      ggplot2::labs(title = paste0("Tree ", tree[i]))
    if (var_names != "hide") {
      p <- p + ggraph::geom_node_text(ggplot2::aes(label = name),
        fontface = "bold",
        repel = TRUE
      )
    }
    if (var_names == "legend") {
      p <- p + ggplot2::labs(caption = paste(x$names, names,
        sep = " = ",
        collapse = ", "
      ))
    }
    plots[[i]] <- p +
      ggplot2::theme(plot.margin = ggplot2::margin(5, 5, 5, 5, "pt"))
  }

  if (length(tree) > 3) {
    invisible(multiplot(plotlist = plots, cols = 2))
  } else {
    return(invisible(multiplot(plotlist = plots)))
  }
}

#' @rdname plot.vinecop_dist
#' @export
plot.vinecop <- plot.vinecop_dist

## creates a graph object for a tree in a given vinecop_dist
get_graph <- function(tree, vc, edge_labels) {
  M <- get_matrix(vc)
  d <- dim(vc)[1]

  I <- matrix(0, d - tree + 1, d - tree + 1)

  ## extract node and edge labels as numbers
  if (tree > 1) {
    V <- t(sapply(seq.int(d - tree + 1), function(j)
      M[c(d - j + 1, (tree - 1):1), j]))
  } else {
    V <- matrix(diag(M[d:1, ]), ncol = 1)
  }
  E <- t(sapply(seq.int(d - tree), function(j)
    M[c(d - j + 1, tree:1), j]))

  ## build adjacency matrix by matching V and E
  for (i in 1:nrow(E)) {
    ind.i <- which(apply(V, 1, function(x) all(x %in% E[i, ])))
    I[ind.i[1], ind.i[2]] <- I[ind.i[2], ind.i[1]] <- 1
  }

  ## convert to variable names
  if (tree > 1) {
    colnames(I) <- rownames(I) <- sapply(1:(d - tree + 1),
      get_name,
      tree = tree - 1,
      vc = vc
    )
  } else {
    colnames(I) <- rownames(I) <- vc$names[diag(M[d:1, ])]
  }

  ## build graph
  g <- igraph::graph_from_adjacency_matrix(I, mode = "undirected")

  ## add edge labels
  if (!is.null(edge_labels)) {
    igraph::E(g)$name <- sapply(tree, set_edge_labels,
      vc = vc,
      edge_labels = edge_labels
    )
  }

  g
}

## finds appropriate edge labels for the plot
set_edge_labels <- function(tree, vc, edge_labels) {
  d <- dim(vc)[1]
  get_edge_label <- switch(edge_labels,
    family = get_family,
    tau = function(vc, tree, edge)
      round(get_ktau(vc, tree, edge),
        digits = 2
      ),
    family_tau = get_family_tau,
    pair = get_name
  )
  sapply(1:(d - tree), function(edge) get_edge_label(vc, tree, edge))
}


## get info for a pair-copula
get_name <- function(vc, tree, edge) {
  M <- get_matrix(vc)
  d <- nrow(M)
  # conditioned set
  bef <- paste0(vc$names[M[c(d - edge + 1, tree), edge]],
    collapse = ","
  )
  # conditioning set
  aft <- ifelse(tree > 1,
    paste0(vc$names[M[(tree - 1):1, edge]], collapse = ","),
    ""
  )
  # paste together
  sep <- ifelse(tree > 1, ";", "")
  paste(bef, aft, sep = sep, collapse = "")
}

get_family_tau <- function(vc, tree, edge) {
  paste0(
    get_family(vc, tree, edge), "(",
    round(get_ktau(vc, tree, edge), digits = 2), ")"
  )
}

#' @method contour vinecop_dist
#' @rdname plot.vinecop_dist
#' @importFrom graphics par plot.new plot.window abline polygon
#' @importFrom graphics strheight strwidth text
#' @export
contour.vinecop_dist <- function(x, tree = "ALL", cex.nums = 1, ...) {
  ## check input
  d <- dim(x)[1]
  trunc_lvl <- dim(x)[2]
  assert_that(is.number(cex.nums))
  assert_that((length(tree) == 1 && tree == "ALL") || all(tree <= trunc_lvl),
    msg = "Selected tree does not exist."
  )

  if (any(tree == "ALL")) {
    tree <- seq_len(trunc_lvl)
  }

  n_tree <- length(tree)

  if (!is.null(list(...)$var_names)) {
    stop("Only contour plots allowed. Don't use the var_names argument!")
  }
  if (!is.null(list(...)$margins)) {
    margins <- list(...)$margins
    assert_that(in_set(margins, c("unif", "norm", "exp", "flexp")))
  } else {
    margins <- "norm"
  }
  if (is.null(list(...)$xlim) & is.null(list(...)$ylim)) {
    xylim <- switch(margins,
      "unif" = c(1e-2, 1 - 1e-2),
      "norm" = c(-3, 3),
      "exp" = c(0, 6),
      "flexp" = c(-6, 0)
    )
    xlim <- ylim <- xylim
  } else {
    xlim <- list(...)$xlim
    ylim <- list(...)$ylim
    xylim <- range(c(list(...)$xlim, list(...)$ylim))
  }

  if (is.null(x$names)) {
    x$names <- as.character(1:d)
  }

  ## set up for plotting windows (restore settings on exit)
  usr <- par(mfrow = c(n_tree, d - min(tree)), mar = rep(0, 4))
  on.exit(par(usr))

  # contours: adjust limits for headings
  offs <- 0.25
  mult <- 1.35
  ylim[2] <- ylim[2] + offs * diff(ylim)

  ## run through trees -----------------------------------------------
  # initialize check variables
  cnt <- 0
  k <- d
  maxnums <- get_name(x, max(tree, length(x$pair_copulas)), 1)
  for (i in rev(tree)) {
    for (j in 1:(d - min(tree))) {
      if (j <= d - i) {
        if (length(x$pair_copulas) >= i) {
          pcfit <- x$pair_copulas[[i]][[j]]
        } else {
          pcfit <- bicop_dist()
        }

        # set up list of contour arguments
        args <- list(
          x = pcfit,
          drawlabels = FALSE,
          xlab = "",
          ylab = "",
          xlim = xlim,
          ylim = ylim,
          xaxt = "n",
          yaxt = "n",
          add = TRUE
        )

        # create empty plot
        plot.new()
        plot.window(
          xlim = xlim, ylim = ylim,
          xaxs = "i", yaxs = "i"
        )

        # call contour.bicop with ... arguments
        do.call(contour, modifyList(args, list(...)))

        # draw area for headings
        abline(h = ylim[2] - diff(ylim) / mult * offs)
        polygon(
          x = c(
            xlim[1] - diff(xlim),
            xlim[1] - diff(xlim),
            xlim[2] + diff(xlim),
            xlim[2] + diff(xlim)
          ),
          y = c(
            ylim[2] + diff(ylim) / mult * offs,
            ylim[2] - diff(ylim) / mult * offs,
            ylim[2] - diff(ylim) / mult * offs,
            ylim[2] + diff(ylim) / mult * offs
          ),
          col = "grey"
        )

        # add separating lines
        abline(v = xlim)
        abline(h = ylim)

        # add pair-copula ID
        cx1 <- 0.75 * diff(xlim) / strwidth(maxnums)
        ty <- ylim[2] - diff(ylim) / mult * offs
        cx2 <- 0.75 * (ylim[2] - ty) / strheight(maxnums)
        cx <- min(cx1, cx2)
        text(
          x = sum(xlim) / 2,
          y = ty + 0.225 / cex.nums * (ylim[2] - ty),
          cex = cex.nums * cx,
          labels = get_name(x, i, j),
          pos = 3,
          offset = 0
        )
      } else {
        plot.new()
      }
    }
  }
}

#' @rdname plot.vinecop_dist
#' @export
contour.vinecop <- contour.vinecop_dist
