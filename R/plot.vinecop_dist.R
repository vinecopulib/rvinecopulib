#' Plotting \code{vinecop_dist} objects.
#'
#' There are two plotting generics for \code{vinecop_dist} objects.
#' \code{plot.vinecop_dist} plots one or all trees of a given R-vine copula
#' model. Edges can be labeld with information about the corresponding
#' pair-copula. \code{contour.vinecop_dist} produces a matrix of contour plots
#' (using \code{\link[rvinecopulib:plot.bicop]{plot.bicop}}).
#'
#' If you want the contour boxes to be perfect squares, the plot height should
#' be \code{1.25/length(tree)*(d - min(tree))} times the plot width.
#'
#' @method plot vinecop_dist
#'
#' @aliases vinecop_dist vinecop_dist plot.vinecop contour.vinecop
#'
#' @inheritParams plot.bicop_dist
#' @param x \code{vinecop_dist} object.
#' @param tree \code{"ALL"} or integer vector; specifies which trees are
#' plotted.
#' @param type integer; specifies how to make use of variable names: \cr
#' \code{0} = variable names are ignored, \cr \code{1} = variable names are
#' used to annotate vertices, \cr \code{2} = uses numbers in plot and adds a
#' legend for variable names.
#' @param edge.labels character; either a vector of edge labels or one of the
#' following: \cr \code{"family"} = pair-copula family (see
#' \code{\link[rvinecopulib:bicop_dist]{bicop_dist}}), \cr \code{"tau"} = 
#' pair-copula Kendall's tau\cr \code{"family_tau"} = pair-copula family and 
#' Kendall's tau, \cr \code{"pair"} = for the name of the involved variables.
#' @param \dots Arguments passed to 
#' \code{\link[rvinecopulib:contour.bicop]{contour.bicop}} respectively.
#'
#' @author Thomas Nagler, Thibault Vatter
#'
#' @seealso \code{\link[rvinecopulib:vinecop_dist]{vinecop_dist}},
#' \code{\link[rvinecopulib:plot.bicop]{plot.bicop}}
#'
#' @keywords plot
#'
#' @examples
#' # specify pair-copulas
# pcs <- list(
#     lapply(runif(3),  # pair-copulas in first tree
#           function(cor) bicop_dist("gaussian", 0, cor)),
#     lapply(runif(2),  # pair-copulas in second tree
#           function(cor) bicop_dist("gaussian", 0, cor)),
#     lapply(runif(1),  # pair-copulas in first tree
#           function(cor) bicop_dist("gaussian", 0, cor))
# )
#' 
#' # specify R-vine matrix
#' mat <- matrix(c(1, 2, 3, 4, 1, 2, 3, 0, 1, 2, 0, 0, 1, 0, 0, 0), 4, 4)
#' 
#' # set up vine copula model
#' vc <- vinecop_dist(pcs, mat)
#' 
#'
#' # plot trees
#' \dontrun{plot(vc)}
#'
#' @export
plot.vinecop_dist <- function(x, tree = 1, type = 0, edge.labels = NULL,
                              ...) {
    if (!requireNamespace("ggraph", quietly = TRUE))
        stop("The 'ggraph' package must be installed to plot.")
    if (!requireNamespace("grid", quietly = TRUE))
        stop("The 'grid' package must be installed to plot.")
    if (!requireNamespace("ggplot2", quietly = TRUE))
        stop("The 'ggplot2' package must be installed to plot.")
    if (!requireNamespace("igraph", quietly = TRUE))
        stop("The 'igraph' package must be installed to plot.")
    
    M <- vc$matrix
    d <- nrow(M)
    
    ## sanity checks
    if (!inherits(x, "vinecop_dist"))
        stop("'x' has to be an vinecop_dist object.")
    if (tree != "ALL" && any(tree > d - 1))
        stop("Selected tree does not exist.")
    if (any(tree == "ALL")) {
        if (d > 5) {
            warning(paste("tree = 'ALL' is not recommended for d > 5 and",
                          " it is set as c(1,2), please use tree = 1:d"))
            tree <- c(1,2)
        } else {
            tree <- 1:(d - 1)
        }
    }
    if (!all(type %in% c(0, 1, 2)))
        stop("type not implemented")
    if (!(is.null(edge.labels) || 
          any(edge.labels %in% c("pair","tau","family","family_tau"))))
        stop("edge.labels not implemented")
    
    ## set names if empty
    if (is.null(x$names))
        x$names <- as.character(1:d)
    if (type %in% c(0, 2)) {
        names <- x$names
        x$names <- as.character(1:d)
    }
    
    #### loop through the trees and create graph objects
    g <- lapply(tree, get_graph, vc = x, edge.labels = edge.labels, type = type)
    
    plots <- vector("list", length(tree))
    for (i in seq_along(tree)) {
        p <- ggraph::ggraph(g[[i]], 'igraph', algorithm = 'tree')
        if (!is.null(edge.labels)) {
            p <- p + 
                ggraph::geom_edge_link(ggplot2::aes(label = name),
                                       angle_calc = 'along',
                                       label_dodge = grid::unit(5e-2, "native"))
        } else {
            p <- p + ggraph::geom_edge_link()
        }
        p <- p +
            ggraph::geom_node_label(ggplot2::aes(label = name), 
                                    bg = 'steelblue', col = 'white') + 
            ggplot2::theme_void() +
            ggplot2::labs(title = paste0("Tree ", tree[i]))
        if (type == 2)
            p <- p + ggplot2::labs(caption = paste(x$names, names, 
                                                   sep = " = ", 
                                                   collapse = ", "))
        plots[[i]] <- p + ggplot2::theme(plot.margin = ggplot2::margin(5, 5, 
                                                                       5, 5, 
                                                                       "pt"))
    }
    
    if (length(tree) > 3) {
        invisible(multiplot(plotlist = plots, cols = 2))
    } else {
        return(invisible(multiplot(plotlist = plots)))
    }
}


## creates a graph object for a tree in a given vinecop_dist
get_graph <- function(tree, vc, edge.labels, type) {
    M <- vc$matrix
    d <- ncol(M)
    
    I <- matrix(0, d - tree + 1, d - tree + 1)
    
    ## extract node and edge labels as numbers
    if (tree > 1) {
        V <- t(sapply(seq.int(d - tree +1), function(j) 
            M[c(d - j + 1, (tree-1):1), j]))
    } else {
        V <- matrix(diag(M[d:1,]), ncol = 1)
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
                                             vc = vc)
    } else {
        colnames(I) <- rownames(I) <- vc$names[diag(M[d:1,])]
    }
    
    ## build graph
    g <- igraph::graph_from_adjacency_matrix(I, mode = "undirected")
    
    ## add edge labels
    if (!is.null(edge.labels)) {
        igraph::E(g)$name <- sapply(tree, set_edge_labels,
                                    vc = vc,
                                    edge.labels = edge.labels,
                                    type = type)
    }
    
    g
}

## finds appropriate edge labels for the plot
set_edge_labels <- function(tree, vc, edge.labels, type) {
    d <- nrow(vc$matrix)
    get_edge_label <- switch(edge.labels,
                             family = get_family,
                             tau = get_tau,
                             family_tau = get_family_tau,
                             pair = get_name)
    sapply(1:(d - tree),
           get_edge_label,
           tree = tree,
           vc = vc)
}


## get info for a pair-copula
get_name <-  function(j, tree, vc) {
    M <- vc$matrix
    d <- nrow(M)
    # conditioned set
    bef <- paste0(vc$names[M[c(d - j + 1, tree), j]], 
                  collapse = ",")
    # conditioning set
    aft <- ifelse (tree > 1, paste0(vc$names[M[(tree-1):1, j]], 
                                    collapse = ","), "")
    # paste together
    sep <- ifelse(tree > 1, ";", "")
    paste(bef, aft, sep = sep, collapse = "")
}

get_family <- function(j, tree, vc) {
    vc$pair_copulas[[tree]][[j]]$family
}

get_tau <- function(j, tree, vc) {
    round(par_to_tau(vc$pair_copulas[[tree]][[j]]), digits = 2)
}

get_family_tau <- function(j, tree, vc) {
    paste0(get_family(j, tree, vc), "(", get_tau(j, tree, vc), ")")
}

