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
#' @aliases contour.vinecop_dist plot.vinecop contour.vinecop
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
#' pcs <- list(
#'     lapply(runif(3),  # pair-copulas in first tree
#'           function(cor) bicop_dist("gaussian", 0, cor)),
#'     lapply(runif(2),  # pair-copulas in second tree
#'           function(cor) bicop_dist("gaussian", 0, cor)),
#'     lapply(runif(1),  # pair-copulas in first tree
#'           function(cor) bicop_dist("gaussian", 0, cor))
#' )
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
    
    M <- x$matrix
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

#' @method contour vinecop_dist
#' @rdname plot.vinecop_dist
#' @importFrom graphics par plot.new plot.window abline polygon
#' @importFrom graphics strheight strwidth text
#' @export
contour.vinecop_dist <- function(x, tree = "ALL", cex.nums = 1, data = NULL, ...) {
    
    ## check input
    d <- nrow(x$matrix)
    if (!inherits(x, "vinecop_dist"))
        stop("'x' has to be an vinecop_dist object.")
    if (tree != "ALL" && any(tree > d - 1))
        stop("Selected tree does not exist.")
    if (any(tree == "ALL")) {
        tree <- 1:(d - 1)
    }
    n.tree <- length(tree)
    
    if (!is.null(list(...)$type))
        stop("Only contour plots allowed. Don't use the type argument!")
    if (!is.null(list(...)$margins)) {
        if (!(margins %in% c("unif", "norm", "exp", "flexp")))
            stop("'margins' has to be one of 'unif', 'norm', 'exp', or 'flexp'.")
    } else {
        margins <- "norm"
    }
    if (is.null(list(...)$xlim) & is.null(list(...)$ylim)) {
        xylim <- switch(margins,
                        "unif"  = c(1e-2, 1 - 1e-2),
                        "norm"  = c(-3, 3),
                        "exp"   = c(0, 6),
                        "flexp" = c(-6, 0))
        xlim <- ylim <- xylim
    } else {
        xlim <- list(...)$xlim
        ylim <- list(...)$ylim
        xylim <- range(c(list(...)$xlim, list(...)$ylim))
    }
    
    ## set up for plotting windows (restore settings on exit)
    usr <- par(mfrow = c(n.tree, d - min(tree)), mar = rep(0, 4))
    on.exit(par(usr))
    
    ## calculate pseudo-observations (if necessary)
    #psobs <- if (!is.null(data)) vine_psobs(data, x) else NULL
    psobs <- NULL
    
    # contours: adjust limits for headings
    offs <- 0.25
    mult <- 1.35
    ylim[2] <- ylim[2] + offs*diff(ylim)
    
    
    ## run through trees -----------------------------------------------
    # initialize check variables
    cnt <- 0
    k <- d
    e <- numeric(0)
    class(e) <- "try-error"
    
    while ("try-error" %in% class(e)) {
        e <- try({
            maxnums <- get_name(1, max(tree), x)
            for (i in rev(tree)) {
                for (j in 1:(d - min(tree))) {
                    if (j <= d - i) {
                        if (is.null(psobs)) {
                            pcfit <- x$pair_copulas[[i]][[j]]
                        } else {
                            pcfit <- bicop(psobs[[i]][[j]], "tll")
                        }
                        
                        # set up list of contour arguments
                        args <- list(x = pcfit,
                                     drawlabels = FALSE,
                                     xlab = "",
                                     ylab = "",
                                     xlim = xlim,
                                     ylim = ylim,
                                     xaxt = "n",
                                     yaxt = "n",
                                     add  = TRUE)
                        
                        # create empty plot
                        plot.new()
                        plot.window(xlim = xlim, ylim = ylim,
                                    xaxs = "i",  yaxs = "i")
                        
                        # call contour.bicop with ... arguments
                        do.call(contour, modifyList(args, list(...)))
                        
                        # draw area for headings
                        abline(h = ylim[2] - diff(ylim)/mult*offs)
                        polygon(x = c(xlim[1] - diff(xlim),
                                      xlim[1] - diff(xlim),
                                      xlim[2] + diff(xlim),
                                      xlim[2] + diff(xlim)),
                                y = c(ylim[2] + diff(ylim)/mult*offs,
                                      ylim[2] - diff(ylim)/mult*offs,
                                      ylim[2] - diff(ylim)/mult*offs,
                                      ylim[2] + diff(ylim)/mult*offs),
                                col = "grey")
                        
                        # add separating lines
                        abline(v = xlim)
                        abline(h = ylim)
                        
                        # add pair-copula ID
                        cx1 <- 0.75 * diff(xlim) / strwidth(maxnums)
                        ty <- ylim[2] - diff(ylim)/mult*offs
                        cx2 <- 0.75 * (ylim[2] - ty) / strheight(maxnums)
                        cx <- min(cx1, cx2)
                        text(x = sum(xlim)/2,
                             y = ty + 0.225 / cex.nums * (ylim[2] - ty),
                             cex    = cex.nums * cx,
                             labels = get_name(j, i, x),
                             pos    = 3,
                             offset = 0)
                    } else {
                        plot.new()
                    }
                }
            }
        }
        , silent = TRUE)
        
        ## adjust to figure margins if necessary
        if (length(tree) < 1)
            stop("Error in plot.new() : figure margins too large")
        if ("try-error" %in% class(e)) {
            cnt <- cnt + 1
            tree <- tree[-which(tree == max(tree))]
            par(mfrow = c(n.tree - cnt, d - min(tree)))
        }
    }
    
    ## message for the user if not all trees could be plotted -----------
    if (length(tree) != n.tree) {
        nmbr.msg <- as.character(tree[1])
        if (length(tree) > 2) {
            for (i in tree[-c(1, length(tree))]) {
                nmbr.msg <- paste(nmbr.msg, i, sep=", ")
            }
        }
        if (length(tree) > 1) {
            s.msg <- "s "
            nmbr.msg <- paste(nmbr.msg,
                              "and",
                              tree[length(tree)],
                              "were plotted. ")
        } else {
            s.msg <- " "
            nmbr.msg <- paste(nmbr.msg, "was plotted. ", sep=" ")
        }
        msg.space <- "There is not enough space."
        msg.tree <- paste("Only Tree",
                          s.msg,
                          nmbr.msg,
                          "Use the 'tree' argument or enlarge figure margins",
                          " to see the others.",
                          sep = "")
        message(paste(msg.space, msg.tree))
    }
}

# vine_psobs <- function (uev, object) {
#     uev <- as.matrix(uev)
#     if (ncol(uev) == 1)
#         uev <- matrix(uev, 1, nrow(uev))
#     if (any(uev > 1) || any(uev < 0))
#         stop("Data has be in the interval [0,1].")
#     n <- ncol(uev)
#     N <- nrow(uev)
#     if (ncol(uev) != ncol(object$matrix))
#         stop("Dimensions of 'data' and 'object' do not match.")
#     if (!is(object,"vinecop"))
#         stop("'object' has to be an vinecop object")
#     
#     o <- diag(object$Matrix)
#     oldobject <- object
#     if (any(o != length(o):1)) {
#         object <- normalizevinecop(object)
#         uev <- matrix(uev[, o[length(o):1]], N, n)
#     }
#     
#     ## initialize objects
#     CondDistr <- neededCondDistr(object$Matrix)
#     val <- array(1, dim = c(n, n, N))
#     out <- lapply(1:(n - 1), list)
#     V <- list()
#     V$direct <- array(NA, dim = c(n, n, N))
#     V$indirect <- array(NA, dim = c(n, n, N))
#     V$direct[n, , ] <- t(uev[, n:1])
#     
#     for (i in (n - 1):1) {
#         for (k in n:(i + 1)) {
#             ## extract data for current tree
#             m <- object$MaxMat[k, i]
#             zr1 <- V$direct[k, i, ]
#             if (m == object$Matrix[k, i]) {
#                 zr2 <- V$direct[k, (n - m + 1), ]
#             } else {
#                 zr2 <- V$indirect[k, (n - m + 1), ]
#             }
#             
#             ## store data
#             out[[n - k + 1]][[i]] <- cbind(zr2, zr1)
#             
#             ## calculate pseudo-observations for next tree
#             if (CondDistr$direct[k - 1, i])
#                 V$direct[k - 1, i, ] <- bicopHfunc1(zr2, zr1,
#                                                     object$family[k, i],
#                                                     object$par[k, i],
#                                                     object$par2[k, i],
#                                                     check.pars = FALSE)
#             if (CondDistr$indirect[k - 1, i])
#                 V$indirect[k - 1, i, ] <- bicopHfunc2(zr2, zr1,
#                                                       object$family[k, i],
#                                                       object$par[k, i],
#                                                       object$par2[k, i],
#                                                       check.pars = FALSE)
#         }
#     }
#     
#     ## return list of pseudo-observations
#     out
# }

