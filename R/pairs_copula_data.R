#' Exploratory pairs plot for copula data
#'
#' This function provides pair plots for copula data. It shows
#' bivariate contour plots on the lower panel, scatter plots and
#' correlations on the upper panel and histograms on the diagonal panel.
#'
#' @param data the data (must lie in the unit hypercube).
#' @param \dots other parameters passed to `pairs.default()`,
#' `contour.bicop()`, `points.default()`, `hist.default()`, or `bicop()`.
#'
#' @importFrom stats cor
#' @importFrom graphics box hist.default pairs.default points.default text.default
#' @export
#'
#' @examples
#' u <- replicate(3, runif(100))
#' pairs_copula_data(u)
pairs_copula_data <- function(data, ...) {
  assert_that(is.matrix(data) || is.data.frame(data))
  assert_that(all(data < 1) && all(data > 0))

  old_par <- par(xaxt = "n", yaxt = "n")
  on.exit(par(old_par))

  labels <- colnames(data)
  if (is.null(labels))
    labels <- paste0("var", seq_len(ncol(data)))
  args <- list(x = data,
               labels = labels,
               lower.panel = function(x, y) lp_pairs_copula_data(x, y, ...),
               diag.panel = function(x) dp_pairs_copula_data(x, ...),
               upper.panel = function(x, y)  up_pairs_copula_data(x, y, ...),
               label.pos = 0.85,
               cex.labels = 1,
               gap = 0)
  call_with_dots(pairs.default, args, ...)
}

call_with_dots <- function(fun, args, ..., except = "") {
  dots <- list(...)
  dots <- dots[names(dots) %in% setdiff(names(formals(fun)), except)]
  args <- modifyList(args, dots)
  do.call(fun, args)
}

## lower panel: empirical contour plot
lp_pairs_copula_data <- function(x, y, ...) {
  old_par <- par(usr = c(-3, 3, -3, 3), new = TRUE)
  on.exit(par(old_par))

  args <- list(data = cbind(x, y), family_set = "tll")
  cop <- call_with_dots(bicop, args, ...)

  args <- list(x = cop,
               size = 100,
               axes = FALSE,
               drawlabels = FALSE)
  call_with_dots(contour.bicop, args, ...)
}

## upper panel: scatter plot (copula data) and correlation
up_pairs_copula_data <- function(x, y, ...) {
  old_par <- par(usr = c(0, 1, 0, 1), new = TRUE)
  on.exit(par(old_par))

  args <- list(x = x, y = y, pch = ".", cex = 2, col = "gray40")
  call_with_dots(points.default, args, ...)

  cor <- call_with_dots(cor, list(x = x, y = y, method = "kendall"), ...)
  txt <- format(x = cor, digits = 2, nsmall = 2)[1]
  args <- list(x = 0.5, y = 0.5,
               labels = txt,
               cex = 1 + abs(cor) * 2)
  call_with_dots(text.default, args, ...)
}

## diagonal panel: histograms (copula data)
dp_pairs_copula_data <- function(x, ...) {
  old_par <- par(usr = c(0, 1, 0, 1.6), new = TRUE)
  on.exit(par(old_par))

  args <- list(x = x,
               freq = FALSE,
               add = TRUE,
               border = "white",
               main = "",
               col = "grey")
  call_with_dots(hist.default, args, ..., except = c("col", "freq"))
  box()
  abline(h = 1, col = "black", lty = 3)
}
