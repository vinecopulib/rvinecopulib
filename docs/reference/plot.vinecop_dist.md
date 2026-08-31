# Plotting `vinecop_dist` and `vinecop` objects.

There are two plotting generics for `vinecop_dist` objects.
`plot.vinecop_dist` plots one or all trees of a given R-vine copula
model. Edges can be labeled with information about the corresponding
pair-copula. `contour.vinecop_dist` produces a matrix of contour plots
(using `plot.bicop`).

## Usage

``` r
# S3 method for class 'vinecop_dist'
plot(x, tree = 1, var_names = "ignore", edge_labels = NULL, ...)

# S3 method for class 'vinecop'
plot(x, tree = 1, var_names = "ignore", edge_labels = NULL, ...)

# S3 method for class 'vinecop_dist'
contour(x, tree = "ALL", cex.nums = 1, ...)

# S3 method for class 'vinecop'
contour(x, tree = "ALL", cex.nums = 1, ...)
```

## Arguments

  - x:
    
    `vinecop_dist` object.

  - tree:
    
    `"ALL"` or integer vector; specifies which trees are plotted.

  - var\_names:
    
    integer; specifies how to make use of variable names:
    
      - \`"ignore"“ = variable names are ignored,
    
      - \`"use"“ = variable names are used to annotate vertices,
    
      - \`"legend"“ = uses numbers in plot and adds a legend for
        variable names,
    
      - \`"hide"“ = no numbers or names, just the node.

  - edge\_labels:
    
    character; options are:
    
      - `"family"` = pair-copula family (see `[bicop_dist()]`),
    
      - \`"tau"“ = pair-copula Kendall's tau
    
      - \`"family\_tau"“ = pair-copula family and Kendall's tau,
    
      - \`"pair"“ = the name of the involved variables.

  - ...:
    
    Unused for `plot` and passed to `contour.bicop` for `contour`.

  - cex.nums:
    
    numeric; expansion factor for font of the numbers.

## Details

If you want the contour boxes to be perfect squares, the plot height
should be `1.25/length(tree)*(d - min(tree))` times the plot width.

The `plot()` method returns an object that (among other things) contains
the `igraph` representation of the graph; see *Examples*.

## See also

`vinecop_dist`, `plot.bicop`

## Author

Thomas Nagler, Thibault Vatter

## Examples

``` r
# set up vine copula model
u <- matrix(runif(20 * 10), 20, 10)
vc <- vinecop(u, family_set = "indep")

# plot
plot(vc, tree = c(1, 2))

plot(vc, edge_labels = "pair")


# extract igraph representation
plt <- plot(vc, edge_labels = "family_tau")

igr_obj <- get("g", plt$plot_env)[[1]]
igr_obj  # print object
#> IGRAPH 116f02d UN-- 10 9 -- 
#> + attr: name (v/c), name (e/c)
#> + edges from 116f02d (vertex names):
#> [1] 7 --10 10--8  4 --2  9 --8  3 --8  2 --8  8 --1  5 --6  6 --1 
igraph::E(igr_obj)$name  # extract edge labels
#>       [,1]      
#>  [1,] "indep(0)"
#>  [2,] "indep(0)"
#>  [3,] "indep(0)"
#>  [4,] "indep(0)"
#>  [5,] "indep(0)"
#>  [6,] "indep(0)"
#>  [7,] "indep(0)"
#>  [8,] "indep(0)"
#>  [9,] "indep(0)"

# set up another vine copula model
pcs <- lapply(1:3, function(j) # pair-copulas in tree j
  lapply(runif(4 - j), function(cor) bicop_dist("gaussian", 0, cor)))
mat <- rvine_matrix_sim(4)
vc <- vinecop_dist(pcs, mat)

# contour plot
contour(vc)
```
